#!/usr/bin/env python3
"""
EEH H5 → Cloud Optimized GeoTIFF (COG)
========================================
Supports: TES (LSTE), STIC (ET), GPP, WUE

Workflow per file:
  1. Read GEO H5 for lat/lon/land_fraction
  2. Optionally read CLOUD H5 for cloud mask (TES only for filtering)
  3. kd_tree nearest-neighbor resampling onto a regular grid
  4. Write intermediate GTiff
  5. gdal_translate → COG (deflate, tiled)
  6. Cleanup intermediate files

Cloud mask convention (aligned with GPP/WASDI):
  cloud_mask = 1 → clear,  cloud_mask = 0 → cloudy

WUE note:
  WUEd is stored inside GPP H5 files (EEH2STIC_L3_GPP_*.h5).
  Use force_product="WUE" to extract only the WUEd band.

Usage (standalone):
    python tools/h5_to_cog.py --eeh-dir ./output/tes --geo-dir ./data/geo --out-dir ./output/cog/tes
    python tools/h5_to_cog.py --eeh-dir ./output/gpp --geo-dir ./data/geo --out-dir ./output/cog/wue --force-product WUE
"""

from __future__ import annotations

import gc
import os
import re
import shutil
import subprocess
import sys
import argparse
import logging
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

import h5py
import numpy as np
import rasterio
from rasterio.transform import from_origin
from rasterio.crs import CRS
from pyresample import geometry, kd_tree

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
)
log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
NODATA = -9999.0
RADIUS_INFLUENCE = 500   # metres
RES = 0.0006             # ~67 m at equator
TMP_BASE = "/tmp/eeh_cog"

# ---------------------------------------------------------------------------
# Band definitions per product
# ---------------------------------------------------------------------------
PRODUCT_BANDS = {
    "TES": [
        {"name": "LST",           "scale": 0.02,  "offset": 0.0,      "nodata": 0,     "units": "K",             "long_name": "Land Surface Temperature"},
        {"name": "Emis2",         "scale": 0.002, "offset": 0.489999, "nodata": 0,     "units": "n/a",           "long_name": "Band 2 Emissivity"},
        {"name": "Emis4",         "scale": 0.002, "offset": 0.489999, "nodata": 0,     "units": "n/a",           "long_name": "Band 4 Emissivity"},
        {"name": "Emis5",         "scale": 0.002, "offset": 0.489999, "nodata": 0,     "units": "n/a",           "long_name": "Band 5 Emissivity"},
        {"name": "BBE",           "scale": 0.002, "offset": 0.489999, "nodata": 0,     "units": "n/a",           "long_name": "Broad Band Emissivity"},
        {"name": "qa",            "scale": 1,     "offset": 0.0,      "nodata": -9999, "units": "N/A",           "long_name": "Quality Assurance"},
        {"name": "land_fraction", "scale": 1,     "offset": 0.0,      "nodata": -9999, "units": "dimensionless", "long_name": "Land fraction (from GEO)"},
        {"name": "cloud_mask",    "scale": 1,     "offset": 0.0,      "nodata": -9999, "units": "dimensionless", "long_name": "Cloud mask (1=clear 0=cloudy)"},
    ],
    "STIC": [
        {"name": "LE",  "scale": 1, "offset": 0, "nodata": -9999, "units": "W.m-2",    "long_name": "Latent heat flux"},
        {"name": "H",   "scale": 1, "offset": 0, "nodata": -9999, "units": "W.m-2",    "long_name": "Sensible heat flux"},
        {"name": "G",   "scale": 1, "offset": 0, "nodata": -9999, "units": "W.m-2",    "long_name": "Soil heat flux"},
        {"name": "Rn",  "scale": 1, "offset": 0, "nodata": -9999, "units": "W.m-2",    "long_name": "Net radiation"},
        {"name": "gah", "scale": 1, "offset": 0, "nodata": -9999, "units": "m.s-1",    "long_name": "Aerodynamic conductance"},
        {"name": "gsc", "scale": 1, "offset": 0, "nodata": -9999, "units": "m.s-1",    "long_name": "Surface conductance"},
        {"name": "Ms",  "scale": 1, "offset": 0, "nodata": -9999, "units": "unitless", "long_name": "Water stress"},
        {"name": "Mrz", "scale": 1, "offset": 0, "nodata": -9999, "units": "unitless", "long_name": "Water stress root zone"},
        {"name": "ETD", "scale": 1, "offset": 0, "nodata": -9999, "units": "mm.day-1", "long_name": "Daily ET"},
    ],
    "GPP": [
        {"name": "GPPd",    "scale": 1, "offset": 0, "nodata": -9999, "units": "gC.m-2.day-1", "long_name": "Daily gross primary productivity"},
        {"name": "GPPi",    "scale": 1, "offset": 0, "nodata": -9999, "units": "gC.m-2.s-1",   "long_name": "STIC-gs calibrated BL-LUE GPP"},
        {"name": "GPP_gs",  "scale": 1, "offset": 0, "nodata": -9999, "units": "gC.m-2.s-1",   "long_name": "instantaneous STIC-gs derived GPP"},
        {"name": "GPP_LUE", "scale": 1, "offset": 0, "nodata": -9999, "units": "gC.m-2.s-1",   "long_name": "instantaneous big-leaf model derived GPP"},
        {"name": "LUCC",    "scale": 1, "offset": 0, "nodata": 0,     "units": "unitless",      "long_name": "GLC_FCS30D land cover class"},
    ],
    "WUE": [
        {"name": "WUEd", "scale": 1, "offset": 0, "nodata": -9999, "units": "gC.m-2/mm.H2O", "long_name": "Daily water use efficiency"},
    ],
}

# ---------------------------------------------------------------------------
# Filename patterns
# ---------------------------------------------------------------------------
EEH_RE = re.compile(
    r"EEH2(?:TES_L2_LSTE|STIC_L3_ET|STIC_L3_GPP|STIC_L4_WUE)_(\d+_\d+_\d{8}T\d{6})_"
)
GEO_RE = re.compile(r"ECOv002_L1B_GEO_(\d+_\d+_\d{8}T\d{6})_")
CLOUD_RE = re.compile(r"ECOv002_L2_CLOUD_(\d+_\d+_\d{8}T\d{6})_")


def detect_product(filename: str) -> str | None:
    if "TES_L2_LSTE" in filename:
        return "TES"
    if "STIC_L3_ET" in filename:
        return "STIC"
    if "STIC_L4_WUE" in filename:
        return "WUE"
    if "STIC_L3_GPP" in filename:
        return "GPP"
    return None


def extract_key(filename: str, pattern) -> str | None:
    m = pattern.search(Path(filename).name)
    return m.group(1) if m else None


# ---------------------------------------------------------------------------
# Index builders (local only)
# ---------------------------------------------------------------------------

def build_geo_index(geo_dir: str) -> dict[str, str]:
    index = {}
    if not os.path.isdir(geo_dir):
        log.warning(f"GEO directory not found: {geo_dir}")
        return index
    for fname in os.listdir(geo_dir):
        if fname.endswith(".h5"):
            key = extract_key(fname, GEO_RE)
            if key:
                index[key] = os.path.join(geo_dir, fname)
    log.info(f"GEO index: {len(index)} files in {geo_dir}")
    return index


def build_cloud_index(cloud_dir: str) -> dict[str, str]:
    index = {}
    if not cloud_dir or not os.path.isdir(cloud_dir):
        return index
    for fname in os.listdir(cloud_dir):
        if fname.endswith(".h5"):
            key = extract_key(fname, CLOUD_RE)
            if key:
                index[key] = os.path.join(cloud_dir, fname)
    log.info(f"CLOUD index: {len(index)} files in {cloud_dir}")
    return index


# ---------------------------------------------------------------------------
# GTiff → COG via gdal_translate
# ---------------------------------------------------------------------------

def gtiff_to_cog(input_path: str, output_path: str) -> bool:
    cmd = [
        "gdal_translate", input_path, output_path,
        "--config", "GDAL_CACHEMAX", "512",
        "-of", "COG",
        "-co", "COMPRESS=DEFLATE",
        "-co", "BLOCKSIZE=512",
        "-co", "INTERLEAVE=BAND",
        "-co", "NUM_THREADS=ALL_CPUS",
        "-co", "BIGTIFF=IF_SAFER",
        "-co", "OVERVIEW_RESAMPLING=NEAREST",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        log.warning(f"gdal_translate failed:\n{result.stderr}")
        return False
    return True


# ---------------------------------------------------------------------------
# Single file conversion
# ---------------------------------------------------------------------------

def convert_one(
    eeh_path: str,
    geo_index: dict,
    out_dir: str,
    cloud_index: dict | None = None,
    max_cloud_cover: float | None = None,
    force_product: str | None = None,
    force: bool = False,
) -> str | None:
    fname = Path(eeh_path).name
    product = force_product or detect_product(fname)
    if product is None:
        return None

    if force_product == "WUE":
        cog_name = fname.replace("STIC_L3_GPP", "STIC_L4_WUE").replace(".h5", ".tif")
    else:
        cog_name = fname.replace(".h5", ".tif")

    cog_out = os.path.join(out_dir, cog_name)
    if not force and os.path.isfile(cog_out):
        log.info(f"  Skip (already exists): {cog_name}")
        return cog_name

    eeh_key = extract_key(eeh_path, EEH_RE)
    if eeh_key is None:
        log.warning(f"Cannot extract orbit key from {fname}")
        return None

    geo_path = geo_index.get(eeh_key)
    if geo_path is None:
        log.warning(f"GEO file missing for key '{eeh_key}' ({fname})")
        return None

    work_dir = Path(TMP_BASE) / eeh_key
    work_dir.mkdir(parents=True, exist_ok=True)
    local_gtiff = work_dir / cog_name.replace(".tif", "_tmp.tif")

    try:
        log.info(f"  Processing: {fname} (product={product})")

        # -- Read GEO -------------------------------------------------------
        with h5py.File(geo_path, "r") as f_geo:
            lat = f_geo["Geolocation/latitude"][:].astype("float32")
            lon = f_geo["Geolocation/longitude"][:].astype("float32")
            land_fraction = f_geo["Geolocation/land_fraction"][:].astype("float32")
            day_night_flag = f_geo["StandardMetadata/DayNightFlag"][()].decode().strip()
        log.info(f"    day_night_flag: {day_night_flag}")

        # -- Cloud mask -----------------------------------------------------
        cloud_mask = None
        cloud_path = (cloud_index or {}).get(eeh_key)
        if cloud_path is not None:
            with h5py.File(cloud_path, "r") as f_cld:
                cld_raw = np.array(f_cld["SDS/Cloud_final"])
            cloud_mask = (cld_raw == 0).astype("float32")

            if max_cloud_cover is not None and product == "TES":
                valid_cld = cld_raw != 255
                if valid_cld.any():
                    cloud_pct = float((1 - np.mean(cloud_mask[valid_cld])) * 100)
                    if cloud_pct > max_cloud_cover:
                        log.info(f"  Skip cloud ({cloud_pct:.1f}% > {max_cloud_cover}%): {fname}")
                        return None
                    log.info(f"    cloud_cover: {cloud_pct:.1f}%")
            del cld_raw
            gc.collect()

        # -- Bounding box ---------------------------------------------------
        valid = np.isfinite(lat) & np.isfinite(lon)
        west = float(np.nanmin(lon[valid]))
        east = float(np.nanmax(lon[valid]))
        south = float(np.nanmin(lat[valid]))
        north = float(np.nanmax(lat[valid]))

        # -- Snapped grid ---------------------------------------------------
        RES64 = np.float64(RES)
        west_snap = round(west / RES) * RES64
        north_snap = round(north / RES) * RES64
        east_snap = round(east / RES) * RES64
        south_snap = round(south / RES) * RES64
        ncols = int(round((east_snap - west_snap) / RES64))
        nrows = int(round((north_snap - south_snap) / RES64))
        log.info(f"    grid: cols={ncols} rows={nrows}")

        lons_1d = (west_snap + np.arange(ncols, dtype="float64") * RES64).astype("float32")
        lats_1d = (north_snap - np.arange(nrows, dtype="float64") * RES64).astype("float32")
        grid_lon, grid_lat = np.meshgrid(lons_1d, lats_1d)
        transform = from_origin(west_snap, north_snap, RES64, RES64)

        # -- kd_tree --------------------------------------------------------
        swath_def = geometry.SwathDefinition(lons=lon, lats=lat)
        grid_def = geometry.GridDefinition(lons=grid_lon, lats=grid_lat)

        valid_in, valid_out, idx_arr, dist_arr = kd_tree.get_neighbour_info(
            swath_def, grid_def,
            radius_of_influence=RADIUS_INFLUENCE,
            neighbours=1,
        )
        grid_shape = grid_def.shape

        del lat, lon, valid, lons_1d, lats_1d
        del grid_lon, grid_lat, swath_def, grid_def, dist_arr
        gc.collect()

        # -- Write intermediate GTiff ---------------------------------------
        bands_cfg = PRODUCT_BANDS[product]
        log.info(f"    writing GTiff ({len(bands_cfg)} bands)...")

        with h5py.File(eeh_path, "r") as f_h5:
            with rasterio.open(
                str(local_gtiff), "w",
                driver="GTiff",
                tiled=True,
                height=nrows,
                width=ncols,
                count=len(bands_cfg),
                dtype="float32",
                crs=CRS.from_epsg(4326),
                transform=transform,
                nodata=NODATA,
            ) as dst:
                for i, bc in enumerate(bands_cfg, start=1):
                    if bc["name"] == "land_fraction":
                        raw = land_fraction
                    elif bc["name"] == "cloud_mask":
                        if cloud_mask is not None:
                            raw = cloud_mask
                        else:
                            dst.write(np.full((nrows, ncols), NODATA, "float32"), i)
                            dst.set_band_description(i, bc["name"])
                            dst.update_tags(i, name=bc["name"], long_name=bc["long_name"], units=bc["units"])
                            continue
                    elif bc["name"] not in f_h5:
                        log.warning(f"    Band '{bc['name']}' missing in {fname}")
                        dst.write(np.full((nrows, ncols), NODATA, "float32"), i)
                        dst.set_band_description(i, bc["name"])
                        dst.update_tags(i, name=bc["name"], long_name=bc["long_name"], units=bc["units"])
                        continue
                    else:
                        raw = f_h5[bc["name"]][:].astype("float32")
                        raw[raw == bc["nodata"]] = np.nan
                        raw = raw * bc["scale"] + bc["offset"]

                    resampled = kd_tree.get_sample_from_neighbour_info(
                        "nn", grid_shape, raw,
                        valid_in, valid_out, idx_arr,
                        fill_value=np.nan,
                    ).astype("float32")
                    resampled[~np.isfinite(resampled)] = NODATA

                    dst.write(resampled, i)
                    dst.set_band_description(i, bc["name"])
                    dst.update_tags(i, name=bc["name"], long_name=bc["long_name"], units=bc["units"])

                    del raw, resampled
                    gc.collect()

                dst.update_tags(
                    product=product,
                    source=fname,
                    eeh_key=eeh_key,
                    epsg="4326",
                    bands=",".join(b["name"] for b in bands_cfg),
                    bbox=f"{west},{south},{east},{north}",
                    day_night_flag=day_night_flag,
                )

        del valid_in, valid_out, idx_arr, land_fraction, cloud_mask
        gc.collect()

        # -- GTiff → COG ----------------------------------------------------
        log.info(f"    converting to COG via gdal_translate...")
        os.makedirs(out_dir, exist_ok=True)
        if not gtiff_to_cog(str(local_gtiff), cog_out):
            return None

        log.info(f"    wrote: {cog_out}")
        return cog_name

    except Exception as e:
        log.error(f"ERROR {fname}: {e}", exc_info=True)
        return None

    finally:
        if work_dir.exists():
            shutil.rmtree(str(work_dir), ignore_errors=True)
        gc.collect()


# ---------------------------------------------------------------------------
# Batch conversion (called from run_pipeline.py or CLI)
# ---------------------------------------------------------------------------

def convert_batch(
    eeh_files: list[str],
    geo_dir: str,
    out_dir: str,
    cloud_dir: str | None = None,
    max_cloud_cover: float | None = None,
    force_product: str | None = None,
    force: bool = False,
    workers: int = 4,
) -> tuple[int, int, int]:
    """Convert a list of EEH H5 files to COGs. Returns (success, skipped, errors)."""
    Path(TMP_BASE).mkdir(parents=True, exist_ok=True)
    os.makedirs(out_dir, exist_ok=True)

    geo_index = build_geo_index(geo_dir)
    cloud_index = build_cloud_index(cloud_dir) if cloud_dir else None

    if force_product:
        log.info(f"--force-product active: treating all files as '{force_product}'")

    success = errors = skipped = 0

    with ProcessPoolExecutor(max_workers=workers) as executor:
        futures = {
            executor.submit(
                convert_one, f, geo_index, out_dir,
                cloud_index, max_cloud_cover, force_product, force,
            ): f
            for f in eeh_files
        }
        for future in as_completed(futures):
            try:
                result = future.result()
                if result:
                    success += 1
                else:
                    skipped += 1
            except Exception as e:
                log.error(f"Future exception: {e}")
                errors += 1

    log.info(f"COG conversion done — ok {success}  skip {skipped}  err {errors}")
    return success, skipped, errors


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Convert EEH H5 files to Cloud Optimized GeoTIFFs.",
    )
    parser.add_argument("--eeh-dir", required=True,
                        help="Directory containing EEH H5 files.")
    parser.add_argument("--geo-dir", required=True,
                        help="Directory containing GEO H5 files.")
    parser.add_argument("--cloud-dir", default=None,
                        help="Directory containing CLOUD H5 files (for TES cloud masking).")
    parser.add_argument("--out-dir", required=True,
                        help="Output directory for COG files.")
    parser.add_argument("--force-product", default=None,
                        choices=["TES", "STIC", "GPP", "WUE"],
                        help="Override detected product (e.g. WUE from GPP H5 files).")
    parser.add_argument("--max-cloud-cover", type=float, default=None,
                        help="Skip TES scenes exceeding this cloud cover %% (e.g. 70).")
    parser.add_argument("--force", action="store_true",
                        help="Re-convert even if COG already exists.")
    parser.add_argument("--workers", type=int, default=4,
                        help="Number of parallel workers (default: 4).")
    parser.add_argument("--limit", type=int, default=None,
                        help="Process only the first N files (for testing).")
    args = parser.parse_args()

    eeh_files = sorted(
        str(Path(args.eeh_dir) / f)
        for f in os.listdir(args.eeh_dir)
        if f.endswith(".h5") and detect_product(f)
    )
    if args.limit:
        eeh_files = eeh_files[: args.limit]

    log.info(f"EEH files to process: {len(eeh_files)}")
    convert_batch(
        eeh_files,
        args.geo_dir,
        args.out_dir,
        cloud_dir=args.cloud_dir,
        max_cloud_cover=args.max_cloud_cover,
        force_product=args.force_product,
        force=args.force,
        workers=args.workers,
    )


if __name__ == "__main__":
    main()
