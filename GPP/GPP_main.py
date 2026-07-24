#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Script to run the hybrid model for GPP and WUE estimation
Created on May 1 2025
@author: Ziyu Lin at LIST
© 2025 - Luxembourg Institute of Science and Technology
Authors : Ziyu Lin, Kaniska Mallick, Tian Hu (tian.hu@list.lu)
Code licensed under MIT
SPDX-License-Identifier: MIT

Modified to run fully locally (no WASDI platform dependency):
- Removed all wasdi.* calls -> replaced by plain Python (print) and an
  argparse-based CLI, mirroring STIC_main.py / TES_main.py.
- Removed cache_S3_with_pattern() / cache_S3_simple() / S3_cache.py
  dependency:
    * cache_S3_with_pattern() -> build_local_pattern_map() (local, regex
      based directory scan, {captured_group: filename} mapping).
    * cache_S3_simple()       -> list_local_files_cached() (local, glob
      based directory scan, returns a plain sorted list of full paths).
- Performance: directory scans use os.scandir() and check the
  pattern/fnmatch BEFORE calling is_file() (avoids a stat() syscall for
  every non-matching entry on slow/networked storage); the 8 independent
  directory scans (GEO, CLOUD, PARH, LAI, FVC, ERA5, OCO2, GLC30) are run
  in parallel threads via build_all_caches(); results are cached to disk
  as JSON so that slow/networked storage only needs to be scanned once
  (use --force-rebuild-cache to force a re-scan).
- GEO/CLOUD default patterns are anchored with '$' to correctly exclude
  sidecar files (.xml, .dmrpp, etc.) that would otherwise collide with the
  real .h5 file under the same lookup key.

NOTE: the original WASDI script also cached MOTA (blue-sky albedo) and CI
(clumping index) file listings, but never actually passed them to
run_GPP() -- they were unused dead code. They are intentionally NOT
reproduced here. Reintroduce them (following the same
list_local_files_cached() pattern) if a future version of the model
needs them.
"""
import argparse
import fnmatch
import glob
import hashlib
import json
import os
import re
import time
import traceback
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pandas as pd
import xarray as xr

from ReadData_EEH_GPP_final_TH import (
    regex_extract_ecostress, Create_cloud_mask, read_geo, Extract_STIC,
    extract_PAR_from_global_data, wrap_to_ECOSTRESS, interpolate_LAI_value,
    interpolate_FVC_value, interpolate_ERA5_value, read_oco2,
    filter_tile_bounds, read_and_wrap_GLC30, big_leaf_calculator,
    gpp_gs_calculator, generate_hdf5_file_GPP, GLC30_mapper
)

import matplotlib
matplotlib.use('Agg')


# ---------------------------------------------------------------------------
# Local file-name caching (replaces S3_cache.py)
# ---------------------------------------------------------------------------
def _cache_file_path(cache_dir, label, directory, pattern):
    """Build a stable, unique cache filename for a given (directory, pattern) pair."""
    key = hashlib.sha1(f"{directory}|{pattern}".encode('utf-8')).hexdigest()[:16]
    safe_label = re.sub(r'[^A-Za-z0-9_-]', '_', label or 'cache')
    return os.path.join(cache_dir, f'pattern_cache_{safe_label}_{key}.json')


def build_local_pattern_map(directory, pattern, label=None,
                             cache_dir=None, force_rebuild=False):
    """
    Local, pure-Python replacement for cache_S3_with_pattern().
    Scans `directory`, applies re.search(pattern, filename) to every entry,
    and returns {captured_group_1: actual_filename}. Results are optionally
    cached to disk (JSON) to avoid rescanning slow/networked storage on
    every run.

    When multiple files share the same key (e.g. a .h5 file and its
    .h5.xml / .h5.dmrpp sidecars, if the pattern is not anchored with '$'),
    the file whose name is lexicographically greatest is kept
    deterministically, with a log message.

    Parameters
    ----------
    directory : str
        Local folder to scan.
    pattern : str
        Regular expression with exactly one capturing group, applied with
        re.search(). Anchor with '$' if you need to match a specific file
        extension exactly (e.g. to exclude '.h5.xml'/'.h5.dmrpp' sidecars).
    label : str, optional
        Cosmetic label used in log messages and in the cache filename.
    cache_dir : str, optional
        Directory where a JSON cache of the resulting {key: filename} map
        is stored/read. If None, no caching is performed.
    force_rebuild : bool
        If True, ignore any existing cache file and rescan the directory.

    Returns
    -------
    dict : {captured_group: actual_filename}
    """
    cache_path = None
    if cache_dir:
        os.makedirs(cache_dir, exist_ok=True)
        cache_path = _cache_file_path(cache_dir, label, directory, pattern)
        if not force_rebuild and os.path.isfile(cache_path):
            try:
                with open(cache_path, 'r') as f:
                    file_map = json.load(f)
                print(f'[{label}] Using cached file-name map '
                      f'({len(file_map)} entries) from {cache_path}')
                return file_map
            except Exception as e:
                print(f'[{label}] Warning: failed to read cache {cache_path} '
                      f'({e}) -> rescanning directory')

    if not os.path.isdir(directory):
        raise FileNotFoundError(f"Directory not found for '{label}': {directory}")
    print(f'[{label}] Scanning {directory} (this may take a while on slow/networked storage)...')
    regex = re.compile(pattern)
    with os.scandir(directory) as it:
        entries = sorted(it, key=lambda e: e.name)
    file_map = {}
    for entry in entries:
        fname = entry.name
        # Cheap in-memory regex check first: avoids a stat() syscall for
        # every entry that doesn't even match the naming pattern.
        m = regex.search(fname)
        if not m:
            continue
        if not entry.is_file(follow_symlinks=False):
            continue
        key = m.group(1)
        if key in file_map:
            if fname > file_map[key]:
                print(f'[{label}] Note: preferring more recent "{fname}" over '
                      f'"{file_map[key]}" for key "{key}"')
                file_map[key] = fname
            continue
        file_map[key] = fname
    if not file_map:
        print(f'[{label}] Warning: no files matching pattern found in {directory}')
    else:
        print(f'[{label}] Found {len(file_map)} matching files in {directory}')
    if cache_path:
        try:
            with open(cache_path, 'w') as f:
                json.dump(file_map, f)
            print(f'[{label}] Cached {len(file_map)} entries to {cache_path}')
        except Exception as e:
            print(f'[{label}] Warning: failed to write cache {cache_path}: {e}')
    return file_map


def list_local_files_cached(directory, pattern='*', label=None,
                             cache_dir=None, force_rebuild=False):
    """
    Local, pure-Python replacement for cache_S3_simple(). Lists the
    absolute paths of files in `directory` matching a glob-style `pattern`
    (fnmatch), with the same caching/scandir optimisations as
    build_local_pattern_map().

    Parameters
    ----------
    directory : str
        Local folder to scan.
    pattern : str
        Glob/fnmatch pattern (e.g. '*.nc'). Default '*' matches every file.
    label : str, optional
        Cosmetic label used in log messages and in the cache filename.
    cache_dir : str, optional
        Directory where a JSON cache of the resulting file list is
        stored/read. If None, no caching is performed.
    force_rebuild : bool
        If True, ignore any existing cache file and rescan the directory.

    Returns
    -------
    list of str : sorted absolute file paths matching `pattern`.
    """
    cache_path = None
    if cache_dir:
        os.makedirs(cache_dir, exist_ok=True)
        cache_path = _cache_file_path(cache_dir, label, directory, pattern)
        if not force_rebuild and os.path.isfile(cache_path):
            try:
                with open(cache_path, 'r') as f:
                    files = json.load(f)
                print(f'[{label}] Using cached file list '
                      f'({len(files)} entries) from {cache_path}')
                return files
            except Exception as e:
                print(f'[{label}] Warning: failed to read cache {cache_path} '
                      f'({e}) -> rescanning directory')

    if not os.path.isdir(directory):
        raise FileNotFoundError(f"Directory not found for '{label}': {directory}")
    print(f'[{label}] Scanning {directory} (this may take a while on slow/networked storage)...')
    with os.scandir(directory) as it:
        entries = sorted(it, key=lambda e: e.name)
    files = []
    for entry in entries:
        fname = entry.name
        # Cheap in-memory fnmatch check first: avoids a stat() syscall for
        # every entry that doesn't even match the pattern.
        if not fnmatch.fnmatch(fname, pattern):
            continue
        if not entry.is_file(follow_symlinks=False):
            continue
        files.append(os.path.join(directory, fname))
    if not files:
        print(f'[{label}] Warning: no files matching pattern found in {directory}')
    else:
        print(f'[{label}] Found {len(files)} matching files in {directory}')
    if cache_path:
        try:
            with open(cache_path, 'w') as f:
                json.dump(files, f)
            print(f'[{label}] Cached {len(files)} entries to {cache_path}')
        except Exception as e:
            print(f'[{label}] Warning: failed to write cache {cache_path}: {e}')
    return files


def build_all_caches(geo_dir, cld_dir, parh_dir, lai_dir, fvc_dir, era5_dir, oco2_dir, glc30_dir,
                      geo_pattern, cloud_pattern, parh_pattern, lai_pattern, fvc_pattern,
                      era5_pattern, oco2_pattern, glc30_pattern,
                      cache_dir=None, force_rebuild=False):
    """
    Build the 8 independent file-name/listing lookups (GEO, CLOUD, PARH,
    LAI, FVC, ERA5, OCO2, GLC30) in parallel threads, since each is a
    purely I/O-bound directory scan with no dependency on the others.

    Returns
    -------
    dict with keys 'geo', 'cloud', 'parh', 'lai', 'fvc', 'era5', 'oco2',
    'glc30'. 'geo' and 'cloud' are {key: filename} dicts; the rest are
    plain sorted lists of full file paths.
    """
    jobs = {
        'geo':   (build_local_pattern_map, geo_dir, geo_pattern, 'L1B_GEO'),
        'cloud': (build_local_pattern_map, cld_dir, cloud_pattern, 'CLOUD'),
        'parh':  (list_local_files_cached, parh_dir, parh_pattern, 'PARH'),
        'lai':   (list_local_files_cached, lai_dir, lai_pattern, 'LAI'),
        'fvc':   (list_local_files_cached, fvc_dir, fvc_pattern, 'FCOVER'),
        'era5':  (list_local_files_cached, era5_dir, era5_pattern, 'ERA5'),
        'oco2':  (list_local_files_cached, oco2_dir, oco2_pattern, 'OCO2'),
        'glc30': (list_local_files_cached, glc30_dir, glc30_pattern, 'GLC30'),
    }
    results = {}
    with ThreadPoolExecutor(max_workers=len(jobs)) as executor:
        futures = {
            name: executor.submit(func, directory, pattern, label=label,
                                   cache_dir=cache_dir, force_rebuild=force_rebuild)
            for name, (func, directory, pattern, label) in jobs.items()
        }
        for name, fut in futures.items():
            results[name] = fut.result()  # re-raises any exception from the thread
    return results


# ---------------------------------------------------------------------------
# Core processing (unchanged logic, wasdi.* -> print / return)
# ---------------------------------------------------------------------------
def run_GPP(stic_file, writepath, map_cld, map_geo, files_PARH, files_lai, files_fvc,
            files_era5, files_oco2, files_GLC30, directory_cld, directory_geo, df_LUT):
    print(stic_file)
    # Extract orbit & time details
    orbit_str, scene_str, year_str, month_str, day_str, hour_str, min_str, sec_str = regex_extract_ecostress(stic_file)
    year, month, day, hour, minute = int(year_str), int(month_str), int(day_str), int(hour_str), int(min_str)
    doy = pd.Timestamp(year=year, month=month, day=day).dayofyear
    doy_str = str(int(doy))

    key_eco = f"{orbit_str}_{scene_str}_{year_str}{month_str}{day_str}T{hour_str}{min_str}{sec_str}"
    orbit_id = f"{orbit_str}_{scene_str}"

    geo_key = 'ECOv002_L1B_GEO_' + key_eco
    cld_key = 'ECOv002_L2_CLOUD_' + key_eco
    if geo_key not in map_geo or cld_key not in map_cld:
        print("Corresponding GEO/CLOUD data not available")
        return
    path_geo = os.path.join(directory_geo, map_geo[geo_key])
    path_cld = os.path.join(directory_cld, map_cld[cld_key])

    if not all([path_geo, path_cld, stic_file]):
        print("Corresponding data not available")
        return

    # exit if cloudy image (>75%) to save time
    try:
        start_time = time.time()
        mask_cld = Create_cloud_mask(path_cld)
        print(f"Reading Cloud file done in {time.time() - start_time:.2f} seconds")
        cld_percent = 1 - np.mean(mask_cld)
        if cld_percent >= 0.75:
            print(f'too many cloud cover {int(cld_percent * 100)}%, pass')
            return
    except Exception as e:
        print(f"Error processing Cloud file: {e}")
        return

    # get lat, lon, SZA from ECOSTRESS metadata
    eco_lat, eco_lon, eco_sza = read_geo(path_geo)
    sza_max = np.max(eco_sza)
    if sza_max > 90:
        print(f'Skip night time with high SZA:{sza_max}')
        return

    lat_max, lat_min, lon_max, lon_min = (np.nanmax(eco_lat), np.nanmin(eco_lat),
                                           np.nanmax(eco_lon), np.nanmin(eco_lon))
    eco_bound = (lat_max + 0.5, lat_min - 0.5, lon_max + 0.5, lon_min - 0.5)

    # STIC-gs and ET from EEH2 product
    try:
        start_time = time.time()
        (ETD, gt) = Extract_STIC(stic_file, mask_cld)
        valid_percent = np.mean(ETD > 0)
        print(valid_percent)
        if valid_percent == 0:
            print("There is no valid ETD>0, GPP output will be full of NaN")
            return
        print(f"{stic_file} \nExtracting STIC results done in {time.time() - start_time:.2f} seconds")
    except Exception as e:
        print(f"Error extracting STIC result: {e}")
        return

    #####################################################################################
    # STEP2 read and rescale other inputs
    #####################################################################################
    path_PARH_matches = [p for p in files_PARH if f"PARin{year_str}{month_str}{day_str}" in p]
    if not path_PARH_matches:
        print(f"No PARH file found for {year_str}{month_str}{day_str}")
        return
    path_PARH = path_PARH_matches[0]
    PARH, PARmean = extract_PAR_from_global_data(path_PARH, eco_bound, hour, minute)
    PARH_70m = wrap_to_ECOSTRESS(PARH, eco_lat, eco_lon)
    PARDmean_70m = wrap_to_ECOSTRESS(PARmean, eco_lat, eco_lon)
    print("PARH done")

    LAI_300m = interpolate_LAI_value(year, month, day, files_lai, eco_bound)
    LAI_70m = wrap_to_ECOSTRESS(LAI_300m, eco_lat, eco_lon)
    FVC_300m = interpolate_FVC_value(year, month, day, files_fvc, eco_bound)
    FVC_70m = wrap_to_ECOSTRESS(FVC_300m, eco_lat, eco_lon)
    print("LAI FVC done")

    files_era5_single = [f for f in files_era5 if 'single' in f]
    t2m_ERA5_10km, vpd_ERA5_10km = interpolate_ERA5_value(year, month, day, hour, minute, files_era5_single, eco_bound)
    temp_70m = wrap_to_ECOSTRESS(t2m_ERA5_10km, eco_lat, eco_lon)
    vpd_70m = wrap_to_ECOSTRESS(vpd_ERA5_10km, eco_lat, eco_lon)
    print("ERA5 done")

    oco2_50km = read_oco2(year, month, day, eco_bound, files_oco2)
    oco2_70m = wrap_to_ECOSTRESS(oco2_50km, eco_lat, eco_lon)
    print("OCO2 done")

    list_path_GLC30 = [p for p in files_GLC30 if filter_tile_bounds(p, eco_bound)]
    glc70m = read_and_wrap_GLC30(list_path_GLC30, eco_lat, eco_lon, eco_bound, year)
    glc70m = glc70m.fillna(0).astype(int)
    print("GLC30 done")

    inputs = {
        "ETD": ETD,
        'gt': gt,
        "LAI": LAI_70m.data,
        "FVC": FVC_70m.data,
        'Tair': temp_70m.data,
        'VPD': vpd_70m.data,
        "PARH": PARH_70m.data,
        "PARDmean": PARDmean_70m.data,
        "CO2": oco2_70m.data,
        "LUCC": glc70m.data
    }

    code_to_pft = {v: k for k, v in GLC30_mapper.items()}
    unique_codes = np.unique(glc70m.values[~np.isnan(glc70m.values)]).astype(int)
    if len(unique_codes) == 1 and unique_codes[0] == -9999:
        print("No vegetation tile -> ignored")
        return
    elif -9999 in unique_codes:
        unique_codes = np.where(unique_codes == -9999, 0, unique_codes)

    valid_pfts = [code_to_pft.get(code) for code in unique_codes if code in code_to_pft]
    df_LUT_filtered = df_LUT[df_LUT['type'].isin(valid_pfts)]
    ds_BL = {}
    BL_cols = [c for c in df_LUT.columns if 'BL' in c and 'gsFULL' in c and 'train' not in c]

    for col in BL_cols:
        lookup = np.zeros(max(unique_codes) + 1)
        for code in unique_codes:
            pft = code_to_pft.get(code)
            if pft in df_LUT_filtered.set_index('type').index:
                value = df_LUT_filtered.set_index('type').loc[pft][col]
                lookup[code] = value
            else:
                lookup[code] = np.nan
        mapped_array = xr.apply_ufunc(
            lambda x: lookup[x],
            glc70m,
            vectorize=True,
            dask='allowed'
        )
        col_short = col.split('_')[-1]
        da = mapped_array
        da.name = col_short
        ds_BL[col_short] = da.data

    #####################################################################################
    # STEP3 calculate hourly GPP and daily GPP&WUE
    #####################################################################################
    s_to_day = 60 * 60 * 24

    vector_BL = (ds_BL['LUEmax'], ds_BL['Tmax'], ds_BL['Tmin'], ds_BL['Topt'], ds_BL['k'])
    GPP_LUE = big_leaf_calculator(inputs, vector_BL).clip(0, 0.001)
    GPP_gs = gpp_gs_calculator(inputs['gt'], inputs['CO2'], inputs['Tair'], inputs['VPD']).clip(0, 0.001)

    vpd_mask = (vpd_70m > 20)
    GPPmin = np.where(vpd_mask, GPP_gs, GPP_LUE).clip(0, 0.001)

    veg_p25_mask = (inputs['FVC'] > 0.25)
    eco_valid_mask = (ETD > 0.01) & (~np.isnan(ETD)) & veg_p25_mask
    GPP_LUE, GPP_gs, GPPmin, ETD_safe = (
        np.where(eco_valid_mask, GPP_LUE, np.nan),
        np.where(eco_valid_mask, GPP_gs, np.nan),
        np.where(eco_valid_mask, GPPmin, np.nan),
        np.where(eco_valid_mask, ETD, np.nan),
    )

    PAR_scaler_raw = inputs['PARDmean'] / inputs['PARH']
    GPPd_gC = ((GPPmin * PAR_scaler_raw) * s_to_day).clip(0, 50)

    WUE_raw = (GPPd_gC / ETD_safe)
    LUCC = inputs['LUCC']

    outputs = {
        'GPPd': np.nan_to_num(GPPd_gC, nan=-9999),
        'WUEd': np.nan_to_num(WUE_raw, nan=-9999),
        'GPP_gs': np.nan_to_num(GPP_gs, nan=-9999),
        'GPP_LUE': np.nan_to_num(GPP_LUE, nan=-9999),
        'GPPmin': np.nan_to_num(GPPmin, nan=-9999),
        'LUCC': np.nan_to_num(LUCC, nan=-9999)
    }
    generate_hdf5_file_GPP(writepath, outputs, key_eco)


# ---------------------------------------------------------------------------
# Batch driver
# ---------------------------------------------------------------------------
def run_batch(stic_files, output_path, geo_dir, cld_dir, parh_dir, lai_dir, fvc_dir,
              era5_dir, oco2_dir, glc30_dir, df_LUT,
              geo_pattern, cloud_pattern, parh_pattern, lai_pattern, fvc_pattern,
              era5_pattern, oco2_pattern, glc30_pattern,
              cache_dir=None, force_rebuild_cache=False, error_log_dir='.'):
    """
    Build all 8 file-name/listing caches once (in parallel), then loop over
    the STIC files (mirrors run_STIC's per-file try/except and JSON error
    log pattern).
    """
    map_error = {}
    print('Building local file-name caches (parallel scan)...')
    caches = build_all_caches(
        geo_dir, cld_dir, parh_dir, lai_dir, fvc_dir, era5_dir, oco2_dir, glc30_dir,
        geo_pattern, cloud_pattern, parh_pattern, lai_pattern, fvc_pattern,
        era5_pattern, oco2_pattern, glc30_pattern,
        cache_dir=cache_dir, force_rebuild=force_rebuild_cache,
    )
    print('Done caching.')

    map_geo = caches['geo']
    map_cld = caches['cloud']
    parh_files = caches['parh']
    lai_files = caches['lai']
    fvc_files = caches['fvc']
    era5_files = caches['era5']
    oco2_files = caches['oco2']
    glc30_files = caches['glc30']

    os.makedirs(output_path, exist_ok=True)

    for stic_file in stic_files:
        try:
            run_GPP(stic_file, output_path, map_cld, map_geo, parh_files, lai_files, fvc_files,
                    era5_files, oco2_files, glc30_files, cld_dir, geo_dir, df_LUT)
        except Exception:
            traceback.print_exc()
            print(f'Error while processing {stic_file} -> skipping')
            map_error[stic_file] = traceback.format_exc()

    if map_error:
        pid = str(os.getpid())
        os.makedirs(error_log_dir, exist_ok=True)
        error_log_path = os.path.join(error_log_dir, f'execution_{pid}_errors.txt')
        with open(error_log_path, 'w') as f:
            f.write(json.dumps(map_error))
        print(f'Some files failed to process, see {error_log_path}')


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def build_arg_parser():
    parser = argparse.ArgumentParser(
        prog='GPP_main.py',
        description='Run the hybrid model for GPP and WUE estimation on ECOSTRESS STIC files.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        '--stic-files', '-i', nargs='+', metavar='STIC_FILE',
        help='One or more EEH2STIC L3 ET/STIC .h5 files to process (STIC output).'
    )
    input_group.add_argument(
        '--stic-dir', metavar='DIR',
        help='Directory to scan for STIC files (used together with --stic-pattern).'
    )
    parser.add_argument(
        '--stic-pattern', default='*EEH2STIC*.h5',
        help='Glob pattern used to select STIC files when --stic-dir is given.'
    )
    parser.add_argument('--output-dir', required=True,
                         help='Directory where the output GPP/WUE HDF5 files are written.')
    parser.add_argument('--geo-dir', required=True,
                         help='Directory containing the ECOSTRESS L1B_GEO .h5 files.')
    parser.add_argument('--cloud-dir', required=True,
                         help='Directory containing the ECOSTRESS L2_CLOUD .h5 files.')
    parser.add_argument('--parh-dir', required=True,
                         help='Directory containing the half-hourly PAR (PARin*) files.')
    parser.add_argument('--lai-dir', required=True,
                         help='Directory containing the LAI (300m, 10-day) files.')
    parser.add_argument('--fvc-dir', required=True,
                         help='Directory containing the FCOVER (300m, 10-day) files.')
    parser.add_argument('--era5-dir', required=True,
                         help='Directory containing the ERA5 single-level files.')
    parser.add_argument('--oco2-dir', required=True,
                         help='Directory containing the OCO-2 XCO2 files.')
    parser.add_argument('--glc30-dir', required=True,
                         help='Directory containing the GLC30 land-cover tiles.')
    parser.add_argument('--lut-file', required=True,
                         help='Path to the LookUpTable_LUE_HH_gsFULL-globe.csv file.')

    pattern_group = parser.add_argument_group(
        'File-name matching patterns (override only if your local filenames '
        'differ from the standard ECOSTRESS/CGLS/MODIS conventions). GEO/CLOUD '
        'patterns are anchored to the end of the filename ($) to correctly '
        'exclude sidecar files such as .xml/.dmrpp.')
    pattern_group.add_argument(
        '--geo-pattern', default=r'(ECOv002_L1B_GEO.*)_\d{4}_\d{2}\.h5$',
        help='Regex (one capturing group) for ECOSTRESS L1B_GEO filenames.')
    pattern_group.add_argument(
        '--cloud-pattern', default=r'(ECOv002_L2_CLOUD.*)_\d{4}_\d{2}\.h5$',
        help='Regex (one capturing group) for ECOSTRESS L2_CLOUD filenames.')
    pattern_group.add_argument(
        '--parh-pattern', default='*',
        help='Glob pattern for hourly PAR (PARin*) files.')
    pattern_group.add_argument(
        '--lai-pattern', default='*',
        help='Glob pattern for LAI files.')
    pattern_group.add_argument(
        '--fvc-pattern', default='*',
        help='Glob pattern for FCOVER files.')
    pattern_group.add_argument(
        '--era5-pattern', default='*',
        help='Glob pattern for ERA5 files.')
    pattern_group.add_argument(
        '--oco2-pattern', default='*',
        help='Glob pattern for OCO-2 XCO2 files.')
    pattern_group.add_argument(
        '--glc30-pattern', default='*',
        help='Glob pattern for GLC30 land-cover tiles.')

    cache_group = parser.add_argument_group('File-name cache (avoids rescanning slow/networked storage)')
    cache_group.add_argument(
        '--cache-dir', default=None,
        help=('Directory used to persist the GEO/CLOUD/PARH/LAI/FCOVER/ERA5/OCO2/GLC30 '
              'file-name lookup maps to disk (JSON). If not given, defaults to '
              '<output-dir>/.gpp_filename_cache. Pass an empty string ("") to '
              'disable caching entirely and always rescan.')
    )
    cache_group.add_argument(
        '--force-rebuild-cache', action='store_true',
        help='Ignore any existing file-name cache and rescan all directories.'
    )
    parser.add_argument('--error-log-dir', default='.',
                         help='Directory to write the execution error JSON log to, if any.')
    return parser


def main(argv=None):
    parser = build_arg_parser()
    args = parser.parse_args(argv)

    if args.stic_files:
        stic_files = args.stic_files
    else:
        stic_files = sorted(glob.glob(os.path.join(args.stic_dir, args.stic_pattern)))
        if not stic_files:
            parser.error(f"No files matching '{args.stic_pattern}' found in {args.stic_dir}")

    df_LUT = pd.read_csv(args.lut_file)
    if 'nonV' not in df_LUT.index:
        new_row = pd.DataFrame([0] * len(df_LUT.columns), index=df_LUT.columns).T
        new_row['type'] = 'nonV'
        df_LUT = pd.concat([df_LUT, new_row])

    os.makedirs(args.output_dir, exist_ok=True)

    if args.cache_dir is None:
        cache_dir = os.path.join(args.output_dir, '.gpp_filename_cache')
    elif args.cache_dir == '':
        cache_dir = None  # caching explicitly disabled
    else:
        cache_dir = args.cache_dir

    run_batch(
        stic_files=stic_files,
        output_path=args.output_dir,
        geo_dir=args.geo_dir,
        cld_dir=args.cloud_dir,
        parh_dir=args.parh_dir,
        lai_dir=args.lai_dir,
        fvc_dir=args.fvc_dir,
        era5_dir=args.era5_dir,
        oco2_dir=args.oco2_dir,
        glc30_dir=args.glc30_dir,
        df_LUT=df_LUT,
        geo_pattern=args.geo_pattern,
        cloud_pattern=args.cloud_pattern,
        parh_pattern=args.parh_pattern,
        lai_pattern=args.lai_pattern,
        fvc_pattern=args.fvc_pattern,
        era5_pattern=args.era5_pattern,
        oco2_pattern=args.oco2_pattern,
        glc30_pattern=args.glc30_pattern,
        cache_dir=cache_dir,
        force_rebuild_cache=args.force_rebuild_cache,
        error_log_dir=args.error_log_dir,
    )


if __name__ == '__main__':
    main()
