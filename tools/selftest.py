#!/usr/bin/env python3
"""
EEH2 installation self-test
===========================
Two independent stages:

  1. environment — Python version, every third-party module the three steps
     import, the RTTOV installation and its ECOSTRESS coefficient file, the GPP
     look-up table, the data and output directories, and which credentials are
     configured. Needs no data and runs in seconds.

  2. products — statistics of the pipeline outputs compared against the
     reference case shipped in tools/reference/. Answers "did my installation
     reproduce the validated results?", which the environment stage cannot.

Usage (in the container):
    docker compose run selftest
    docker compose run selftest --products

Usage (standalone):
    python tools/selftest.py
    python tools/selftest.py --products
    python tools/selftest.py --products --write-reference
"""

from __future__ import annotations

import argparse
import importlib
import json
import logging
import os
import sys
from pathlib import Path

import h5py
import numpy as np

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
)
log = logging.getLogger(__name__)

REFERENCE_FILE = Path(__file__).parent / "reference" / "reference_36798_20250101.json"

MIN_PYTHON = (3, 10)

# Third-party imports of the three steps, with the step that needs each one.
REQUIRED_MODULES = [
    ("numpy", "shared"),
    ("h5py", "shared"),
    ("scipy", "shared"),
    ("cfgrib", "shared (ERA5)"),
    ("cv2", "TES"),
    ("netCDF4", "STIC"),
    ("pyhdf.SD", "STIC (MODIS albedo)"),
    ("pandas", "GPP"),
    ("xarray", "GPP"),
    ("rioxarray", "GPP"),
    ("rasterio", "GPP, COG"),
    ("shapely", "GPP"),
    ("geopandas", "GPP"),
    ("pyresample", "COG"),
    ("dotenv", "orchestrator"),
]

# Directories the steps read from, by the environment variable that names them.
DATA_DIRS = [
    "INPUT_RAD_DIR", "GEO_DIR", "CLOUD_DIR", "ERA5_DIR", "FCOVER_DIR",
    "ALBEDO_MOTA_DIR", "LULC_DIR", "PARH_DIR", "LAI_DIR", "OCO2_DIR",
    "GLC30_DIR",
]

OUTPUT_DIRS = ["OUTPUT_TES", "OUTPUT_STIC", "OUTPUT_GPP"]

CREDENTIALS = [
    ("EARTHDATA_USERNAME", "ECOSTRESS, MODIS albedo"),
    ("CDS_API_KEY", "ERA5"),
    ("CDSE_USERNAME", "FCOVER, LAI"),
    ("EUMETSAT_CONSUMER_KEY", "PAR"),
    ("S3_ENDPOINT_URL", "optional pre-staged data"),
]

# Products are float32 with fill -9999, or integers with fill 0.
FILL_FLOAT = -9999.0


class Report:
    """Collected check results, rendered as one line per check."""

    def __init__(self):
        self.rows = []

    def add(self, ok, label, detail=""):
        self.rows.append((ok, label, detail))
        log.info("  [%s] %-38s %s", "PASS" if ok else "FAIL", label, detail)
        return ok

    def note(self, label, detail=""):
        log.info("  [ -- ] %-38s %s", label, detail)

    def failures(self):
        return [label for ok, label, _ in self.rows if not ok]


def check_environment(report):
    log.info("Environment")

    v = sys.version_info
    report.add(
        (v.major, v.minor) >= MIN_PYTHON,
        "Python >= %d.%d" % MIN_PYTHON,
        "found %d.%d.%d" % (v.major, v.minor, v.micro),
    )

    for name, needed_by in REQUIRED_MODULES:
        try:
            importlib.import_module(name)
            report.add(True, name, needed_by)
        except Exception as exc:
            report.add(False, name, "%s — needed by %s" % (exc, needed_by))

    installdir = os.environ.get("RTTOV_INSTALLDIR", "")
    if not report.add(bool(installdir), "RTTOV_INSTALLDIR set", installdir or "unset"):
        return
    report.add(Path(installdir).is_dir(), "RTTOV installation", installdir)

    coef = Path(installdir) / "rtcoef_rttov12" / "rttov8pred54L" / "rtcoef_iss_1_ecostres.dat"
    report.add(coef.is_file(), "RTTOV ECOSTRESS coefficients", str(coef))

    wrapper = Path(installdir) / "wrapper"
    sys.path.append(str(wrapper))
    try:
        importlib.import_module("pyrttov")
        report.add(True, "pyrttov", str(wrapper))
    except Exception as exc:
        report.add(False, "pyrttov", "%s — TES cannot run" % exc)

    lut = os.environ.get("LUT_FILE", "")
    report.add(bool(lut) and Path(lut).is_file(), "GPP look-up table", lut or "LUT_FILE unset")


def check_directories(report):
    log.info("Directories")

    for var in DATA_DIRS:
        path = os.environ.get(var)
        if not path:
            report.note(var, "unset — falls back to the code default")
            continue
        p = Path(path)
        if p.is_dir():
            count = sum(1 for _ in p.iterdir())
            report.note(var, "%s (%d entr%s)" % (path, count, "y" if count == 1 else "ies"))
        else:
            report.note(var, "%s — absent, the downloader creates it" % path)

    for var in OUTPUT_DIRS:
        path = os.environ.get(var)
        if not path:
            report.note(var, "unset")
            continue
        p = Path(path)
        # Writability is worth a hard check: a read-only bind mount only fails
        # at the end of a run, after hours of processing.
        try:
            p.mkdir(parents=True, exist_ok=True)
            probe = p / ".selftest"
            probe.write_bytes(b"")
            probe.unlink()
            report.add(True, var + " writable", path)
        except Exception as exc:
            report.add(False, var + " writable", "%s — %s" % (path, exc))


def check_credentials(report):
    log.info("Credentials (presence only, never values)")
    for var, used_for in CREDENTIALS:
        state = "set" if os.environ.get(var) else "unset — downloads for %s will fail" % used_for
        report.note(var, state)


def _sig(value):
    # Significant digits, not decimals: GPP_LUE means are around 1e-5, and
    # rounding those to six decimals would leave a 4% error for the relative
    # tolerance to trip over.
    return float("%.9g" % value)


def dataset_stats(dset):
    data = dset[...]
    if np.issubdtype(data.dtype, np.floating):
        valid = np.isfinite(data) & (data != FILL_FLOAT)
    else:
        valid = data != 0
    n_valid = int(valid.sum())
    stats = {
        "shape": list(data.shape),
        "dtype": str(data.dtype),
        "valid_fraction": round(n_valid / data.size, 6),
    }
    if n_valid:
        v = data[valid].astype(np.float64)
        stats.update(
            minimum=_sig(v.min()),
            maximum=_sig(v.max()),
            mean=_sig(v.mean()),
        )
    return stats


def product_stats(path):
    out = {}
    with h5py.File(path, "r") as f:
        def visit(name, obj):
            if isinstance(obj, h5py.Dataset):
                out[name] = dataset_stats(obj)
        f.visititems(visit)
    return out


def compare(report, name, expected, actual, rel_tol, frac_tol):
    if expected["shape"] != actual["shape"]:
        return report.add(False, name, "shape %s, expected %s" % (actual["shape"], expected["shape"]))
    if expected["dtype"] != actual["dtype"]:
        return report.add(False, name, "dtype %s, expected %s" % (actual["dtype"], expected["dtype"]))

    delta = abs(actual["valid_fraction"] - expected["valid_fraction"])
    if delta > frac_tol:
        return report.add(
            False, name,
            "valid fraction %.4f, expected %.4f" % (actual["valid_fraction"], expected["valid_fraction"]),
        )

    for key in ("minimum", "maximum", "mean"):
        if key not in expected:
            continue
        if key not in actual:
            return report.add(False, name, "no valid pixel, expected %s %.4f" % (key, expected[key]))
        scale = max(abs(expected[key]), 1e-6)
        if abs(actual[key] - expected[key]) / scale > rel_tol:
            return report.add(
                False, name,
                "%s %.4f, expected %.4f (tolerance %.2g relative)" % (key, actual[key], expected[key], rel_tol),
            )

    return report.add(True, name, "mean %.4f, %.1f%% valid" % (
        actual.get("mean", float("nan")), 100 * actual["valid_fraction"]))


def check_products(report, output_dir, write_reference, rel_tol, frac_tol):
    if write_reference:
        return write_reference_file(output_dir)

    if not REFERENCE_FILE.is_file():
        report.add(False, "reference case", "%s is missing" % REFERENCE_FILE)
        return

    reference = json.loads(REFERENCE_FILE.read_text(encoding="utf-8"))
    case = reference["case"]
    log.info("Products — reference case %s, orbit %s scene %s",
             case["date"], case["orbit"], case["scene"])

    for relative, expected in reference["products"].items():
        path = Path(output_dir) / relative
        if not path.is_file():
            report.add(False, relative, "missing — run the pipeline on the reference case first")
            continue
        actual = product_stats(path)
        log.info("%s", relative)
        for dataset, expected_stats in expected.items():
            if dataset not in actual:
                report.add(False, dataset, "dataset absent from the product")
                continue
            compare(report, dataset, expected_stats, actual[dataset], rel_tol, frac_tol)


def write_reference_file(output_dir):
    reference = json.loads(REFERENCE_FILE.read_text(encoding="utf-8"))
    for relative in reference["products"]:
        path = Path(output_dir) / relative
        if not path.is_file():
            log.error("%s is missing, reference not written", path)
            return
        reference["products"][relative] = product_stats(path)
        log.info("Read %s", relative)
    REFERENCE_FILE.write_text(json.dumps(reference, indent=2) + "\n", encoding="utf-8")
    log.info("Wrote %s", REFERENCE_FILE)


def main(argv=None):
    parser = argparse.ArgumentParser(description="Verify an EEH2 installation.")
    parser.add_argument(
        "--products", action="store_true",
        help="Also compare the pipeline outputs against the reference case.",
    )
    parser.add_argument(
        "--output-dir", default=os.environ.get("EEH2_OUTPUT_ROOT", "/output"), metavar="DIR",
        help="Root of the output tree holding tes/, stic/ and gpp/ (default: /output).",
    )
    parser.add_argument(
        "--write-reference", action="store_true",
        help="Recompute the shipped reference from the outputs in --output-dir.",
    )
    parser.add_argument(
        "--relative-tolerance", type=float, default=1e-3, metavar="R",
        help="Accepted relative deviation on minimum, maximum and mean (default: 1e-3).",
    )
    parser.add_argument(
        "--valid-fraction-tolerance", type=float, default=5e-3, metavar="F",
        help="Accepted absolute deviation on the fraction of valid pixels (default: 5e-3).",
    )
    args = parser.parse_args(argv)

    report = Report()
    if not args.write_reference:
        check_environment(report)
        check_directories(report)
        check_credentials(report)

    if args.products or args.write_reference:
        check_products(report, args.output_dir, args.write_reference,
                       args.relative_tolerance, args.valid_fraction_tolerance)

    if args.write_reference:
        return 0

    failures = report.failures()
    total = len(report.rows)
    if failures:
        log.error("%d of %d checks failed: %s", len(failures), total, ", ".join(failures))
        return 1
    log.info("All %d checks passed.", total)
    return 0


if __name__ == "__main__":
    sys.exit(main())
