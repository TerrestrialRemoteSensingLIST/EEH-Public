#!/usr/bin/env python3
"""
EEH2 Pipeline Orchestrator
===========================
Runs the full TES → STIC → GPP processing chain, with optional COG conversion.

Configuration is read from a **.env** file (see .env.example).
Each step can be enabled/disabled individually via RUN_TES, RUN_STIC, RUN_GPP.

Usage
-----
    uv run python run_pipeline.py                  # full pipeline, all files
    uv run python run_pipeline.py --orbit 36798    # single orbit only
    uv run python run_pipeline.py --orbit 36798 --cog  # + COG conversion
    uv run python run_pipeline.py --start-date 2026-08-01 --end-date 2026-08-03  # date range
    uv run python run_pipeline.py --start-date 2026-08-01 --end-date 2026-08-01 --orbit 45031 --no-s3
    uv run python run_pipeline.py --step stic      # single step, all files
    uv run python run_pipeline.py --step cog        # COG conversion only
    uv run python run_pipeline.py --env my.env     # custom env file
"""

from __future__ import annotations

import argparse
import glob
import os
import re
import subprocess
import time
import sys
from datetime import datetime, timedelta
from pathlib import Path

from dotenv import dotenv_values

# Resolved from this file, not the cwd, so the pipeline can be launched from
# anywhere (the Docker entrypoint happens to run in /app, but callers may not).
REPO_ROOT = Path(__file__).resolve().parent
DOWNLOAD_SCRIPT = str(REPO_ROOT / "tools" / "download_sample.py")


def _env_bool(val: str | None) -> bool:
    if val is None:
        return False
    return val.strip().lower() in ("1", "true", "yes")


def _optional(env: dict, key: str) -> list[str] | None:
    val = env.get(key, "")
    if val and val.strip():
        return val.strip()
    return None


def _run_step(name: str, script: str, argv: list[str]) -> None:
    print(f"\n{'='*60}")
    print(f"  EEH2 Pipeline — {name}")
    print(f"{'='*60}\n")
    cmd = [sys.executable, str(REPO_ROOT / script)] + argv
    print(f"  Command: {' '.join(cmd)}\n")
    result = subprocess.run(cmd)
    if result.returncode != 0:
        print(f"\n[ERROR] {name} exited with code {result.returncode}")
        sys.exit(result.returncode)
    print(f"\n[OK] {name} finished (exit code 0).\n")


def run_tes(env: dict, input_files: list[str] | None = None) -> None:
    argv = [
        "--geo-dir",   env["GEO_DIR"],
        "--era5-dir",  env["ERA5_DIR"],
        "--output-dir", env["OUTPUT_TES"],
        "--rttov-installdir", env["RTTOV_INSTALLDIR"],
        # Not the default '.' (= /app, which is not a mounted volume):
        # the traceback of every failed granule was lost on container exit.
        "--error-log-dir", env["OUTPUT_TES"],
    ]
    if input_files is not None:
        argv += ["--input-files"] + input_files
    else:
        argv += ["--input-dir", env["INPUT_RAD_DIR"]]
    wrapper = _optional(env, "RTTOV_WRAPPER_DIR")
    if wrapper:
        argv += ["--rttov-wrapper-dir", wrapper]
    coef = _optional(env, "RTTOV_COEF_FILE")
    if coef:
        argv += ["--rttov-coef-file", coef]
    preload = _optional(env, "RTTOV_LIB_PRELOAD")
    if preload:
        argv += ["--rttov-lib-preload", preload]

    _run_step("TES (Land Surface Temperature)", "TES/TES_main.py", argv)


def run_stic(env: dict, input_files: list[str] | None = None) -> None:
    argv = [
        "--geo-dir",        env["GEO_DIR"],
        "--cloud-dir",      env["CLOUD_DIR"],
        "--fcover-dir",     env["FCOVER_DIR"],
        "--albedo-mota-dir", env["ALBEDO_MOTA_DIR"],
        "--lulc-dir",       env["LULC_DIR"],
        "--era5-dir",       env["ERA5_DIR"],
        "--output-dir",     env["OUTPUT_STIC"],
        # Not the default '.' (= /app, which is not a mounted volume):
        # the traceback of every failed granule was lost on container exit.
        "--error-log-dir", env["OUTPUT_STIC"],
    ]
    if input_files is not None:
        argv += ["--input-files"] + input_files
    else:
        argv += ["--input-dir", env.get("INPUT_LSTE_DIR", env["OUTPUT_TES"])]
    _run_step("STIC (Evapotranspiration)", "STIC/STIC_main.py", argv)


def run_gpp(env: dict, input_files: list[str] | None = None) -> None:
    lut = env.get("LUT_FILE") or str(
        REPO_ROOT / "GPP" / "LookUpTable_LUE_HH_gsFULL-globe.csv")
    argv = [
        "--geo-dir",   env["GEO_DIR"],
        "--cloud-dir", env["CLOUD_DIR"],
        "--parh-dir",  env["PARH_DIR"],
        "--lai-dir",   env["LAI_DIR"],
        "--fvc-dir",   env["FVC_DIR"],
        "--era5-dir",  env["ERA5_DIR"],
        "--oco2-dir",  env["OCO2_DIR"],
        "--glc30-dir", env["GLC30_DIR"],
        "--lut-file",  lut,
        "--output-dir", env["OUTPUT_GPP"],
        # Not the default '.' (= /app, which is not a mounted volume):
        # the traceback of every failed granule was lost on container exit.
        "--error-log-dir", env["OUTPUT_GPP"],
    ]
    if input_files is not None:
        argv += ["--stic-files"] + input_files
    else:
        argv += ["--stic-dir", env.get("INPUT_STIC_DIR", env["OUTPUT_STIC"])]
    _run_step("GPP (Gross Primary Productivity)", "GPP/GPP_main.py", argv)


def run_cog(env: dict, input_files: list[str] | None = None,
            orbit: str | None = None, force: bool = False) -> None:
    from tools.h5_to_cog import convert_batch

    geo_dir = env["GEO_DIR"]
    cloud_dir = env.get("CLOUD_DIR", "")
    cog_base = env.get("OUTPUT_COG", "./output/cog")

    tes_dir = env.get("OUTPUT_TES", "./output/tes")
    stic_dir = env.get("OUTPUT_STIC", "./output/stic")
    gpp_dir = env.get("OUTPUT_GPP", "./output/gpp")

    pattern = f"*{orbit}*" if orbit else "*.h5"

    conversions = [
        ("TES",  tes_dir,  os.path.join(cog_base, "tes"),  cloud_dir, None),
        ("STIC", stic_dir, os.path.join(cog_base, "stic"), None,      None),
        ("GPP",  gpp_dir,  os.path.join(cog_base, "gpp"),  None,      None),
        ("WUE",  gpp_dir,  os.path.join(cog_base, "wue"),  None,      "WUE"),
    ]

    for product, src_dir, out_dir, cld_dir, force_product in conversions:
        h5_files = sorted(glob.glob(os.path.join(src_dir, pattern)))
        if not h5_files:
            print(f"[COG] No {product} H5 files found in {src_dir}")
            continue

        print(f"\n{'='*60}")
        print(f"  EEH2 Pipeline — COG conversion: {product} ({len(h5_files)} files)")
        print(f"{'='*60}\n")

        convert_batch(
            h5_files, geo_dir, out_dir,
            cloud_dir=cld_dir if cld_dir else None,
            force_product=force_product,
            force=force,
            workers=min(4, len(h5_files)),
        )

    print(f"\n[OK] COG conversion completed.\n")


CMR_URL = "https://cmr.earthdata.nasa.gov/search/granules.json"
RAD_CONCEPT = "C2076116385-LPCLOUD"


def _list_orbits(start_date: str, end_date: str) -> None:
    import requests

    dt_start = datetime.strptime(start_date, "%Y-%m-%d")
    dt_end = datetime.strptime(end_date, "%Y-%m-%d")
    dt_cur = dt_start
    while dt_cur <= dt_end:
        day_str = dt_cur.strftime("%Y-%m-%d")
        temporal = f"{day_str}T00:00:00Z,{day_str}T23:59:59Z"
        # Paginated: CMR returns one page per request, ordered by start time, so
        # a single call listed only the first 200 of a 283-granule day and hid
        # every orbit of the late evening -- including the ones an operator would
        # pick to exercise the ERA5 midnight rollover.
        entries = []
        search_after = None
        while True:
            headers = {"Accept": "application/json"}
            if search_after:
                headers["CMR-Search-After"] = search_after
            resp = requests.get(CMR_URL, params={
                "concept_id": RAD_CONCEPT,
                "temporal": temporal,
                "page_size": 500,
            }, headers=headers)
            resp.raise_for_status()
            page = resp.json().get("feed", {}).get("entry", [])
            if not page:
                break
            entries += page
            search_after = resp.headers.get("CMR-Search-After")
            if not search_after:
                break
        orbits: dict[str, int] = {}
        for e in entries:
            m = re.search(r'_(\d{5})_\d{3}_', e.get("title", ""))
            if m:
                orbits[m.group(1)] = orbits.get(m.group(1), 0) + 1
        print(f"\n  {day_str}: {len(entries)} granules, {len(orbits)} orbits")
        for orb in sorted(orbits):
            print(f"    {orb}  ({orbits[orb]} scenes)")
        dt_cur += timedelta(days=1)
    print()


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Run the EEH2 processing pipeline (TES → STIC → GPP).",
    )
    parser.add_argument(
        "--env", default=".env", metavar="FILE",
        help="Path to the .env configuration file (default: .env).",
    )
    parser.add_argument(
        "--step", choices=["tes", "stic", "gpp", "cog"],
        help="Run only one step instead of the full pipeline.",
    )
    parser.add_argument(
        "--cog", action="store_true",
        help="Convert H5 outputs to Cloud Optimized GeoTIFF (COG) after processing.",
    )
    parser.add_argument(
        "--orbit", metavar="NUMBER",
        help="Process only files matching this orbit number (e.g. 36798). "
             "If no matching RAD files exist locally, attempts S3 download.",
    )
    parser.add_argument(
        "--input-files", nargs="+", metavar="FILE",
        help="Process specific RAD files (passed to TES). Outputs are propagated to STIC and GPP.",
    )
    parser.add_argument(
        "--download", action="store_true",
        help="Download sample data before running the pipeline.",
    )
    parser.add_argument(
        "--force", action="store_true",
        help="Re-run steps even if output files already exist.",
    )
    parser.add_argument(
        "--start-date", metavar="YYYY-MM-DD",
        help="Start of date range to process (inclusive). Discovers and downloads "
             "all ECOSTRESS orbits within the range. Requires --end-date.",
    )
    parser.add_argument(
        "--end-date", metavar="YYYY-MM-DD",
        help="End of date range to process (inclusive). Requires --start-date.",
    )
    parser.add_argument(
        "--no-s3", action="store_true",
        help="Skip S3 sync in downloads (passed to tools/download_sample.py).",
    )
    parser.add_argument(
        "--list-orbits", action="store_true",
        help="List available ECOSTRESS orbits for the given date range, then exit. "
             "Requires --start-date/--end-date.",
    )
    args = parser.parse_args(argv)

    if (args.start_date is None) != (args.end_date is None):
        parser.error("--start-date and --end-date must be used together.")
    if args.list_orbits and not args.start_date:
        parser.error("--list-orbits requires --start-date and --end-date.")
    if args.start_date and args.input_files:
        parser.error("--start-date/--end-date and --input-files are mutually exclusive.")

    if args.list_orbits:
        _list_orbits(args.start_date, args.end_date)
        return

    env_path = Path(args.env)
    if env_path.exists():
        env = {**dotenv_values(env_path), **os.environ}
    else:
        env = dict(os.environ)

    if args.download:
        dl_cmd = [sys.executable, DOWNLOAD_SCRIPT,
                  "--output-dir", env.get("INPUT_RAD_DIR", "./data").rsplit("/", 1)[0],
                  "--env", args.env]
        if args.no_s3:
            dl_cmd.append("--no-s3")
        print("="*60)
        print("  EEH2 Pipeline — Downloading sample data")
        print("="*60 + "\n")
        result = subprocess.run(dl_cmd)
        if result.returncode == 2:
            print("\n[WARN] Some downloads failed after retries — continuing with available data")
        elif result.returncode != 0:
            print(f"\n[ERROR] download exited with code {result.returncode}")
            sys.exit(result.returncode)

    os.makedirs(env.get("OUTPUT_TES", ""), exist_ok=True)
    os.makedirs(env.get("OUTPUT_STIC", ""), exist_ok=True)
    os.makedirs(env.get("OUTPUT_GPP", ""), exist_ok=True)

    # ── Resolve input files ──────────────────────────────────────────
    rad_dir = env.get("INPUT_RAD_DIR", "./data/rad")
    output_root = str(Path(rad_dir).parent)
    rad_files = args.input_files

    if not rad_files and args.orbit and not args.start_date:
        rad_files = sorted(glob.glob(os.path.join(rad_dir, f"*{args.orbit}*")))
        if not rad_files:
            print(f"[INFO] No RAD files for orbit {args.orbit} in {rad_dir}, downloading ECOSTRESS data...")
            dl_eco = [sys.executable, DOWNLOAD_SCRIPT,
                      "--orbit", args.orbit,
                      "--skip-era5", "--skip-clms", "--skip-lulc",
                      "--skip-glc30", "--skip-par",
                      "--output-dir", output_root,
                      "--env", args.env]
            if args.no_s3:
                dl_eco.append("--no-s3")
            result = subprocess.run(dl_eco)
            if result.returncode == 2:
                print("[WARN] Some ECOSTRESS downloads failed (retries exhausted)")
            rad_files = sorted(glob.glob(os.path.join(rad_dir, f"*{args.orbit}*")))
            if not rad_files:
                print(f"[ERROR] Still no RAD files for orbit {args.orbit} after download.")
                sys.exit(1)
        print(f"[INFO] Processing orbit {args.orbit}: {len(rad_files)} RAD file(s)")

    # ── Date range: download all orbits for each date ───────────────
    if args.start_date:
        dt_start = datetime.strptime(args.start_date, "%Y-%m-%d")
        dt_end = datetime.strptime(args.end_date, "%Y-%m-%d")
        dt_cur = dt_start
        while dt_cur <= dt_end:
            day_str = dt_cur.strftime("%Y-%m-%d")
            print(f"\n{'='*60}")
            print(f"  EEH2 Pipeline — Downloading data for {day_str}")
            print(f"{'='*60}\n")
            dl_cmd = [sys.executable, DOWNLOAD_SCRIPT,
                      "--date", day_str,
                      "--output-dir", output_root,
                      "--env", args.env]
            if args.orbit:
                dl_cmd += ["--orbit", args.orbit]
            if args.no_s3:
                dl_cmd.append("--no-s3")
            result = subprocess.run(dl_cmd)
            if result.returncode == 2:
                print(f"[WARN] Some downloads failed for {day_str} (retries exhausted)")
            dt_cur += timedelta(days=1)
        if args.orbit:
            rad_files = sorted(glob.glob(os.path.join(rad_dir, f"*{args.orbit}*")))
            if not rad_files:
                print(f"[ERROR] No RAD files for orbit {args.orbit} after download.")
                sys.exit(1)
            print(f"[INFO] Processing orbit {args.orbit}: {len(rad_files)} RAD file(s)")
        else:
            rad_files = None

    # Ensure ancillary data is present for each date found in RAD files
    source_files = rad_files or glob.glob(os.path.join(rad_dir, "*.h5"))
    dates_seen = set()
    for f in source_files:
        m = re.search(r'_(\d{4})(\d{2})(\d{2})T', Path(f).name)
        if m:
            dates_seen.add(f"{m.group(1)}-{m.group(2)}-{m.group(3)}")
    for orbit_date in sorted(dates_seen):
        print(f"[INFO] Ensuring ancillary data for {orbit_date}...")
        dl_anc = [sys.executable, DOWNLOAD_SCRIPT,
                  "--date", orbit_date,
                  "--skip-ecostress",
                  "--output-dir", output_root,
                  "--env", args.env]
        # Pass the orbit so ERA5 only fetches the hours that orbit needs
        if args.orbit:
            dl_anc += ["--orbit", args.orbit]
        if args.no_s3:
            dl_anc.append("--no-s3")
        result = subprocess.run(dl_anc)
        if result.returncode == 2:
            print(f"[WARN] Some ancillary downloads failed for {orbit_date} "
                  "(retries exhausted)")

    steps = {
        "tes":  (run_tes,  "RUN_TES"),
        "stic": (run_stic, "RUN_STIC"),
        "gpp":  (run_gpp,  "RUN_GPP"),
    }

    # ── Output directories per step ────────────────────────────────
    step_output = {
        "tes":  env.get("OUTPUT_TES",  "./data/output/tes"),
        "stic": env.get("OUTPUT_STIC", "./data/output/stic"),
        "gpp":  env.get("OUTPUT_GPP",  "./data/output/gpp"),
    }

    def _matches_selection(name: str) -> bool:
        if args.orbit and f"_{args.orbit}_" not in name:
            return False
        if dates_seen:
            m = re.search(r"_(\d{4})(\d{2})(\d{2})T", name)
            if not m or f"{m.group(1)}-{m.group(2)}-{m.group(3)}" not in dates_seen:
                return False
        return True

    # Steps that ran but wrote nothing, so the closing banner cannot claim success.
    empty_steps: list[str] = []

    def _outputs_written_since(step_name: str, since: float) -> list[str]:
        """Outputs of a step written since `since` (a time.time() stamp).

        By mtime rather than by set difference against a listing taken before the
        run, so that a --force re-run overwriting its own outputs still counts as
        having produced them. The 2 s slack absorbs the coarse mtime granularity
        of a bind-mounted volume.
        """
        written = []
        for f in glob.glob(os.path.join(step_output[step_name], "*.h5")):
            try:
                if os.path.getmtime(f) >= since - 2:
                    written.append(f)
            except OSError:
                pass
        return sorted(written)

    def _report_step_output(step_name: str, since: float) -> None:
        """Say what the step actually produced, not merely that it exited 0.

        TES, STIC and GPP catch per-granule failures so they can carry on with the
        rest, and so exit 0 even when every granule failed -- which used to be
        reported as success with nothing to show for it.
        """
        written = _outputs_written_since(step_name, since)
        if written:
            print(f"[INFO] {step_name.upper()} wrote {len(written)} output file(s).")
        else:
            empty_steps.append(step_name.upper())
            print(f"[WARN] {step_name.upper()} exited cleanly but wrote no output: "
                  "every granule was skipped or failed. See the error log in "
                  f"{step_output[step_name]}")

    def _step_has_output(step_name: str) -> list[str] | None:
        """Existing outputs of a step that match the requested selection.

        Returns None when nothing narrows the selection, i.e. the legacy "process
        the whole directory" mode. Callers must distinguish that from [], which
        means a selection was applied and matched nothing: passing [] on as an
        input file list would make the step fall back to its whole input
        directory and silently reprocess unrelated orbits and dates.
        """
        out_dir = step_output[step_name]
        if not args.orbit and not dates_seen:
            return None
        files = sorted(glob.glob(os.path.join(out_dir, "*.h5")))
        return [f for f in files if _matches_selection(os.path.basename(f))]

    if args.step:
        if args.step == "cog":
            run_cog(env, orbit=args.orbit, force=args.force)
        else:
            fn, _ = steps[args.step]
            existing = _step_has_output(args.step)
            if existing and not args.force:
                print(f"[SKIP] {args.step.upper()} — {len(existing)} output file(s) already exist (use --force to re-run)")
            else:
                upstream = {"stic": "tes", "gpp": "stic"}.get(args.step)
                step_files = _step_has_output(upstream) if upstream else rad_files
                if step_files is not None and not step_files:
                    print(f"[SKIP] {args.step.upper()} — no matching input for this "
                          "selection, nothing to process")
                else:
                    started = time.time()
                    fn(env, input_files=step_files)
                    _report_step_output(args.step, started)
    else:
        tes_files = rad_files
        # Seeded from what is already on disk so that a disabled or skipped
        # upstream step still hands the selection down instead of letting the
        # next one loose on its whole input directory.
        stic_files = _step_has_output("tes")
        gpp_files = _step_has_output("stic")

        for name, (fn, env_key) in steps.items():
            if not _env_bool(env.get(env_key, "true")):
                print(f"[SKIP] {name.upper()} (disabled in {env_path})")
                continue

            step_files = {"tes": tes_files, "stic": stic_files,
                          "gpp": gpp_files}[name]
            if step_files is not None and not step_files:
                previous = {"tes": "the selection", "stic": "TES",
                            "gpp": "STIC"}[name]
                print(f"[SKIP] {name.upper()} — {previous} produced no matching "
                      "input for this selection, nothing to process")
                continue

            existing = _step_has_output(name)
            if existing and not args.force:
                print(f"[SKIP] {name.upper()} — {len(existing)} output file(s) already exist (use --force to re-run)")
                if name == "tes":
                    stic_files = existing
                elif name == "stic":
                    gpp_files = existing
                continue

            started = time.time()
            if name == "tes":
                fn(env, input_files=tes_files)
                stic_files = _step_has_output("tes")
            elif name == "stic":
                fn(env, input_files=stic_files)
                gpp_files = _step_has_output("stic")
            elif name == "gpp":
                fn(env, input_files=gpp_files)
            _report_step_output(name, started)

    if args.cog:
        run_cog(env, orbit=args.orbit, force=args.force)

    print("\n" + "="*60)
    if empty_steps:
        print("  EEH2 Pipeline — finished, but "
              f"{', '.join(empty_steps)} produced no output.")
    else:
        print("  EEH2 Pipeline — All steps completed.")
    print("="*60)


if __name__ == "__main__":
    main()
