#!/usr/bin/env python3
"""
Download a minimal sample dataset for the EEH2 pipeline.

Reference orbit: 36798_018_20250101T165942

Downloads:
  - ECOSTRESS L1B_RAD, L1B_GEO, L2_CLOUD  (NASA Earthdata)
  - MCD43C3 MODIS albedo                    (NASA Earthdata)
  - OCO-2 XCO2                              (NASA Earthdata / GES DISC)
  - ERA5 single-levels + pressure-levels    (Copernicus CDS)
  - FCOVER 300m, LAI 300m                   (Copernicus CDSE)
  - PAR half-hourly                         (EUMETSAT Data Store, optional)
  - LULC PROBAV LC100                       (Zenodo)
  - GLC30 land-cover tiles                  (Zenodo)

Prerequisites — set in .env (or as environment variables):

  EARTHDATA_USERNAME=<your Earthdata user>
  EARTHDATA_PASSWORD=<your Earthdata password>
  CDS_API_URL=https://cds.climate.copernicus.eu/api
  CDS_API_KEY=<your-uid>:<your-api-key>
  CDSE_USERNAME=<your Copernicus Data Space username>
  CDSE_PASSWORD=<your Copernicus Data Space password>
  EUMETSAT_CONSUMER_KEY=<your key>    (optional, for PAR)
  EUMETSAT_CONSUMER_SECRET=<your secret>

If EARTHDATA_* are not set, falls back to ~/.netrc.
If CDS_* are not set, falls back to ~/.cdsapirc.
If EUMETSAT_* are not set, PAR is skipped (manual download from CM SAF).

Usage:
    uv run python tools/download_sample.py [--output-dir ./data]
"""

from __future__ import annotations

import argparse
import calendar
import io
import logging
import os
import re
import sys
import tempfile
import time
import zipfile
from datetime import date, datetime, timedelta
from pathlib import Path

import requests

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Reference orbit
# ---------------------------------------------------------------------------
ORBIT = "36798"
SCENE = "018"
DATE = "2025-01-01"
TEMPORAL = f"{DATE}T16:00:00Z,{DATE}T17:59:59Z"

# NASA Earthdata CMR
CMR_URL = "https://cmr.earthdata.nasa.gov/search/granules.json"
GEO_CONCEPT = "C2076087338-LPCLOUD"
RAD_CONCEPT = "C2076116385-LPCLOUD"
CLOUD_CONCEPT = "C2076115306-LPCLOUD"
MOTA_CONCEPT = "C2532068039-LPCLOUD"

# ERA5 — the hours actually needed are derived from the ECOSTRESS acquisition
# times on disk (see _era5_required_times); these are the fallback when no
# granule is available to derive from, and match the 16:57 UTC sample orbit.
ERA5_SINGLE_HOURS = [11, 12, 13, 16, 17]
ERA5_PRESSURE_HOURS = [16, 17]
# STIC derives the daily max Ta from single-level t2m at these fixed hours
ERA5_TAMAX_HOURS = [11, 12, 13]
# The processing domain. ERA5 is only fetched over this box, so a granule whose
# footprint falls outside it has no atmospheric profile to be corrected with:
# the ERA5 crop degenerates to a single grid cell, the bilinear weights collapse
# to zero and RTTOV is handed a 0 K 2m temperature. Such granules are therefore
# excluded at download time. Override with EEH2_DOMAIN="N,W,S,E".
DOMAIN_DEFAULT = (75.0, -20.0, -35.0, 60.0)  # N, W, S, E


def processing_domain() -> tuple[float, float, float, float]:
    raw = os.environ.get("EEH2_DOMAIN", "")
    if raw:
        try:
            n, w, s, e = (float(v) for v in raw.split(","))
            return n, w, s, e
        except ValueError:
            logger.warning(f"cannot parse EEH2_DOMAIN={raw!r} — using default domain")
    return DOMAIN_DEFAULT


ERA5_AREA = list(processing_domain())  # N, W, S, E, as the CDS API expects

ERA5_SINGLE_VARS = [
    "10m_u_component_of_wind", "10m_v_component_of_wind",
    "2m_dewpoint_temperature", "2m_temperature", "skin_temperature",
    "surface_pressure", "toa_incident_solar_radiation",
    "surface_solar_radiation_downwards", "total_column_water_vapour",
    "total_sky_direct_solar_radiation_at_surface", "boundary_layer_height",
]
ERA5_PRESSURE_VARS = [
    "temperature", "specific_humidity",
    "ozone_mass_mixing_ratio", "relative_humidity",
]
ERA5_PRESSURE_LEVELS = [
    "1", "2", "3", "5", "7", "10", "20", "30", "50", "70",
    "100", "125", "150", "175", "200", "225", "250", "300",
    "350", "400", "450", "500", "550", "600", "650", "700",
    "750", "775", "800", "825", "850", "875", "900", "925",
    "950", "975", "1000",
]

# OCO-2 GEOS L3 CO2 daily (NASA GES DISC, same Earthdata auth)
# Coverage ends ~2022; the GPP code falls back to closest year with same DOY
OCO2_CONCEPT = "C2240248762-GES_DISC"
# OCO-2 GEOS L3 CO2 exists for these days only; read_oco2 falls back to the
# same month and day in the nearest year it covers.
OCO2_COVERAGE = (date(2018, 8, 9), date(2022, 2, 28))

# CDSE (Copernicus Data Space Ecosystem) — FCOVER, LAI
CDSE_ODATA_URL = "https://catalogue.dataspace.copernicus.eu/odata/v1/Products"
CDSE_TOKEN_URL = ("https://identity.dataspace.copernicus.eu/auth/realms/CDSE"
                  "/protocol/openid-connect/token")


# ── Earthdata helpers ──────────────────────────────────────────────────────

def _ensure_netrc():
    user = os.environ.get("EARTHDATA_USERNAME", "")
    pwd = os.environ.get("EARTHDATA_PASSWORD", "")
    if not user or not pwd:
        return
    netrc_path = Path.home() / ".netrc"
    if netrc_path.exists():
        content = netrc_path.read_text()
        if "urs.earthdata.nasa.gov" in content:
            return
    with open(netrc_path, "a") as f:
        f.write(f"machine urs.earthdata.nasa.gov login {user} password {pwd}\n")
    netrc_path.chmod(0o600)


def _earthdata_session() -> requests.Session:
    _ensure_netrc()
    s = requests.Session()
    s.headers.update({"Accept": "application/json"})
    return s


CMR_PAGE_SIZE = 500  # CMR caps a page at 2000; 500 keeps each response small


def _cmr_all(params: dict) -> list[dict]:
    """Every granule matching `params`, following CMR pagination to the end.

    One request returns one page, and CMR orders granules by start time, so a day
    busier than the page size silently loses its late hours. That is what hid
    orbit 44173 of 2026-04-20 (21:57-23:16 UTC) outright: 283 L1B_RAD granules
    that day, a single page of 200, and the orbit filter then matched nothing and
    reported "No L1B_RAD granules found" as if the orbit did not exist. The
    2025-01-01 reference orbit passes at 16:57 and fell inside the first page, so
    nothing ever looked wrong.

    CMR-Search-After is used instead of page_num: it is the documented way past
    the offset limit and it does not re-scan the pages already read.
    """
    entries: list[dict] = []
    search_after = None
    while True:
        headers = {"Accept": "application/json"}
        if search_after:
            headers["CMR-Search-After"] = search_after
        resp = requests.get(CMR_URL, params={**params, "page_size": CMR_PAGE_SIZE},
                            headers=headers)
        resp.raise_for_status()
        page = resp.json().get("feed", {}).get("entry", [])
        if not page:
            return entries
        entries += page
        search_after = resp.headers.get("CMR-Search-After")
        if not search_after:
            return entries


def _search_cmr(session: requests.Session, concept_id: str, temporal: str,
                orbit_key: str, bbox: tuple[float, float, float, float] | None = None
                ) -> list[dict]:
    """Search CMR, optionally restricted to granules intersecting bbox (N, W, S, E).

    The bbox is applied server-side so CMR matches against each granule's own
    footprint polygon rather than a crude corner test.
    """
    params: dict = {
        "concept_id": concept_id,
        "temporal": temporal,
    }
    if bbox:
        n, w, s, e = bbox
        params["bounding_box"] = f"{w},{s},{e},{n}"
    entries = _cmr_all(params)
    return [e for e in entries if orbit_key in e.get("title", "")]


MAX_RETRIES = 4
RETRY_BACKOFF = 5  # seconds; doubles each retry → 5, 10, 20, 40


def _download_file(session: requests.Session, url: str, dest: Path,
                   max_retries: int = MAX_RETRIES) -> bool:
    """Download a single file with retry on transient HTTP errors.

    Returns True on success, False if all retries exhausted.
    """
    if dest.exists():
        logger.info(f"  [skip] {dest.name}")
        return True
    for attempt in range(1, max_retries + 1):
        try:
            logger.info(f"  [download] {dest.name} ...")
            resp = session.get(url, allow_redirects=True, timeout=300)
            resp.raise_for_status()
            dest.parent.mkdir(parents=True, exist_ok=True)
            dest.write_bytes(resp.content)
            logger.info(f"  [ok] {dest.name} ({len(resp.content) / 1e6:.1f} MB)")
            return True
        except (requests.exceptions.HTTPError,
                requests.exceptions.ConnectionError,
                requests.exceptions.Timeout) as exc:
            status = getattr(getattr(exc, "response", None), "status_code", None)
            if status and status < 500 and status != 429:
                logger.info(f"  [FAIL] {dest.name} — HTTP {status} (not retryable)")
                return False
            wait = RETRY_BACKOFF * (2 ** (attempt - 1))
            logger.info(f"  [retry {attempt}/{max_retries}] {dest.name} — {exc} "
                  f"(waiting {wait}s)")
            time.sleep(wait)
    logger.info(f"  [FAIL] {dest.name} — gave up after {max_retries} retries")
    return False


# ── ECOSTRESS ──────────────────────────────────────────────────────────────

def _scene_key(title: str) -> str | None:
    """The {orbit}_{scene} key shared by the RAD, GEO and CLOUD granule names."""
    m = re.search(r"_(\d{5}_\d{3})_", title)
    return m.group(1) if m else None


def _domain_scene_keys(session: requests.Session, orbit_key: str,
                       bbox: tuple[float, float, float, float]) -> set[str]:
    """Scene keys whose footprint intersects bbox (N, W, S, E).

    Derived from the L1B_GEO collection: of the three ECOSTRESS collections it is
    the only one CMR indexes spatially (a bounding_box search on L1B_RAD matches
    nothing at all), so GEO decides which scenes are in the domain and RAD and
    CLOUD are then matched on the scene key they share with it.
    """
    entries = _search_cmr(session, GEO_CONCEPT, TEMPORAL, orbit_key, bbox=bbox)
    return {k for k in (_scene_key(e.get("title", "")) for e in entries) if k}


def download_ecostress(session: requests.Session, output_dir: Path,
                       orbit: str | None = None,
                       whole_orbit: bool = False) -> tuple[int, int]:
    """Returns (total_files, failed_count).

    By default only granules intersecting the processing domain are fetched: an
    ISS orbit circles the globe, and granules outside the domain have no ERA5
    profile to be atmospherically corrected with, so downloading them only for
    TES to reject them later wastes bandwidth and runtime. Pass whole_orbit=True
    to fetch every granule regardless.
    """
    orbit_key = orbit if orbit else ""
    products = [
        ("L1B_RAD",  RAD_CONCEPT,   "rad"),
        ("L1B_GEO",  GEO_CONCEPT,   "geo"),
        ("L2_CLOUD", CLOUD_CONCEPT, "cloud"),
    ]
    filter_msg = f"orbit {orbit_key}" if orbit_key else "all orbits"
    bbox = None if whole_orbit else processing_domain()
    if bbox:
        n, w, s, e = bbox
        logger.info(f"\n  Processing domain: lat [{s}, {n}], lon [{w}, {e}] — granules "
              "outside it are skipped (--whole-orbit to download them anyway)")
    domain_keys = _domain_scene_keys(session, orbit_key, bbox) if bbox else None
    if domain_keys is not None and not domain_keys:
        logger.warning("no granule intersects the processing domain on this date")
    total, failed = 0, 0
    for label, concept_id, subdir in products:
        logger.info(f"\n--- {label} ({filter_msg}) ---")
        entries = _search_cmr(session, concept_id, TEMPORAL, orbit_key)
        if domain_keys is not None:
            inside = [e for e in entries
                      if _scene_key(e.get("title", "")) in domain_keys]
            outside = len(entries) - len(inside)
            if outside:
                logger.info(f"  {outside} of {len(entries)} granule(s) outside the "
                      "processing domain — not downloaded")
            entries = inside
        if not entries:
            logger.warning(f"No {label} granules found ({filter_msg})")
            continue
        logger.info(f"  Found {len(entries)} {label} granule(s)")
        dest_dir = output_dir / subdir
        for entry in entries:
            href = entry["links"][0]["href"]
            fname = href.split("/")[-1]
            total += 1
            if not _download_file(session, href, dest_dir / fname):
                failed += 1
    return total, failed


def download_mota(session: requests.Session, output_dir: Path) -> None:
    logger.info("\n--- MCD43C3 (MOTA) ---")
    dt = datetime.strptime(DATE, "%Y-%m-%d")
    target_doy = dt.timetuple().tm_yday
    target_key = f"A{dt.year}{target_doy:03d}"
    # MCD43C3 has ~8-day production lag; widen window to find the right DOY
    dt_start = dt - timedelta(days=10)
    dt_end = dt + timedelta(days=10)
    mota_temporal = (f"{dt_start.strftime('%Y-%m-%d')}T00:00:00Z,"
                     f"{dt_end.strftime('%Y-%m-%d')}T23:59:59Z")
    entries = _cmr_all({
        "concept_id": MOTA_CONCEPT,
        "temporal": mota_temporal,
    })
    if not entries:
        logger.warning(f"No MOTA granules for {DATE}")
        return
    dest_dir = output_dir / "mota"
    # Prefer the granule matching the exact DOY
    exact = [e for e in entries if target_key in e.get("title", "")]
    to_download = exact[:1] if exact else entries[:1]
    logger.info(f"  Target DOY: {target_key}, found {len(entries)} granules, "
          f"exact match: {len(exact)}")
    for entry in to_download:
        href = entry["links"][0]["href"]
        fname = href.split("/")[-1]
        _download_file(session, href, dest_dir / fname)


# ── OCO-2 (GES DISC, Earthdata auth) ─────────────────────────────────────

def _oco2_target_date(target_date: str) -> date | None:
    """The OCO-2 day that read_oco2 will actually select for `target_date`.

    The product covers 2018-08-09 to 2022-02-28 only, so read_oco2 falls back to
    the same month and day in the nearest year available (see read_oco2 in
    GPP/ReadData_EEH_GPP_final_TH.py). The download window has to follow that
    rule: the fixed range picked for the 2025-01-01 reference orbit leaves every
    other month and day with no file at all, and GPP then dies inside
    wrap_to_ECOSTRESS on a None dataset instead of reporting a missing input.
    """
    ymd = target_date.replace("-", "")
    year, month, day = int(ymd[:4]), int(ymd[4:6]), int(ymd[6:8])
    first, last = OCO2_COVERAGE
    candidates = []
    for y in range(first.year, last.year + 1):
        try:
            cand = date(y, month, day)
        except ValueError:
            continue  # 29 February outside a leap year
        if first <= cand <= last:
            candidates.append(cand)
    if not candidates:
        return None
    return min(candidates, key=lambda d: (abs(d.year - year), d.year))


def _oco2_on_disk(dest_dir: Path, target_date: str) -> Path | None:
    """A file already present that read_oco2 would accept for `target_date`.

    Mirrors read_oco2's own test -- the exact date, else any year with the same
    month and day -- so the skip decision cannot disagree with the consumer.
    """
    ymd = target_date.replace("-", "")
    mmdd = ymd[4:8]
    for f in sorted(dest_dir.glob("*.nc4")):
        if ymd in f.name or f"{mmdd}_" in f.name:
            return f
    return None


def download_oco2(session: requests.Session, output_dir: Path) -> None:
    logger.info("\n--- OCO-2 GEOS L3 CO2 ---")
    dest_dir = output_dir / "oco2"
    dest_dir.mkdir(parents=True, exist_ok=True)

    # Date-aware: "any .nc4 present" let the 2021-12/2022-01 sample files stand
    # in for every date, so nothing was ever downloaded and read_oco2 returned
    # None for any month and day outside that fortnight.
    have = _oco2_on_disk(dest_dir, DATE)
    if have:
        logger.info(f"  [skip] {have.name} already covers {DATE}")
        return

    wanted = _oco2_target_date(DATE)
    if wanted is None:
        first, last = OCO2_COVERAGE
        logger.warning(f"{DATE} has no counterpart in OCO-2 coverage "
              f"({first} to {last}) - GPP will have no CO2 input")
        return
    if wanted.strftime("%Y%m%d") != DATE.replace("-", ""):
        logger.info(f"  Outside OCO-2 coverage; using the same day of {wanted.year}: "
              f"{wanted}")

    entries = _cmr_all({
        "concept_id": OCO2_CONCEPT,
        "temporal": f"{wanted}T00:00:00Z,{wanted}T23:59:59Z",
    })
    if not entries:
        logger.warning("No OCO-2 granules found")
        return
    logger.info(f"  Found {len(entries)} OCO-2 granules")
    for entry in entries:
        links = [lk["href"] for lk in entry.get("links", [])
                 if lk.get("href", "").endswith(".nc4")]
        if not links:
            continue
        href = links[0]
        fname = href.split("/")[-1]
        _download_file(session, href, dest_dir / fname)


# ── FCOVER + LAI (CDSE) ──────────────────────────────────────────────────

def _get_cdse_token() -> str | None:
    username = os.environ.get("CDSE_USERNAME", "")
    password = os.environ.get("CDSE_PASSWORD", "")
    if not username or not password:
        logger.info("  [skip] CDSE_USERNAME/CDSE_PASSWORD not set — skipping FCOVER/LAI")
        return None
    resp = requests.post(CDSE_TOKEN_URL, data={
        "grant_type": "password",
        "client_id": "cdse-public",
        "username": username,
        "password": password,
    })
    resp.raise_for_status()
    return resp.json()["access_token"]


def _clms_composite_date(name: str) -> str | None:
    """Extract YYYYMMDD composite date from CLMS product name."""
    m = re.search(r'_(\d{12})_', name)
    return m.group(1)[:8] if m else None


def _clms_rt_score(name: str) -> int:
    """Score CLMS product quality: final (no RT) = 100, RT6 = 6, ..., RT0 = 0."""
    rt_match = re.search(r'[-_]RT(\d)', name)
    return int(rt_match.group(1)) if rt_match else 100


def _clms_expected_date(target_date: str) -> str:
    """The composite date STIC and GPP will look for, as YYYYMMDD.

    CLMS 300m FCOVER/LAI are 10-day composites labelled by the end of the period:
    day 1-10 -> YYYYMM10, day 11-20 -> YYYYMM20, day 21+ -> YYYYMM{lastday}.
    STIC's Read_FVC indexes its product map with exactly this key, so a composite
    from a neighbouring period is of no use to it.
    """
    ymd = target_date.replace("-", "")
    year, month, day = int(ymd[:4]), int(ymd[4:6]), int(ymd[6:8])
    if day <= 10:
        sday = 10
    elif day <= 20:
        sday = 20
    else:
        sday = calendar.monthrange(year, month)[1]
    return f"{year:04d}{month:02d}{sday:02d}"


def _clms_required_dates(target_date: str) -> list[str]:
    """Every 10-day composite date needed on `target_date`, oldest first.

    Two consumers with two different needs. STIC indexes its product map with
    the composite that *contains* the date (_clms_expected_date). GPP instead
    interpolates LAI and FCOVER between the two composites that *bracket* it,
    so it needs a second, earlier one -- see interpolate_LAI_value in
    GPP/ReadData_EEH_GPP_final_TH.py, whose grid of observation days (the 10th,
    the 20th and the last day of each month) is the one mirrored here. On a
    composite boundary the two needs do not even overlap: 2026-08-10 sends STIC
    to 20260810 and GPP to 20260731 and 20260820.

    Fetching only STIC's composite is what made GPP raise StopIteration on
    every granule of orbit 45770, while the 2025-01-01 reference orbit passed:
    day-of-year 1 falls in GPP's <=10 special case, which reads a single file
    and never interpolates, so the bracketing path was never exercised.
    """
    ymd = target_date.replace("-", "")
    year, month, day = int(ymd[:4]), int(ymd[4:6]), int(ymd[6:8])
    target_doy = datetime(year, month, day).timetuple().tm_yday

    required = {_clms_expected_date(target_date)}

    if target_doy <= 10:
        required.add(f"{year:04d}0110")
    elif target_doy >= 364:
        # GPP takes the latest December composite it can find; that is 1231.
        required.add(f"{year:04d}1231")
    else:
        grid = []
        for m in range(1, 13):
            for d in (10, 20, calendar.monthrange(year, m)[1]):
                dt = datetime(year, m, d)
                grid.append((dt.timetuple().tm_yday, dt))
        before = max((item for item in grid if item[0] < target_doy),
                     default=None)
        after = min((item for item in grid if item[0] > target_doy),
                    default=None)
        for item in (before, after):
            if item is not None:
                dt = item[1]
                required.add(f"{dt.year:04d}{dt.month:02d}{dt.day:02d}")

    return sorted(required)


def _clms_missing(output_dir: Path, target_date: str) -> list[str]:
    """Which required composites are absent from fcover/ and lai/, as "sub/date".

    Used to decide whether CDSE still has work to do. The S3 sync may well have
    delivered one of GPP's two brackets and not the other, so "the sync touched
    this directory" is not evidence that the date is covered.
    """
    missing = []
    for subdir in ("fcover", "lai"):
        have = {_clms_composite_date(f.name)
                for f in (output_dir / subdir).glob("*.nc")}
        missing += [f"{subdir}/{d}" for d in _clms_required_dates(target_date)
                    if d not in have]
    return missing


def _download_clms_product(token: str, product: dict, dest_dir: Path) -> None:
    """Fetch one CLMS product into `dest_dir`, unwrapping its ZIP if there is one."""
    product_name = product["Name"]
    filename = (product_name[:-3] + ".nc"
                if product_name.endswith("_nc") else product_name)
    score = _clms_rt_score(product_name)
    rt_label = f"RT{score}" if score < 100 else "final"
    logger.info(f"  [download] {filename} ({rt_label}) ...")

    dl_url = (f"https://zipper.dataspace.copernicus.eu/odata/v1"
              f"/Products({product['Id']})/$value")
    dl_resp = requests.get(
        dl_url, headers={"Authorization": f"Bearer {token}"}, stream=True,
    )
    dl_resp.raise_for_status()
    content = dl_resp.content

    try:
        with zipfile.ZipFile(io.BytesIO(content)) as zf:
            nc_files = [n for n in zf.namelist() if n.endswith(".nc")]
            if not nc_files:
                logger.warning("no .nc inside the archive")
                return
            data = zf.read(nc_files[0])
    except zipfile.BadZipFile:
        data = content

    (dest_dir / filename).write_bytes(data)
    logger.info(f"  [ok] {filename} ({len(data) / 1e6:.1f} MB)")


def download_fcover_lai(output_dir: Path) -> None:
    token = _get_cdse_token()
    if not token:
        return

    required = _clms_required_dates(DATE)

    for product_key, subdir in [("FCOVER300", "fcover"), ("LAI300", "lai")]:
        logger.info(f"\n--- {product_key} (CDSE) ---")
        dest_dir = output_dir / subdir
        dest_dir.mkdir(parents=True, exist_ok=True)

        # Skip per composite, not per directory: "any .nc present" made the
        # 2025-01 sample files suppress every later download, and STIC then died
        # with KeyError on the missing dekad.
        on_disk: dict[str, str] = {}
        for f in sorted(dest_dir.glob("*.nc")):
            cdate = _clms_composite_date(f.name)
            if cdate:
                on_disk.setdefault(cdate, f.name)

        missing = []
        for want in required:
            if want in on_disk:
                logger.info(f"  [skip] composite {want} already in {subdir}/ "
                      f"({on_disk[want]})")
            else:
                missing.append(want)
        if not missing:
            continue

        # Bound the query by the newest composite wanted, not by DATE: GPP's
        # later bracket can start after the target date, and 'le DATE' dropped
        # it without a word.
        upper = max(missing)
        upper_iso = f"{upper[:4]}-{upper[4:6]}-{upper[6:8]}"
        filter_str = (
            "Collection/Name eq 'CLMS' "
            f"and contains(Name,'{product_key}') "
            "and endswith(Name,'_nc') "
            f"and ContentDate/Start le {upper_iso}T23:59:59.999Z"
        )
        resp = requests.get(CDSE_ODATA_URL, params={
            "$filter": filter_str,
            "$orderby": "ContentDate/Start desc",
            "$top": 100,
        })
        resp.raise_for_status()
        products = resp.json().get("value", [])
        if not products:
            logger.warning(f"no {product_key} products found on CDSE")
            continue

        by_date: dict[str, list[dict]] = {}
        for p in products:
            cdate = _clms_composite_date(p.get("Name", ""))
            if cdate:
                by_date.setdefault(cdate, []).append(p)

        for want in missing:
            if want not in by_date:
                # No neighbouring composite is substituted: GPP selects by name
                # and would only raise StopIteration on it later, further from
                # the cause.
                logger.warning(f"composite {want} is not published by CLMS yet "
                      f"— {product_key} stays incomplete for {DATE}")
                continue
            best = max(by_date[want],
                       key=lambda p: _clms_rt_score(p.get("Name", "")))
            _download_clms_product(token, best, dest_dir)


# ── PAR (CM SAF SARAH, via EUMETSAT Data Store) ──────────────────────────

EUMETSAT_SARAH_COLLECTION = "EO:EUM:DAT:0863"

# GPP reads the half-hourly instantaneous product only: it takes the overpass
# time slice and derives the daily mean from the same 48 time steps (see
# extract_PAR_from_global_data in GPP/ReadData_EEH_GPP_final_TH.py). The daily
# (PARdm) and monthly (PARmm) means are never opened, so they are not fetched.
PAR_PRODUCT_PREFIX = "PARin"
# Variables extract_PAR_from_global_data needs; a file missing any of them
# would fail mid-granule instead of at download time.
PAR_REQUIRED_VARS = ("PAR", "lat", "lon", "time")


def _is_valid_par_netcdf(path: Path) -> bool:
    """Whether `path` is a NetCDF that GPP can actually open.

    The extension alone proves nothing: the EUMETSAT Data Store returns each
    product as a ZIP, which used to be written straight to a .nc name, and GPP
    then died on "NetCDF: Unknown file format" one granule at a time.
    """
    if not path.is_file() or zipfile.is_zipfile(path):
        return False
    try:
        import netCDF4
    except ImportError:
        # Cannot validate here; the archive check above still rules out the
        # case that actually occurred.
        return True
    try:
        with netCDF4.Dataset(path) as ds:
            return all(v in ds.variables for v in PAR_REQUIRED_VARS)
    except OSError:
        return False


def _par_member(zf: zipfile.ZipFile) -> zipfile.ZipInfo | None:
    """The PARin NetCDF inside a Data Store archive, if there is one.

    The archive also holds EOPMetadata.xml and manifest.xml, and a PARdm/PARmm
    archive holds a NetCDF that GPP would never select -- so match on the
    product prefix, not on ".nc" alone. Members carrying a path are ignored:
    nothing legitimate needs one.
    """
    for info in zf.infolist():
        if info.is_dir():
            continue
        name = info.filename
        if name != Path(name).name:
            continue
        if name.startswith(PAR_PRODUCT_PREFIX) and name.endswith(".nc"):
            return info
    return None


def _extract_par_archive(zf: zipfile.ZipFile, dest_dir: Path) -> Path | None:
    """Extract the PARin NetCDF of `zf` into `dest_dir`, validated.

    Written to a temporary file and moved into place only once it opens as a
    NetCDF, for two reasons: the member name equals the name of the archive it
    comes from, so a direct write would truncate the file being read; and a
    half-written file must never be left under a name the skip check trusts.
    """
    info = _par_member(zf)
    if info is None:
        return None
    dest = dest_dir / info.filename
    tmp = None
    try:
        with zf.open(info) as src:
            fd, tmp_name = tempfile.mkstemp(dir=dest_dir, suffix=".nc.part")
            tmp = Path(tmp_name)
            with os.fdopen(fd, "wb") as out:
                while True:
                    chunk = src.read(8 << 20)
                    if not chunk:
                        break
                    out.write(chunk)
        if not _is_valid_par_netcdf(tmp):
            logger.warning(f"{info.filename} in the archive is not a readable "
                  "NetCDF")
            return None
        tmp.replace(dest)
        tmp = None
        return dest
    finally:
        if tmp is not None and tmp.exists():
            tmp.unlink()


def repair_par_archives(output_dir: Path) -> None:
    """Turn any PARin*.nc that is really a ZIP into the NetCDF it contains.

    Runs whatever populated parh/ -- an EUMETSAT download, an S3 sync or a
    manual copy -- because the bad file is indistinguishable from a good one by
    name, and GPP selects it by name. Files that are already NetCDF are left
    strictly alone, so correct data from earlier runs is never rewritten.
    """
    dest_dir = output_dir / "parh"
    if not dest_dir.is_dir():
        return
    for nc in sorted(dest_dir.glob("*.nc")):
        if not zipfile.is_zipfile(nc):
            continue
        try:
            with zipfile.ZipFile(nc) as zf:
                extracted = _extract_par_archive(zf, dest_dir)
        except (OSError, zipfile.BadZipFile) as err:
            logger.warning(f"cannot read {nc.name} as an archive: {err}")
            continue
        if extracted is None:
            # No usable PARin inside (a PARdm/PARmm archive, or a damaged one).
            # Left as *.nc it would keep satisfying the skip check and keep
            # failing GPP, so move it out of the way.
            quarantined = nc.with_suffix(".nc.zip")
            nc.replace(quarantined)
            logger.warning(f"{nc.name} holds no usable {PAR_PRODUCT_PREFIX} "
                  f"NetCDF — set aside as {quarantined.name}")
            continue
        size_mb = extracted.stat().st_size / 1e6
        if extracted.resolve() != nc.resolve():
            nc.unlink()
        logger.info(f"  [fix] {nc.name} was a ZIP — extracted {extracted.name} "
              f"({size_mb:.1f} MB)")


def download_par(output_dir: Path) -> bool:
    """Download half-hourly PAR (PARin) from EUMETSAT Data Store via eumdac.
    Returns True if PAR files are available after this call."""
    logger.info("\n--- PAR half-hourly (CM SAF SARAH) ---")
    dest_dir = output_dir / "parh"
    dest_dir.mkdir(parents=True, exist_ok=True)
    repair_par_archives(output_dir)
    # PAR is per-day and GPP selects it by name ("PARin" + YYYYMMDD), so the
    # presence of some other day's file says nothing about this date. The file
    # also has to be readable: a name match alone used to let a ZIP stand in
    # for the NetCDF and fail every granule downstream.
    ymd = DATE.replace("-", "")
    if any(_is_valid_par_netcdf(f)
           for f in dest_dir.glob(f"{PAR_PRODUCT_PREFIX}{ymd}*.nc")):
        logger.info(f"  [skip] {PAR_PRODUCT_PREFIX} files for {DATE} already present")
        return True

    consumer_key = os.environ.get("EUMETSAT_CONSUMER_KEY", "")
    consumer_secret = os.environ.get("EUMETSAT_CONSUMER_SECRET", "")
    if not consumer_key or not consumer_secret:
        logger.info("  [skip] EUMETSAT_CONSUMER_KEY/SECRET not set")
        return False

    try:
        import eumdac
    except ImportError:
        logger.info("  [skip] eumdac not installed — run: pip install eumdac")
        return False

    token = eumdac.AccessToken((consumer_key, consumer_secret))
    datastore = eumdac.DataStore(token)

    try:
        collection = datastore.get_collection(EUMETSAT_SARAH_COLLECTION)
    except Exception as e:
        logger.warning(f"Cannot access collection {EUMETSAT_SARAH_COLLECTION}: {e}")
        return False

    logger.info(f"  Searching collection {EUMETSAT_SARAH_COLLECTION} ...")
    products = list(collection.search(dtstart=DATE, dtend=DATE))
    logger.info(f"  Found {len(products)} product(s) in SARAH-3 for {DATE}")

    par_products = [p for p in products if PAR_PRODUCT_PREFIX in str(p)]
    if not par_products:
        titles = [str(p) for p in products[:5]]
        logger.warning(f"No {PAR_PRODUCT_PREFIX} products found "
              f"(sample: {titles})")
        return False

    logger.info(f"  {len(par_products)} {PAR_PRODUCT_PREFIX} product(s) to download")
    count = 0
    for product in par_products:
        title = str(product)
        try:
            with product.open() as stream:
                data = stream.read() if hasattr(stream, 'read') else b''.join(stream)
            if not data:
                continue
            # The Data Store answers with an archive (NetCDF + two XML metadata
            # files). Only its PARin member is kept, under the member's own
            # name, which is the CM SAF name GPP matches on.
            try:
                with zipfile.ZipFile(io.BytesIO(data)) as zf:
                    dest = _extract_par_archive(zf, dest_dir)
            except zipfile.BadZipFile:
                # A bare NetCDF is acceptable too, but only once it is shown to
                # be one -- never by trusting the product title.
                dest = _write_par_netcdf(dest_dir, title, data)
            if dest is None:
                logger.warning(f"no usable {PAR_PRODUCT_PREFIX} NetCDF in {title}")
                continue
            count += 1
            logger.info(f"  [ok] {dest.name} ({dest.stat().st_size / 1e6:.1f} MB)")
        except Exception as e:
            logger.warning(f"Failed to download {title}: {e}")

    if count == 0:
        logger.warning(f"No {PAR_PRODUCT_PREFIX} .nc files extracted")
        return False
    logger.info(f"  [ok] Downloaded {count} {PAR_PRODUCT_PREFIX} file(s)")
    return True


def _write_par_netcdf(dest_dir: Path, title: str, data: bytes) -> Path | None:
    """Store an unarchived response, but only if it really is a NetCDF."""
    fd, tmp_name = tempfile.mkstemp(dir=dest_dir, suffix=".nc.part")
    tmp = Path(tmp_name)
    with os.fdopen(fd, "wb") as out:
        out.write(data)
    if not _is_valid_par_netcdf(tmp):
        tmp.unlink()
        logger.warning(f"{title} is neither an archive nor a readable NetCDF")
        return None
    fname = title.replace("/", "_").replace(chr(92), "_")
    if not fname.endswith(".nc"):
        fname += ".nc"
    dest = dest_dir / fname
    tmp.replace(dest)
    return dest


# ── LULC (Copernicus PROBAV LC100, hosted on Zenodo) ─────────────────────

LULC_FILENAME = "PROBAV_LC100_global_v3.0.1_2018-conso_Discrete-Classification-map_EPSG-4326.tif"
LULC_URL = (
    "https://zenodo.org/api/records/3518038/files/"
    "PROBAV_LC100_global_v3.0.1_2018-conso_Discrete-Classification-map_EPSG-4326.tif/content"
)


def download_lulc(output_dir: Path) -> None:
    logger.info("\n--- LULC (PROBAV LC100) ---")
    dest_dir = output_dir / "lulc"
    dest_dir.mkdir(parents=True, exist_ok=True)
    dest = dest_dir / LULC_FILENAME
    if dest.exists():
        logger.info(f"  [skip] {LULC_FILENAME}")
        return
    logger.info(f"  [download] {LULC_FILENAME} (~1.6 GB) ...")
    resp = requests.get(LULC_URL, stream=True)
    resp.raise_for_status()
    with open(dest, "wb") as f:
        for chunk in resp.iter_content(chunk_size=8192 * 1024):
            f.write(chunk)
    logger.info(f"  [ok] {LULC_FILENAME} ({dest.stat().st_size / 1e9:.2f} GB)")


# ── GLC30 / GLC_FCS30D (Zenodo, tiles by longitude strip) ────────────────

GLC30_ZENODO_RECORD = "8239305"
GLC30_API_URL = f"https://zenodo.org/api/records/{GLC30_ZENODO_RECORD}"


# Stride for the fallback pixel scan, used only when a granule has no bounding
# metadata: the box only has to be good to a fraction of a degree to select
# 5-degree tiles. GEO_BBOX_MARGIN_DEG then pads each footprint, covering both the
# stride error and reprojection edge effects -- an extra tile costs a little
# bandwidth, a missing one leaves a hole in the land cover.
GEO_BBOX_STRIDE = 25
GEO_BBOX_MARGIN_DEG = 0.25
# ERA5 is on a 0.25 deg grid and bilinearly interpolated, so a granule needs at
# least one full cell inside the domain in each direction to be usable. Keep in
# sync with DOMAIN_MIN_OVERLAP_DEG in TES/TES_main.py.
DOMAIN_MIN_OVERLAP_DEG = 0.25


def _geo_granule_box(geo_file: Path) -> tuple | None:
    """One granule's footprint as (lat_min, lat_max, lon_min, lon_max).

    Taken from StandardMetadata, which ECOSTRESS GEO products carry, so no pixel
    data is touched: the geolocation arrays are ~100 MB per granule and an HDF5
    hyperslab still reads whole chunks, which made scanning a swath take minutes.
    Falls back to a strided scan of the arrays if the metadata is absent.
    """
    import h5py
    with h5py.File(geo_file, "r") as f:
        meta = f.get("StandardMetadata")
        if meta is not None and "NorthBoundingCoordinate" in meta:
            return (float(meta["SouthBoundingCoordinate"][()]),
                    float(meta["NorthBoundingCoordinate"][()]),
                    float(meta["WestBoundingCoordinate"][()]),
                    float(meta["EastBoundingCoordinate"][()]))
        st = GEO_BBOX_STRIDE
        lat = f["Geolocation/latitude"][::st, ::st]
        lon = f["Geolocation/longitude"][::st, ::st]
    valid = (lat > -900) & (lon > -900)
    if not valid.any():
        return None
    return (float(lat[valid].min()), float(lat[valid].max()),
            float(lon[valid].min()), float(lon[valid].max()))


GLC30_FETCHED_FILE = ".glc30_fetched_zips"


def _glc30_fetched_zips(dest_dir: Path) -> set[str]:
    """Names of the GLC30 ZIPs already downloaded and extracted here.

    Recorded because the tiles on disk cannot answer the question on their own:
    GLC30 has no tile over open water, so a longitude strip the swath crosses at
    sea would look permanently missing and its multi-GB ZIP would be fetched again
    on every run.
    """
    marker = dest_dir / GLC30_FETCHED_FILE
    if not marker.exists():
        return set()
    return {line.strip() for line in
            marker.read_text(encoding="utf-8").splitlines() if line.strip()}


def _record_glc30_zip(dest_dir: Path, zip_name: str) -> None:
    with (dest_dir / GLC30_FETCHED_FILE).open("a", encoding="utf-8") as f:
        f.write(zip_name + chr(10))


def _union_boxes_in_domain(boxes, label: str) -> tuple | None:
    """Union of granule footprints clipped to the domain, or None.

    `boxes` is an iterable of (lat_min, lat_max, lon_min, lon_max), each padded by
    GEO_BBOX_MARGIN_DEG and clipped to the processing domain before being unioned.
    Unioning first and clipping once at the end returned the entire domain -- an
    ISS orbit circles the globe -- so land cover was requested for all of Europe
    and Africa when four granules over Brittany were all that TES would accept.
    """
    n, w, s, e = processing_domain()
    m = GEO_BBOX_MARGIN_DEG
    lat_lo = lon_lo = float("inf")
    lat_hi = lon_hi = float("-inf")
    total = kept = outside = 0
    for box in boxes:
        total += 1
        if box is None:
            continue
        g_lat_lo, g_lat_hi, g_lon_lo, g_lon_hi = box
        if g_lon_lo > g_lon_hi:
            # Crosses the antimeridian. A CDS area is N/W/S/E with W < E, so the
            # processing domain can never span it either: such a granule sits at
            # |lon| ~ 180, outside any valid domain, and must not be allowed to
            # stretch the box across the globe.
            outside += 1
            continue
        g_lat_lo, g_lat_hi = g_lat_lo - m, g_lat_hi + m
        g_lon_lo, g_lon_hi = g_lon_lo - m, g_lon_hi + m
        c_lat_lo, c_lat_hi = max(g_lat_lo, s), min(g_lat_hi, n)
        c_lon_lo, c_lon_hi = max(g_lon_lo, w), min(g_lon_hi, e)
        if (c_lat_hi - c_lat_lo < DOMAIN_MIN_OVERLAP_DEG
                or c_lon_hi - c_lon_lo < DOMAIN_MIN_OVERLAP_DEG):
            outside += 1
            continue
        kept += 1
        lat_lo, lat_hi = min(lat_lo, c_lat_lo), max(lat_hi, c_lat_hi)
        lon_lo, lon_hi = min(lon_lo, c_lon_lo), max(lon_hi, c_lon_hi)

    if total == 0:
        logger.warning(f"no {label} to build the footprint from")
        return None
    if not kept:
        logger.warning(f"none of the {total} {label} intersect "
              "the processing domain")
        return None
    if outside:
        logger.info(f"  footprint from {kept} of {total} {label} "
              f"({outside} outside the processing domain)")
    return (lat_lo, lat_hi, lon_lo, lon_hi)


def _cmr_geo_boxes(orbit: str | None, date: str) -> list[tuple]:
    """GEO footprints for `orbit` on `date`, from CMR metadata alone.

    The S3 sync runs before any granule is downloaded, so on a first run the
    footprint cannot come from geo/*.h5: there are none yet. That silently
    disabled the GLC30 tile filter and made the sync fall back to the global tile
    set -- 43 GB of central and eastern Asian land cover fetched for an orbit over
    Europe, and still counting when it was stopped. The reference orbit never
    showed it because its GEO granules were already on disk.

    CMR carries a bounding box per granule, so the question is answerable before
    downloading anything. Only L1B_GEO is spatially indexed in CMR, which is the
    collection wanted here anyway.

    The two sources agree: on the reference orbit both resolve to
    lat [-35.00, -6.82] lon [19.43, 48.50] and select the same 58 tiles. CMR
    reports the same StandardMetadata bounds rounded to 7 decimals, so its box can
    be up to 5e-8 deg smaller -- six orders of magnitude inside
    GEO_BBOX_MARGIN_DEG, which is there to absorb exactly this.
    """
    n, w, s, e = processing_domain()
    params = {
        "concept_id": GEO_CONCEPT,
        "temporal": f"{date}T00:00:00Z,{date}T23:59:59Z",
        "bounding_box": f"{w},{s},{e},{n}",
    }
    try:
        entries = _cmr_all(params)
    except requests.RequestException as err:
        logger.warning(f"CMR footprint lookup failed: {err}")
        return []

    boxes = []
    for entry in entries:
        if orbit and f"_{orbit}_" not in entry.get("title", ""):
            continue
        for raw in entry.get("boxes") or []:
            parts = raw.split()
            if len(parts) != 4:
                continue
            # CMR orders a box "south west north east".
            b_s, b_w, b_n, b_e = (float(v) for v in parts)
            boxes.append((b_s, b_n, b_w, b_e))
            break
    return boxes


def _get_orbit_bbox(output_dir: Path, orbit: str | None = None,
                    date: str | None = None) -> tuple | None:
    """Footprint of the orbit GEO granules: (lat_min, lat_max, lon_min, lon_max).

    Restricted to the granules of `orbit` / `date` when given, unioned over all of
    them, and clipped to the processing domain -- granules outside it are never
    processed, so no land cover is needed there.

    Read from geo/*.h5 when those are on disk and from CMR otherwise, so the
    answer no longer depends on whether the ECOSTRESS download has already run.

    NB: this used to read geo_files[0] only, so the box was the footprint of
    whichever single granule the glob returned first, from any orbit and any date
    present on disk. With granules of two orbits side by side the GLC30 tiles were
    selected for the wrong part of the world, and even within the right orbit one
    granule is a small fraction of the swath.
    """
    geo_files = sorted((output_dir / "geo").glob("*.h5"))
    if orbit:
        geo_files = [f for f in geo_files if f"_{orbit}_" in f.name]
    if date:
        ymd = date.replace("-", "")
        geo_files = [f for f in geo_files if f"_{ymd}T" in f.name]

    try:
        import h5py  # noqa: F401  (used via _geo_granule_box)
    except ImportError:
        geo_files = []

    if geo_files:
        boxes = []
        for gf in geo_files:
            try:
                boxes.append(_geo_granule_box(gf))
            except (OSError, KeyError) as err:
                logger.warning(f"cannot read {gf.name} for the footprint: {err}")
        return _union_boxes_in_domain(boxes, "GEO granule(s) on disk")

    if not date:
        return None
    return _union_boxes_in_domain(_cmr_geo_boxes(orbit, date),
                                  "GEO granule(s) in CMR")


def _parse_zip_lon_range(filename: str) -> tuple[float, float] | None:
    """Parse longitude range from a Zenodo GLC30 ZIP filename.

    E.g. 'GLC_FCS30D_19852022maps_E10-E15.zip' -> (10, 20)
         'GLC_FCS30D_19852022maps_W95-W100.zip' -> (-100, -90)
    Each ZIP covers 10 degrees (two 5-degree tiles).
    """
    m = re.search(r'_([EW])(\d+)-[EW](\d+)\.zip$', filename)
    if not m:
        return None
    direction, start, end = m.group(1), int(m.group(2)), int(m.group(3))
    if direction == 'E':
        return (start, end + 5)
    else:
        return (-(end + 5), -start)


def download_glc30(output_dir: Path, orbit: str | None = None) -> None:
    logger.info("\n--- GLC30 / GLC_FCS30D (Zenodo) ---")
    dest_dir = output_dir / "glc30"
    dest_dir.mkdir(parents=True, exist_ok=True)

    # The bbox is resolved before any skip decision: which tiles are needed
    # depends on where this orbit flies, so "some .tif is present" cannot answer
    # that question (see the strip check further down).
    bbox = _get_orbit_bbox(output_dir, orbit=orbit, date=DATE)
    if not bbox:
        logger.info("  [skip] Cannot determine orbit bbox (no GEO file yet)")
        logger.info("         Run ECOSTRESS download first, then re-run.")
        return

    lat_min, lat_max, lon_min, lon_max = bbox
    logger.info(f"  Orbit bbox: lat [{lat_min:.1f}, {lat_max:.1f}], lon [{lon_min:.1f}, {lon_max:.1f}]")

    resp = requests.get(GLC30_API_URL)
    resp.raise_for_status()
    all_files = resp.json().get("files", [])
    zip_files = [f for f in all_files if f["key"].endswith(".zip")]

    needed = []
    for zf in zip_files:
        rng = _parse_zip_lon_range(zf["key"])
        if rng is None:
            continue
        zf_lon_min, zf_lon_max = rng
        if zf_lon_max > lon_min and zf_lon_min < lon_max:
            needed.append((zf, rng))

    if not needed:
        logger.warning("No matching longitude strips found")
        return

    # Which 5-degree longitude strips have tiles on disk. Every strip a ZIP covers
    # and the footprint needs must be present before that ZIP can be skipped: the
    # S3 path filters tiles by bbox, so a strip pair can sit on disk half-complete
    # and "any tile from this ZIP" would then skip a genuinely missing strip.
    on_disk_lons = set()
    for tif in dest_dir.glob("*.tif"):
        tile = _parse_glc30_tile_bounds(tif.name)
        if tile is not None:
            on_disk_lons.add(int(tile[2]))
    fetched = _glc30_fetched_zips(dest_dir)

    pending = []
    for zf, (zf_lon_min, zf_lon_max) in needed:
        if zf["key"] in fetched:
            logger.info(f"  [have] {zf['key']} — already downloaded")
            continue
        strips = {x for x in range(int(zf_lon_min), int(zf_lon_max), 5)
                  if x + 5 > lon_min and x < lon_max}
        if strips and strips <= on_disk_lons:
            logger.info(f"  [have] {zf['key']} — tiles for lon "
                  f"{sorted(strips)} already on disk")
            continue
        pending.append(zf)

    if not pending:
        logger.info(f"  [skip] the {len(needed)} longitude strip(s) covering this "
              "footprint are already present")
        return

    logger.info(f"  {len(pending)} ZIP(s) to download: {[z['key'] for z in pending]}")

    for zf_info in pending:
        zf_name = zf_info["key"]
        zf_url = zf_info["links"]["self"]
        zf_size = zf_info.get("size", 0)
        logger.info(f"  [download] {zf_name} ({zf_size / 1e9:.1f} GB) ...")

        with tempfile.NamedTemporaryFile(suffix=".zip", delete=False) as tmp:
            tmp_path = Path(tmp.name)
        try:
            with requests.get(zf_url, stream=True) as r:
                r.raise_for_status()
                downloaded = 0
                # Every GB, not every chunk: this archive is tens of GB and the
                # old in-place \r counter would now be one log record per MB.
                next_mark = 1e9
                with open(tmp_path, "wb") as f:
                    for chunk in r.iter_content(chunk_size=1024 * 1024):
                        f.write(chunk)
                        downloaded += len(chunk)
                        if downloaded >= next_mark:
                            logger.info(f"  [{downloaded / 1e9:.1f} / "
                                        f"{zf_size / 1e9:.1f} GB]")
                            next_mark += 1e9

            with zipfile.ZipFile(tmp_path) as zf:
                annual_tiles = [n for n in zf.namelist() if "Annual" in n and n.endswith(".tif")]
                for tile_path in annual_tiles:
                    tile_name = Path(tile_path).name
                    dest = dest_dir / tile_name
                    if dest.exists():
                        continue
                    dest.write_bytes(zf.read(tile_path))
                logger.info(f"  [ok] Extracted {len(annual_tiles)} Annual tiles from {zf_name}")
            _record_glc30_zip(dest_dir, zf_name)
        finally:
            tmp_path.unlink(missing_ok=True)


# ── ERA5 ───────────────────────────────────────────────────────────────────

def _era5_granule_hours(output_dir: Path, orbit: str | None = None) -> list[int]:
    """UTC acquisition hours of the ECOSTRESS granules on disk for DATE."""
    ymd = datetime.strptime(DATE, "%Y-%m-%d").strftime("%Y%m%d")
    hours = set()
    for sub in ("rad", "geo"):
        for f in (output_dir / sub).glob("*.h5"):
            if orbit and f"_{orbit}_" not in f.name:
                continue
            m = re.search(rf"_{ymd}T(\d{{2}})\d{{4}}_", f.name)
            if m:
                hours.add(int(m.group(1)))
    return sorted(hours)


def _era5_required_times(output_dir: Path,
                         orbit: str | None = None) -> tuple[list[datetime],
                                                            list[datetime]]:
    """ERA5 timestamps needed to process the granules on disk.

    TES and STIC both interpolate between the acquisition hour H and H+1, and
    STIC additionally reads single-level t2m at ERA5_TAMAX_HOURS for the daily
    max Ta. Returns (single_level_times, pressure_level_times); H+1 may roll
    over into the next day for a late-evening acquisition.
    """
    base = datetime.strptime(DATE, "%Y-%m-%d")
    hours = _era5_granule_hours(output_dir, orbit)
    if hours:
        logger.info(f"  acquisition hours on {DATE} — "
                    + ", ".join(f"{h:02d}h" for h in hours))
        pressure = {base + timedelta(hours=h + off) for h in hours for off in (0, 1)}
        single = pressure | {base + timedelta(hours=h) for h in ERA5_TAMAX_HOURS}
    else:
        logger.warning(f"no ECOSTRESS granule on disk for {DATE} — "
              "falling back to default hours")
        single = {base + timedelta(hours=h) for h in ERA5_SINGLE_HOURS}
        pressure = {base + timedelta(hours=h) for h in ERA5_PRESSURE_HOURS}
    return sorted(single), sorted(pressure)


def _era5_path(dest_dir: Path, prefix: str, t: datetime) -> Path:
    return dest_dir / f"{prefix}-{t.strftime('%Y_%m_%d_%H')}:00.grib"


def _split_era5_grib(src: Path, dest_dir: Path, prefix: str,
                     times: list[datetime]) -> bool:
    """Split a multi-hour ERA5 GRIB into one file per hour.

    Messages are bucketed on their own validity time, so any subset of hours
    splits correctly. Returns False (leaving src alone) when a requested hour
    is absent, so the caller can fall back to per-hour requests.
    """
    try:
        import eccodes
    except ImportError:
        logger.warning("eccodes not available — cannot split GRIB")
        return False

    buckets: dict[datetime, list[bytes]] = {}
    with open(src, "rb") as fh:
        while True:
            gid = eccodes.codes_grib_new_from_file(fh)
            if gid is None:
                break
            try:
                stamp = (datetime.strptime(
                            str(eccodes.codes_get(gid, "validityDate")), "%Y%m%d")
                         + timedelta(hours=int(eccodes.codes_get(
                             gid, "validityTime")) // 100))
                buckets.setdefault(stamp, []).append(
                    eccodes.codes_get_message(gid))
            finally:
                eccodes.codes_release(gid)

    missing = [t for t in times if t not in buckets]
    if missing:
        logger.warning("GRIB has no messages for "
                       + ", ".join(t.strftime("%H:00") for t in missing))
        return False

    for t in times:
        _era5_path(dest_dir, prefix, t).write_bytes(b"".join(buckets[t]))
    logger.info(f"  [ok] split into {len(times)} hourly files")
    src.unlink(missing_ok=True)
    return True


def download_era5(output_dir: Path, orbit: str | None = None) -> None:
    logger.info("\n--- ERA5 ---")
    try:
        import cdsapi
    except ImportError:
        logger.info("  [skip] cdsapi not installed — run: pip install cdsapi")
        logger.info("         Then configure ~/.cdsapirc with your CDS API key.")
        _print_era5_filenames(output_dir, orbit)
        return

    dest_dir = output_dir / "era5"
    dest_dir.mkdir(parents=True, exist_ok=True)
    single_times, pressure_times = _era5_required_times(output_dir, orbit)

    cds_url = os.environ.get("CDS_API_URL", "")
    cds_key = os.environ.get("CDS_API_KEY", "")
    if cds_url and cds_key:
        client = cdsapi.Client(url=cds_url, key=cds_key)
    else:
        client = cdsapi.Client()

    _download_era5_dataset(client, dest_dir, "single_levels", single_times)
    _download_era5_dataset(client, dest_dir, "pressure_levels", pressure_times)


# kind -> (CDS dataset, filename prefix, variables, pressure levels)
ERA5_DATASETS = {
    "single_levels": ("reanalysis-era5-single-levels", "era5_single_levels",
                      ERA5_SINGLE_VARS, None),
    "pressure_levels": ("reanalysis-era5-pressure-levels", "era5_37levels",
                        ERA5_PRESSURE_VARS, ERA5_PRESSURE_LEVELS),
}


def _download_era5_dataset(client, dest_dir: Path, kind: str,
                           times: list[datetime]) -> None:
    """Fetch the given hours of one ERA5 dataset, one bulk request per day.

    The bulk GRIB is split locally into the hourly files the algorithms expect;
    only if the split fails do we fall back to one request per hour (much
    slower — each CDS request queues separately).
    """
    dataset, prefix, variables, levels = ERA5_DATASETS[kind]
    missing = [t for t in times if not _era5_path(dest_dir, prefix, t).exists()]
    if not missing:
        logger.info(f"  [skip] {prefix} ({len(times)} hourly files exist)")
        return

    by_day: dict[date, list[datetime]] = {}
    for t in missing:
        by_day.setdefault(t.date(), []).append(t)

    for day, day_times in sorted(by_day.items()):
        req = {
            "product_type": "reanalysis",
            "variable": variables,
            "year": str(day.year),
            "month": f"{day.month:02d}",
            "day": f"{day.day:02d}",
            "area": ERA5_AREA,
            "data_format": "grib",
        }
        if levels:
            req["pressure_level"] = levels
        hours = ", ".join(f"{t.hour:02d}h" for t in day_times)
        logger.info(f"  [download] {prefix} {day} ({hours}) ...")
        bulk = dest_dir / f"{prefix}-{day:%Y_%m_%d}-bulk.grib"
        client.retrieve(
            dataset,
            {**req, "time": [f"{t.hour:02d}:00" for t in day_times]},
            str(bulk),
        )
        if _split_era5_grib(bulk, dest_dir, prefix, day_times):
            continue
        for t in day_times:
            dest = _era5_path(dest_dir, prefix, t)
            if dest.exists():
                continue
            logger.info(f"  [download] {prefix} {day} {t.hour:02d}h individually ...")
            client.retrieve(dataset, {**req, "time": f"{t.hour:02d}:00"},
                            str(dest))
        bulk.unlink(missing_ok=True)
    logger.info(f"  [ok] {prefix}")


def _print_era5_filenames(output_dir: Path, orbit: str | None = None) -> None:
    single_times, pressure_times = _era5_required_times(output_dir, orbit)
    logger.info(f"\n  ERA5 files needed in {output_dir / 'era5'}:")
    for prefix, times in (("era5_single_levels", single_times),
                          ("era5_37levels", pressure_times)):
        for t in times:
            logger.info(f"    {_era5_path(Path(), prefix, t).name}")


# ── S3 sync (alternative to individual downloads) ────────────────────────

_S3_PREFIX_DEFAULTS = {
    "rad": "L1B_RAD_V002",
    "geo": "L1B_GEO_V002",
    "cloud": "L2_CLOUD_V002",
    "era5": "ERA5",
    "mota": "MOTA",
    "fcover": "FCOVER",
    "lai": "LAI",
    "lulc": "LULC",
    "oco2": "OCO2",
    "glc30": "GLC_FCS30D",
    "parh": "PARH",
}

# Override any prefix via S3_PREFIX_RAD=my/custom/path, etc.
_ENV_TO_KEY = {f"S3_PREFIX_{k.upper()}": k for k in _S3_PREFIX_DEFAULTS}


def _get_s3_prefix_map() -> dict[str, str]:
    """Build {s3_prefix: local_subdir} map, applying env overrides."""
    result = {}
    for local_key, default_prefix in _S3_PREFIX_DEFAULTS.items():
        env_var = f"S3_PREFIX_{local_key.upper()}"
        prefix = os.environ.get(env_var, default_prefix)
        result[prefix] = local_key
    return result


def _get_s3_client():
    endpoint = os.environ.get("S3_ENDPOINT_URL", "")
    access_key = os.environ.get("S3_ACCESS_KEY", "")
    secret_key = os.environ.get("S3_SECRET_KEY", "")
    if not endpoint or not access_key or not secret_key:
        return None
    try:
        import boto3
    except ImportError:
        logger.info("  [skip] boto3 not installed — run: pip install boto3")
        return None
    return boto3.client(
        "s3",
        endpoint_url=endpoint,
        aws_access_key_id=access_key,
        aws_secret_access_key=secret_key,
    )


def _s3_date_matches(fname: str, date_str: str, monthly: bool = False,
                     doy: bool = False) -> bool:
    """Check if a filename matches the target date.

    With monthly=True, matches by YYYYMM or previous month
    (for 10-day composites like FCOVER/LAI that may span month boundaries).
    With doy=True, also matches YYYY+DOY format (e.g. A2025001 for Jan 1).
    Otherwise matches exact YYYYMMDD or YYYY_MM_DD.
    """
    ymd = date_str.replace("-", "")
    if monthly:
        dt = datetime.strptime(ymd[:6] + "01", "%Y%m%d")
        prev = (dt - timedelta(days=1)).strftime("%Y%m")
        return ymd[:6] in fname or prev in fname
    ymd_us = date_str.replace("-", "_")
    if ymd in fname or ymd_us in fname:
        return True
    if doy:
        dt = datetime.strptime(ymd[:8], "%Y%m%d")
        doy_key = f"{dt.year}{dt.timetuple().tm_yday:03d}"
        return doy_key in fname
    return False


# Datasets that are global (not date-specific) — sync all files
_S3_GLOBAL_LOCAL_KEYS = {"lulc", "glc30"}
_S3_MONTHLY_LOCAL_KEYS = {"fcover", "lai"}


def _extract_orbit_number(fname: str) -> str | None:
    """Extract 5-digit orbit number from ECOSTRESS filename (e.g. '36798' from '..._36798_018_...')."""
    m = re.search(r'_(\d{5})_\d{3}_', fname)
    return m.group(1) if m else None


def _parse_glc30_tile_bounds(fname: str) -> tuple[float, float, float, float] | None:
    """Parse GLC30 tile bounds from filename like 'GLC_FCS30D_20162022_E0N10_Annual.tif'.

    Returns (lat_min, lat_max, lon_min, lon_max) or None.
    Each tile covers 5° lon × 5° lat.
    """
    m = re.search(r'_([EW])(\d+)([NS])(\d+)_', fname)
    if not m:
        return None
    ew, lon_val, ns, lat_val = m.group(1), int(m.group(2)), m.group(3), int(m.group(4))
    lon = lon_val if ew == 'E' else -lon_val
    lat = lat_val if ns == 'N' else -lat_val
    return (lat - 5, lat, lon, lon + 5)


def _glc30_tile_overlaps_bbox(fname: str, bbox: tuple) -> bool:
    """Check if a GLC30 tile overlaps with orbit bounding box."""
    tile = _parse_glc30_tile_bounds(fname)
    if tile is None:
        return True  # can't parse → don't filter out
    t_lat_min, t_lat_max, t_lon_min, t_lon_max = tile
    lat_min, lat_max, lon_min, lon_max = bbox
    return (t_lat_max > lat_min and t_lat_min < lat_max and
            t_lon_max > lon_min and t_lon_min < lon_max)


def _s3_orbit_matches(fname: str, orbit: str) -> bool:
    """Check if a filename contains the orbit number (e.g. '36798')."""
    return orbit in fname


def _s3_domain_scene_keys(skip: set[str], orbit: str | None,
                          whole_orbit: bool) -> set[str] | None:
    """Scene keys inside EEH2_DOMAIN, or None when no filter applies.

    The bucket carries no footprint metadata, so the domain has to be resolved
    the same way the CMR path resolves it: ask CMR which GEO scenes intersect
    the box, then match S3 filenames on the shared {orbit}_{scene} key. Without
    this, a narrowed EEH2_DOMAIN was silently ignored whenever S3 was reachable
    and the sync pulled the whole orbit — hundreds of GB for a few useful
    granules. Degrades to None (no filter) rather than aborting: an S3-only
    deployment may have no Earthdata credentials at all.
    """
    if whole_orbit or {"rad", "geo", "cloud"} <= skip:
        return None
    bbox = processing_domain()
    if not bbox:
        return None
    try:
        keys = _domain_scene_keys(_earthdata_session(), orbit or "", bbox)
    except Exception as e:
        logger.warning(f"domain filter unavailable ({e}) — S3 sync will not "
                       "apply EEH2_DOMAIN")
        return None
    n, w, s, e_ = bbox
    logger.info(f"\n  Processing domain: lat [{s}, {n}], lon [{w}, {e_}] — "
                f"{len(keys)} scene(s) inside it")
    if not keys:
        logger.warning("no granule intersects the processing domain on this date")
    return keys


def download_from_s3(output_dir: Path, date_str: str,
                     skip: set[str] | None = None,
                     orbit: str | None = None,
                     whole_orbit: bool = False) -> set[str]:
    """Sync data from S3 bucket for a specific date.
    Returns set of local subdirs that were populated."""
    bucket = os.environ.get("S3_BUCKET", "ECOSTRESS")
    s3 = _get_s3_client()
    if s3 is None:
        return set()

    orbit_msg = f", orbit: {orbit}" if orbit else ""
    logger.info(f"\n{'='*60}")
    logger.info(f"  S3 sync from s3://{bucket}/ (date: {date_str}{orbit_msg})")
    logger.info(f"{'='*60}")

    skip = skip or set()
    synced = set()
    rad_orbits: set[str] = set()
    orbit_bbox = _get_orbit_bbox(output_dir, orbit=orbit, date=date_str)
    domain_keys = _s3_domain_scene_keys(skip, orbit, whole_orbit)

    for s3_prefix, local_subdir in _get_s3_prefix_map().items():
        if local_subdir in skip:
            logger.info(f"\n  [skip] {s3_prefix} (--skip flag)")
            continue

        dest_dir = output_dir / local_subdir
        dest_dir.mkdir(parents=True, exist_ok=True)

        if local_subdir == "glc30" and orbit_bbox is None:
            # Without a footprint the tile filter cannot run, and the global
            # GLC30 set is hundreds of GB. Skip rather than fetch all of it.
            logger.info(f"  [skip] {s3_prefix} (footprint unknown)")
            continue

        is_global = local_subdir in _S3_GLOBAL_LOCAL_KEYS
        is_monthly = local_subdir in _S3_MONTHLY_LOCAL_KEYS
        is_doy = local_subdir == "mota"
        filter_bbox = local_subdir == "glc30" and orbit_bbox is not None
        needs_rad_match = local_subdir in ("geo", "cloud") and rad_orbits
        extra = ""
        if filter_bbox:
            lat_min, lat_max, lon_min, lon_max = orbit_bbox
            extra = f" (bbox-filtered: lat [{lat_min:.1f}, {lat_max:.1f}], lon [{lon_min:.1f}, {lon_max:.1f}])"
        elif local_subdir == "lulc":
            extra = " (the one raster STIC reads)"
        elif is_global:
            extra = " (global)"
        elif is_monthly:
            extra = f" (monthly: {date_str[:7]})"
        elif needs_rad_match:
            extra = f" (filtered to {len(rad_orbits)} RAD orbit(s))"
        logger.info(f"\n  --- {s3_prefix} → {local_subdir}/{extra} ---")

        try:
            paginator = s3.get_paginator("list_objects_v2")
            pages = paginator.paginate(Bucket=bucket, Prefix=s3_prefix + "/")

            # Collect candidates (for monthly datasets, we pick best RT after)
            candidates = []
            outside = 0
            for page in pages:
                for obj in page.get("Contents", []):
                    key = obj["Key"]
                    fname = Path(key).name
                    if not fname or fname.startswith("."):
                        continue
                    if fname.endswith((".xml", ".dmrpp")):
                        continue
                    if not is_global and not _s3_date_matches(fname, date_str, monthly=is_monthly, doy=is_doy):
                        continue
                    if filter_bbox and not _glc30_tile_overlaps_bbox(fname, orbit_bbox):
                        continue
                    if local_subdir == "lulc" and fname != LULC_FILENAME:
                        # The bucket carries all four LC100 rasters: the 2018-conso
                        # and 2019-nrt epochs, each as a discrete-classification map
                        # and as its per-class probability stack. STIC opens exactly
                        # one of them -- Read_LULC hardcodes the 2018-conso map -- and
                        # the Zenodo fallback only ever fetches that one. Syncing the
                        # prefix wholesale pulled 17 GB that nothing ever opens, and
                        # deleting the surplus locally did not help: the next run
                        # fetched it again.
                        continue
                    if orbit and local_subdir in ("rad", "geo", "cloud") and not _s3_orbit_matches(fname, orbit):
                        continue
                    if domain_keys is not None and local_subdir in ("rad", "geo", "cloud"):
                        if _scene_key(fname) not in domain_keys:
                            outside += 1
                            continue
                    if needs_rad_match:
                        fnum = _extract_orbit_number(fname)
                        if fnum and fnum not in rad_orbits:
                            continue
                    candidates.append((fname, key, obj))

            if outside:
                logger.info(f"  {outside} granule(s) outside the processing "
                            "domain — not downloaded")

            # For monthly datasets (FCOVER/LAI 10-day composites): take the
            # best RT of every composite the date needs, which is more than one.
            # STIC wants the composite containing the date; GPP interpolates
            # between the two that bracket it (see _clms_required_dates).
            if is_monthly and candidates:
                by_date: dict[str, list[tuple]] = {}
                for item in candidates:
                    cdate = _clms_composite_date(item[0])
                    if cdate:
                        by_date.setdefault(cdate, []).append(item)
                selected = []
                for want in _clms_required_dates(date_str):
                    if want not in by_date:
                        # No substitution of a neighbouring composite: the
                        # consumers select by name, so a stand-in only moves the
                        # failure further from its cause. The CDSE step that
                        # follows can still fill the gap.
                        logger.warning(f"composite {want} not on S3")
                        continue
                    best_item = max(by_date[want],
                                    key=lambda item: _clms_rt_score(item[0]))
                    score = _clms_rt_score(best_item[0])
                    rt_label = f"RT{score}" if score < 100 else "final"
                    logger.info(f"  Best for {want}: {best_item[0]} ({rt_label})")
                    selected.append(best_item)
                candidates = selected

            count = 0
            for fname, key, obj in candidates:
                dest = dest_dir / fname
                if dest.exists() and dest.stat().st_size == obj["Size"]:
                    count += 1
                    logger.info(f"  [{count}] {fname} (exists, skip)")
                    if local_subdir == "rad":
                        onum = _extract_orbit_number(fname)
                        if onum:
                            rad_orbits.add(onum)
                    continue
                s3.download_file(bucket, key, str(dest))
                count += 1
                logger.info(f"  [{count}] {fname} ({obj['Size'] / 1e6:.1f} MB)")
                if local_subdir == "rad":
                    onum = _extract_orbit_number(fname)
                    if onum:
                        rad_orbits.add(onum)
            if count == 0:
                logger.info(f"  [skip] no files found on S3")
            else:
                logger.info(f"  [ok] {count} file(s) synced")
                synced.add(local_subdir)
        except Exception as e:
            logger.warning(f"S3 sync failed for {s3_prefix}: {e}")

    return synced


# ── Manual data instructions ──────────────────────────────────────────────

def print_manual_data(par_ok: bool) -> None:
    if par_ok:
        return
    logger.info("\n" + "="*60)
    logger.info("  MANUAL STEP REQUIRED")
    logger.info("="*60)
    logger.info("""
PAR half-hourly (PARin*.nc files) could not be downloaded automatically.
Order manually from CM SAF SAFIRA: https://wui.cmsaf.eu
  Select SARAH-3, product "PAR instantaneous". Once processed,
  download via the FTP/HTTPS link in the confirmation email:
    wget -r -np -nd -A "*.nc" <order_url>
  Place PARin*.nc files in: data/parh/

Note: TES and STIC run without PAR — it is only needed for GPP.
""")


# ── Main ──────────────────────────────────────────────────────────────────

def main():
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s — %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    parser = argparse.ArgumentParser(
        description="Download sample data for orbit 36798_018_20250101T165942.",
    )
    parser.add_argument("--output-dir", default="./data",
                        help="Root directory for downloaded data (default: ./data).")
    parser.add_argument("--skip-era5", action="store_true",
                        help="Skip ERA5 download (requires cdsapi + CDS account).")
    parser.add_argument("--skip-ecostress", action="store_true",
                        help="Skip ECOSTRESS + OCO-2 download (MOTA is always downloaded).")
    parser.add_argument("--skip-clms", action="store_true",
                        help="Skip FCOVER/LAI download (requires CDSE account).")
    parser.add_argument("--skip-lulc", action="store_true",
                        help="Skip LULC download from Zenodo.")
    parser.add_argument("--skip-glc30", action="store_true",
                        help="Skip GLC30 download from Zenodo (large ZIPs).")
    parser.add_argument("--skip-par", action="store_true",
                        help="Skip PAR download (requires EUMETSAT account + eumdac).")
    parser.add_argument("--date", default=None,
                        help="Target date YYYY-MM-DD (default: 2025-01-01). "
                             "Affects S3 filtering and ancillary data (ERA5, PAR, FCOVER, LAI, MOTA).")
    parser.add_argument("--orbit", default=None,
                        help="Filter S3 ECOSTRESS files by orbit number (e.g. 36798). "
                             "Only affects rad/geo/cloud from S3. Without this flag, "
                             "all orbits for the target date are synced.")
    parser.add_argument("--whole-orbit", action="store_true",
                        help="Download every granule of the orbit, including those "
                             "outside the processing domain (see EEH2_DOMAIN). Those "
                             "granules have no ERA5 coverage and TES will skip them.")
    parser.add_argument("--no-s3", action="store_true",
                        help="Skip S3 sync entirely, even when S3_ENDPOINT_URL is set. "
                             "Use when you know the data is not on S3.")
    parser.add_argument("--env", default=".env", metavar="FILE",
                        help="Path to .env file for credentials (default: .env).")
    args = parser.parse_args()

    # Apply --date override to module-level constants
    global DATE, TEMPORAL
    if args.date is not None:
        DATE = args.date
    TEMPORAL = f"{DATE}T00:00:00Z,{DATE}T23:59:59Z"

    env_path = Path(args.env)
    if env_path.exists():
        try:
            from dotenv import dotenv_values
            for k, v in dotenv_values(env_path).items():
                if v and k not in os.environ:
                    os.environ[k] = v
        except ImportError:
            pass

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # ── S3 sync (if configured) ──────────────────────────────────────
    s3_synced = set()
    if not args.no_s3 and os.environ.get("S3_ENDPOINT_URL"):
        skip_s3 = set()
        if args.skip_ecostress:
            skip_s3.update(["rad", "geo", "cloud", "oco2"])
        if args.skip_era5:
            skip_s3.add("era5")
        if args.skip_clms:
            skip_s3.update(["fcover", "lai"])
        if args.skip_lulc:
            skip_s3.add("lulc")
        if args.skip_glc30:
            skip_s3.add("glc30")
        if args.skip_par:
            skip_s3.add("parh")
        s3_synced = download_from_s3(output_dir, date_str=DATE, skip=skip_s3,
                                     orbit=args.orbit,
                                     whole_orbit=args.whole_orbit)

    # ── Normal downloads (skip what S3 already provided) ─────────────
    errors = []

    session = None
    dl_failed = 0
    if not args.skip_ecostress and "rad" not in s3_synced:
        try:
            session = _earthdata_session()
            total, failed = download_ecostress(session, output_dir,
                                               orbit=args.orbit,
                                               whole_orbit=args.whole_orbit)
            dl_failed += failed
            if failed:
                logger.warning(f"{failed}/{total} ECOSTRESS files failed")
            download_oco2(session, output_dir)
        except Exception as e:
            logger.error(f"ECOSTRESS/OCO2 download failed: {e}")
            errors.append("ecostress")

    # Ancillaries only serve RAD granules; an orbit outside the domain used to
    # pull MODIS, ERA5, CLMS and PAR for a run that then stopped on "no RAD".
    if not args.skip_ecostress:
        ymd = DATE.replace("-", "")
        pattern = f"*_{args.orbit}_*_{ymd}T*.h5" if args.orbit else f"*_{ymd}T*.h5"
        if not list((output_dir / "rad").glob(pattern)):
            logger.error(f"No RAD granule for {DATE}"
                         + (f", orbit {args.orbit}" if args.orbit else "")
                         + " inside the processing domain — ancillary data not downloaded")
            sys.exit(1)

    if "mota" not in s3_synced:
        try:
            if session is None:
                session = _earthdata_session()
            download_mota(session, output_dir)
        except Exception as e:
            logger.error(f"MOTA download failed: {e}")
            errors.append("mota")

    if not args.skip_era5 and "era5" not in s3_synced:
        try:
            download_era5(output_dir, orbit=args.orbit)
        except Exception as e:
            logger.error(f"ERA5 download failed: {e}")
            errors.append("era5")

    # Gated on what is actually missing, not on whether the S3 sync ran: the
    # bucket can hold one of GPP's two brackets and not the other, and the old
    # 'fcover in s3_synced' test then locked CDSE out of filling the gap.
    if not args.skip_clms and _clms_missing(output_dir, DATE):
        try:
            download_fcover_lai(output_dir)
        except Exception as e:
            logger.error(f"FCOVER/LAI download failed: {e}")
            errors.append("clms")

    par_ok = True
    if not args.skip_par and "parh" not in s3_synced:
        try:
            par_ok = download_par(output_dir)
        except Exception as e:
            logger.error(f"PAR download failed: {e}")
            errors.append("par")
            par_ok = False
    elif "parh" in s3_synced:
        # A synced parh/ never reaches download_par(), so validate here too.
        # The bucket itself stores PARin unpacked, but a file that arrived by
        # some other route can still be sitting there, and it would fail every
        # granule with "NetCDF: Unknown file format" rather than fail here.
        logger.info("\n--- PAR half-hourly (from S3) ---")
        repair_par_archives(output_dir)
        ymd = DATE.replace("-", "")
        par_ok = any(_is_valid_par_netcdf(f) for f in
                     (output_dir / "parh").glob(f"{PAR_PRODUCT_PREFIX}{ymd}*.nc"))
        if not par_ok:
            logger.warning(f"no readable {PAR_PRODUCT_PREFIX} NetCDF for {DATE}")

    if not args.skip_lulc and "lulc" not in s3_synced:
        try:
            download_lulc(output_dir)
        except Exception as e:
            logger.error(f"LULC download failed: {e}")
            errors.append("lulc")

    if not args.skip_glc30 and "glc30" not in s3_synced:
        try:
            download_glc30(output_dir, orbit=args.orbit)
        except Exception as e:
            logger.error(f"GLC30 download failed: {e}")
            errors.append("glc30")

    if errors:
        logger.warning(f"Some downloads failed: {', '.join(errors)}")
    if dl_failed:
        logger.warning(f"{dl_failed} file(s) failed after retries")

    print_manual_data(par_ok)

    logger.info("="*60)
    logger.info(f"  Done. Data root: {output_dir.resolve()}")
    logger.info("="*60)

    if dl_failed or errors:
        sys.exit(2)


if __name__ == "__main__":
    main()
