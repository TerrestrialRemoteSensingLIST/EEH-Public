# EEH2 — ECOSTRESS European Hub Processing Pipeline

Algorithms developed for the ESA EEH project to generate land surface temperature (TES), evapotranspiration (STIC), and gross primary productivity / water use efficiency (GPP) from ECOSTRESS satellite data.

More information: https://ecostresshub.eu/

## Pipeline overview

```
ECOSTRESS L1B_RAD ──► TES ──► STIC ──► GPP ──► COG (optional)
                       │        │        │        │
                       ▼        ▼        ▼        ▼
                      LST    ET/LE/H   GPP/WUE   GeoTIFF (4 products)
```

| Step | Input | Output | Key dependencies |
|------|-------|--------|------------------|
| **TES** | L1B_RAD, L1B_GEO, ERA5 | LST, emissivity | RTTOV 12.3 |
| **STIC** | L2_LSTE (TES output), L1B_GEO, CLOUD, FCOVER, MCD43C3, LULC, ERA5 | LE, H, G, Rn, ETD | GDAL |
| **GPP** | STIC output, L1B_GEO, CLOUD, PAR, LAI, FCOVER, ERA5, OCO-2, GLC30, LUT | GPPd, WUEd | — |
| **COG** | TES/STIC/GPP H5 outputs, L1B_GEO, L2_CLOUD | Cloud Optimized GeoTIFF (4 products) | pyresample, GDAL |

### Data dependency matrix

Each pipeline step requires specific ancillary datasets. **TES and STIC can run without PAR, OCO-2, and GLC30** — those are only needed for GPP.

| Dataset | TES | STIC | GPP | Auto-download | Notes |
|---------|:---:|:----:|:---:|:-------------:|-------|
| ECOSTRESS L1B_RAD | **x** | | | yes | Earthdata |
| ECOSTRESS L1B_GEO | **x** | **x** | **x** | yes | Earthdata |
| ECOSTRESS L2_CLOUD | | **x** | **x** | yes | Earthdata |
| ERA5 single-levels | **x** | **x** | **x** | yes | CDS API |
| ERA5 pressure-levels | **x** | | | yes | CDS API |
| MCD43C3 MODIS albedo | | **x** | | yes | Earthdata |
| FCOVER 300m 10-day | | **x** | **x** | yes | CDSE |
| LAI 300m 10-day | | | **x** | yes | CDSE |
| LULC (PROBAV LC100) | | **x** | | yes | Zenodo (free, ~1.6 GB) |
| OCO-2 XCO2 | | | **x** | yes | Earthdata (coverage 2018–2022, falls back to closest year) |
| GLC30 land-cover tiles | | | **x** | yes | Zenodo (free, ~3–6 GB per longitude strip) |
| PAR half-hourly (PARin) | | | **x** | partial | EUMETSAT Data Store or manual order via [CM SAF SAFIRA](https://wui.cmsaf.eu) |
| LUT CSV | | | **x** | — | Shipped with this repository |

**Out-of-coverage behavior (GPP):**
- **OCO-2**: if the scene date is outside 2018–2022, the code falls back to the closest available year with the same day-of-year. Scene is skipped only if no OCO-2 files are present at all.
- **PAR**: scene is skipped if no matching `PARin` file is found for that date.
- **GLC30**: scene is skipped if no tiles overlap the orbit footprint ("No vegetation tile").

## Quick start (Docker — recommended)

Docker is the simplest way to run the pipeline. It includes Python, RTTOV, GDAL, and all dependencies.

### Prerequisites

- [Docker](https://docs.docker.com/get-docker/) installed
- `rttov123.tar.gz` — RTTOV 12.3 source tarball (**provided with this delivery**)
- [NASA Earthdata](https://urs.earthdata.nasa.gov/) account (for ECOSTRESS data)
- [Copernicus CDS](https://cds.climate.copernicus.eu/) account (for ERA5 data)

### 1. Clone and configure

```bash
git clone https://git.list.lu/terrestrial-ecosystem-sensing/eeh-public.git
cd eeh-public
cp /path/to/rttov123.tar.gz .
```

> **RTTOV 12.3** is licensed by [NWP-SAF](https://nwp-saf.eumetsat.int/site/software/rttov/) and cannot be redistributed publicly. Version 12.3 is no longer available on the NWP-SAF website; the tarball is provided alongside this delivery. The Docker build compiles it automatically.

### 2. Configure credentials

```bash
cp .env.example .env
```

Edit `.env` and fill in the credentials section:

```dotenv
EARTHDATA_USERNAME=your_earthdata_user
EARTHDATA_PASSWORD=your_earthdata_pass
CDS_API_URL=https://cds.climate.copernicus.eu/api
CDS_API_KEY=your_uid:your_api_key
CDSE_USERNAME=your_cdse_user
CDSE_PASSWORD=your_cdse_pass
EUMETSAT_CONSUMER_KEY=your_consumer_key
EUMETSAT_CONSUMER_SECRET=your_consumer_secret
```

- **Earthdata** — free account at https://urs.earthdata.nasa.gov/ (ECOSTRESS, MOTA, OCO-2)
- **CDS API** — key from https://cds.climate.copernicus.eu/how-to-api (ERA5). The account must also accept the dataset licence once, logged in, at the bottom of the [ERA5 single levels download page](https://cds.climate.copernicus.eu/datasets/reanalysis-era5-single-levels?tab=download#manage-licences); the same licence covers the pressure-levels dataset. Otherwise every request fails with `403 required licences not accepted`
- **CDSE** — free account at https://dataspace.copernicus.eu/ (FCOVER, LAI)
- **EUMETSAT** — API keys from https://api.eumetsat.int/api-key/ (PAR, optional)
- **S3** — endpoint + credentials for an S3-compatible bucket containing pre-staged data (optional)

### 3. Download sample data

```bash
# Via Docker (recommended — no local Python setup needed)
docker compose run download

# Or locally
pip install requests cdsapi
python tools/download_sample.py --output-dir ./data
```

By default this downloads sample data for orbit **36798_018 on 2025-01-01**. Use `--date` to change the target date:

```bash
docker compose run download --date 2025-03-15
```

**Download flags:**

| Flag | Effect |
|------|--------|
| `--date YYYY-MM-DD` | Target date (default: `2025-01-01`) — affects S3 filtering and ancillary data |
| `--orbit NUMBER` | Filter S3 ECOSTRESS files (rad/geo/cloud) by orbit number (e.g. `36798`) |
| `--no-s3` | Skip S3 sync entirely, even when `S3_ENDPOINT_URL` is set |
| `--skip-ecostress` | Skip ECOSTRESS + OCO-2 (MOTA is always downloaded) |
| `--skip-era5` | Skip ERA5 |
| `--skip-clms` | Skip FCOVER + LAI |
| `--skip-par` | Skip PAR |
| `--skip-lulc` | Skip LULC |
| `--skip-glc30` | Skip GLC30 |
| `--output-dir DIR` | Root output directory (default: `./data`) |

> **Retry on transient errors:** HTTP 5xx and 429 errors are retried up to 4 times with exponential backoff (5s → 10s → 20s → 40s). Files that still fail are reported but don't block other downloads. The script exits with code **2** when some files failed, so the pipeline can continue with available data.

#### Alternative: download from S3

If data is pre-staged on an S3-compatible bucket (MinIO, AWS, etc.), add the S3 credentials to `.env`:

```dotenv
S3_ENDPOINT_URL=https://your-s3-endpoint.example.com
S3_ACCESS_KEY=your_access_key
S3_SECRET_KEY=your_secret_key
S3_BUCKET=ECOSTRESS
```

When `S3_ENDPOINT_URL` is set, the script automatically syncs from S3 first, then falls back to normal downloads for any missing datasets. Use `--no-s3` to skip S3 sync when you know the data isn't there.

Expected S3 bucket layout:

```
s3://ECOSTRESS/
├── L1B_RAD_V002/    # ECOSTRESS radiance
├── L1B_GEO_V002/    # ECOSTRESS geolocation
├── L2_CLOUD_V002/   # ECOSTRESS cloud mask
├── ERA5/            # ERA5 hourly GRIB files
├── MOTA/            # MCD43C3 MODIS albedo
├── FCOVER/          # FCOVER 300m .nc
├── LAI/             # LAI 300m .nc
├── LULC/            # PROBAV LC100 .tif
├── OCO2/            # OCO-2 XCO2 .nc4
├── GLC_FCS30D/      # GLC30 Annual tiles .tif
└── PARH/            # PARin*.nc
```

**ERA5 filenames contain a colon** (`era5_single_levels-2026_04_20_08:00.grib`).
The name is a contract shared by the downloader, the S3 bucket and the
`ReadERA5Data.py` of both STIC and TES, which rebuild it by concatenation. Fine on
Linux, macOS and in the container; a native Windows Python cannot open these files,
since NTFS has no colon and the bind mount maps it to U+F03A.

**S3 smart filtering:**
- **ECOSTRESS** (RAD/GEO/CLOUD): filtered by date; with `--orbit`, only matching orbit files are synced. GEO/CLOUD are further filtered to orbits that have RAD files.
- **FCOVER/LAI**: every 10-day composite the target date needs is selected (see the guards section — one for STIC, two for GPP), best available version each (final > RT6 > ... > RT0). What the bucket does not have, the CDSE step downloads.
- **GLC30**: filtered by orbit bounding box (only tiles covering the orbit footprint). The box comes from the GEO granules on disk, or from CMR metadata when they have not been downloaded yet; if it cannot be resolved, GLC30 is skipped rather than synced in full.
- **LULC**: only the one raster STIC opens. The bucket holds all four LC100 products -- the 2018-conso and 2019-nrt epochs, each as a discrete-classification map and as its per-class probability stack -- but `Read_LULC` hardcodes the 2018-conso map and the Zenodo fallback only fetches that one. The prefix used to be synced wholesale, which cost 17 GB that nothing ever opens; deleting the surplus locally did not help, the next run fetched it again.

If `S3_ENDPOINT_URL` is set in `.env`, S3 sync runs automatically (use `--no-s3` to skip).

If your bucket uses different directory names, override any prefix in `.env`:

```dotenv
S3_PREFIX_RAD=my_custom_radiance_dir
S3_PREFIX_GLC30=land_cover_tiles
```

Available overrides: `S3_PREFIX_RAD`, `S3_PREFIX_GEO`, `S3_PREFIX_CLOUD`, `S3_PREFIX_ERA5`, `S3_PREFIX_MOTA`, `S3_PREFIX_FCOVER`, `S3_PREFIX_LAI`, `S3_PREFIX_LULC`, `S3_PREFIX_OCO2`, `S3_PREFIX_GLC30`, `S3_PREFIX_PARH`.

### 4. Build and run

Both services (`eeh-pipeline` and `download`) run the same code from the same
Dockerfile and share the image tag `eeh2:latest`; only their entrypoint differs.
Always rebuild with plain `docker compose build` (or `--build`): building a single
service used to leave the other one running stale code, which is easy to miss
because it fails silently -- the old code simply prints its old messages.

```bash
# Rebuild both services after any code change
docker compose build

# Full pipeline (TES → STIC → GPP) on ALL files in data/
docker compose up --build

# Single orbit only (filters RAD/TES/STIC files by orbit number)
docker compose run eeh-pipeline --orbit 36798

# Single orbit + single step
docker compose run eeh-pipeline --orbit 36798 --step tes

# Full pipeline + COG conversion (produces .tif files in output/cog/)
docker compose run eeh-pipeline --orbit 36798 --cog

# COG conversion only (on existing H5 outputs)
docker compose run eeh-pipeline --orbit 36798 --step cog

# List available orbits for a date (no download, no processing)
docker compose run eeh-pipeline --start-date 2026-08-01 --end-date 2026-08-01 --list-orbits

# Date range: discover and process all orbits between two dates
docker compose run eeh-pipeline --start-date 2026-08-01 --end-date 2026-08-03
docker compose run eeh-pipeline --start-date 2026-08-01 --end-date 2026-08-03 --cog

# Date range + specific orbit (filters ECOSTRESS to that orbit only)
docker compose run eeh-pipeline --start-date 2026-08-01 --end-date 2026-08-01 --orbit 45031

# Skip S3 sync (when you know the data isn't on S3)
docker compose run eeh-pipeline --start-date 2026-08-01 --end-date 2026-08-01 --no-s3

# If the orbit isn't in data/rad/, --orbit auto-downloads from S3 or Earthdata,
# then downloads ancillary data for the orbit's date
docker compose run eeh-pipeline --orbit 45031

# Specific files
docker compose run eeh-pipeline --input-files /data/rad/ECOv002_L1B_RAD_36798_018_20250101T165942_0712_01.h5

# Single step (all files)
docker compose run eeh-pipeline --step tes
docker compose run eeh-pipeline --step stic
docker compose run eeh-pipeline --step gpp
```

**Pipeline flags:**

| Flag | Effect |
|------|--------|
| `--orbit NUMBER` | Process only files matching this orbit (auto-downloads if not found locally) |
| `--start-date YYYY-MM-DD` | Start of date range (inclusive). Downloads and processes all orbits in the range. Requires `--end-date` |
| `--end-date YYYY-MM-DD` | End of date range (inclusive). Requires `--start-date` |
| `--cog` | Convert H5 outputs to Cloud Optimized GeoTIFF after processing |
| `--step {tes,stic,gpp,cog}` | Run only one step |
| `--force` | Re-run steps even if output files already exist (default: skip) |
| `--download` | Download sample data before running |
| `--input-files FILE...` | Process specific RAD files |
| `--no-s3` | Skip S3 sync in downloads (passed to `tools/download_sample.py`) |
| `--list-orbits` | List available ECOSTRESS orbits for the date range, then exit. Requires `--start-date`/`--end-date` |
| `--env FILE` | Path to .env file (default: `.env`) |

> `--start-date`/`--end-date` and `--input-files` are mutually exclusive. `--orbit` can be combined with `--start-date`/`--end-date` to filter by orbit number within the date range.

> A named selection is followed through all three steps by `{orbit}_{scene}`, so
> `--input-files` on one granule runs STIC and GPP on that granule only.
> Narrowing on the date alone used to pull in every other granule of the same
> date, because that is what the products of a step have in common. A path that
> does not exist is refused before anything starts, rather than counted as one
> failed granule on the way to the completion banner.

**Exit codes**, for callers that script the pipeline:

| Code | Meaning |
|------|---------|
| `0` | Every enabled step produced output |
| `1` | A step ran, exited cleanly, and wrote nothing — every granule was skipped or failed |
| `2` | Bad invocation (unknown flag, missing input file) |
| other | Propagated verbatim from the step that failed |

Code `1` exists because TES, STIC and GPP each catch their per-granule failures
so one bad granule does not abandon the rest — which means they exit 0 even when
every granule failed. The closing banner has always said so; the exit code now
says so too.

## Verifying the installation

Two checks, one cheap and one conclusive.

### Environment check

```bash
docker compose run selftest
```

Seconds, no data and no credentials needed. It reports the Python version, every
third-party module the three steps import, the RTTOV installation with its
ECOSTRESS coefficient file and its `pyrttov` wrapper, the GPP look-up table, the
data directories, whether the output mounts are writable, and which credentials
are configured — by name only, never their values. It exits non-zero on the
first missing piece, so it is usable in CI.

A read-only output mount is worth checking this way: otherwise it only surfaces
at the end of a run, after hours of processing.

### Reference case

The environment check cannot tell whether the installation *computes* the right
thing. For that, run the pipeline on the reference granule and compare against
the shipped statistics:

```bash
docker compose run download --date 2025-01-01 --orbit 36798
docker compose run eeh-pipeline --input-files /data/rad/ECOv002_L1B_RAD_36798_016_20250101T165758_0713_02.h5
docker compose run selftest --products
```

The reference case is orbit **36798, scene 016 of 2025-01-01**, processed over
the default domain. `tools/reference/reference_36798_20250101.json` holds, for
each dataset of the three products, its shape, its dtype, the fraction of valid
pixels and the minimum, maximum and mean over those pixels. Fill values
(`-9999` for the float datasets, `0` for the integer ones) are excluded, so the
statistics describe retrieved pixels rather than the swath's empty corners.

This granule reaches 36.96° S while the default domain stops at 35° S, so the run
logs one `exceeds the ERA5 grid … outside pixels use edge values` warning per
step. That is expected here and part of what the reference statistics describe —
not a sign of a failed installation.

```
LST      PASS   mean 14798.3569, 99.1% valid
ETD      PASS   mean 1.3311, 31.2% valid
GPPd     PASS   mean 5.6875, 9.1% valid
```

A deviation is reported with both values and the tolerance it exceeded:

```
[FAIL] ETD    mean 1.4102, expected 1.3311 (tolerance 1e-03 relative)
```

Defaults are `1e-3` relative on minimum, maximum and mean, and `5e-3` absolute
on the fraction of valid pixels; both are flags
(`--relative-tolerance`, `--valid-fraction-tolerance`). Tolerances rather than
exact equality, because RTTOV and the BLAS the container links against are free
to reassociate floating-point work; a different BLAS build shifts the last
digits without changing the retrieval.

To adopt a new reference after a deliberate algorithm change, run
`tools/selftest.py --products --write-reference` and commit the regenerated
file, so the delivered code and its reference always move together.

## Local development (without Docker)

### Prerequisites

- Python >= 3.10
- [uv](https://docs.astral.sh/uv/) (Python package manager)
- RTTOV 12.3 installed locally (for TES step)
- System GDAL (for STIC step)

### Setup

```bash
uv sync                       # install Python dependencies
cp .env.example .env          # create configuration
```

Edit `.env` to set your local paths:

```dotenv
RTTOV_INSTALLDIR=D:/Projects/EURANUS/rttov123   # your local RTTOV path
INPUT_RAD_DIR=./data/rad                         # adjust to your data location
# ... etc
```

### Run

```bash
# Full pipeline
uv run python run_pipeline.py

# Single orbit (auto-downloads if not found locally)
uv run python run_pipeline.py --orbit 36798

# Date range (downloads + processes all orbits in the range)
uv run python run_pipeline.py --start-date 2026-08-01 --end-date 2026-08-03 --cog

# Single step via orchestrator
uv run python run_pipeline.py --step tes
uv run python run_pipeline.py --step stic
uv run python run_pipeline.py --step gpp

# Or call each script directly
uv run python TES/TES_main.py --input-dir ./data/rad --geo-dir ./data/geo \
    --era5-dir ./data/era5 --output-dir ./data/output/tes \
    --rttov-installdir /path/to/rttov123

uv run python STIC/STIC_main.py --input-dir ./data/output/tes --geo-dir ./data/geo \
    --cloud-dir ./data/cloud --fcover-dir ./data/fcover \
    --albedo-mota-dir ./data/mota --lulc-dir ./data/lulc \
    --era5-dir ./data/era5 --output-dir ./data/output/stic

uv run python GPP/GPP_main.py --stic-dir ./data/output/stic --geo-dir ./data/geo \
    --cloud-dir ./data/cloud --parh-dir ./data/parh --lai-dir ./data/lai \
    --fvc-dir ./data/fcover --era5-dir ./data/era5 --oco2-dir ./data/oco2 \
    --glc30-dir ./data/glc30 --lut-file GPP/LookUpTable_LUE_HH_gsFULL-globe.csv \
    --output-dir ./data/output/gpp
```

## Running on custom data

`tools/download_sample.py` downloads data for a single orbit. For your own orbits, you need to provide the matching data yourself.

### 1. Prepare your data directories

Place data in the expected directory structure (or adjust paths in `.env`):

```
data/
├── rad/       # ECOSTRESS L1B_RAD .h5 files
├── geo/       # ECOSTRESS L1B_GEO .h5 files (matching orbit numbers)
├── cloud/     # ECOSTRESS L2_CLOUD .h5 files
├── era5/      # ERA5 hourly GRIB files (see naming convention below)
├── mota/      # MCD43C3 MODIS albedo .hdf files
├── fcover/    # FCOVER 300m .nc files
├── lai/       # LAI 300m .nc files
├── lulc/      # PROBAV_LC100 .tif (single global file, ~1.6 GB)
├── oco2/      # OCO-2 XCO2 .nc4 files (GPP only)
├── parh/      # PARin*.nc files (GPP only)
└── glc30/     # GLC_FCS30D *_Annual.tif tiles (GPP only)
```

**Global datasets** (download once, reuse for all orbits): LULC, GLC30, OCO-2.
**Per-scene datasets** (must match each orbit's date and location): ERA5, FCOVER, LAI, MOTA, PAR.

### 2. Run the pipeline

```bash
# Full pipeline (TES → STIC → GPP)
docker compose up --build

# Single step
docker compose run eeh-pipeline --step tes
docker compose run eeh-pipeline --step stic
docker compose run eeh-pipeline --step gpp
```

Each step processes **all matching files** in its input directory. GEO/CLOUD files are matched to input files by orbit number automatically.

### 3. Run on specific files

Each script also accepts `--input-files` for single-file processing:

```bash
# Process a specific file
docker compose run eeh-pipeline python TES/TES_main.py \
    --input-files /data/rad/ECOv002_L1B_RAD_37000_005_20250315T120000_0712_01.h5 \
    --geo-dir /data/geo --era5-dir /data/era5 \
    --output-dir /data/output/tes --rttov-installdir /opt/rttov123

# Run only STIC on existing TES output
docker compose run eeh-pipeline python STIC/STIC_main.py \
    --input-dir /data/output/tes --geo-dir /data/geo \
    --cloud-dir /data/cloud --fcover-dir /data/fcover \
    --albedo-mota-dir /data/mota --lulc-dir /data/lulc \
    --era5-dir /data/era5 --output-dir /data/output/stic
```

Use `--help` on any script for the full list of arguments:

```bash
docker compose run eeh-pipeline python TES/TES_main.py --help
docker compose run eeh-pipeline python STIC/STIC_main.py --help
docker compose run eeh-pipeline python GPP/GPP_main.py --help
```

## Data sources

| Dataset | Source | Placement | Scope |
|---------|--------|-----------|-------|
| ECOSTRESS L1B_RAD, L1B_GEO, L2_CLOUD | [NASA Earthdata](https://search.earthdata.nasa.gov/) | `data/rad/`, `data/geo/`, `data/cloud/` | per scene |
| ERA5 single-levels + pressure-levels | [Copernicus CDS](https://cds.climate.copernicus.eu/) | `data/era5/` | per scene |
| MCD43C3 MODIS albedo | [NASA Earthdata (LP DAAC)](https://search.earthdata.nasa.gov/) | `data/mota/` | per scene |
| FCOVER 300m 10-day | [CDSE](https://dataspace.copernicus.eu/) | `data/fcover/` | per scene |
| LAI 300m 10-day | [CDSE](https://dataspace.copernicus.eu/) | `data/lai/` | per scene |
| LULC (PROBAV LC100) | [Zenodo 3518038](https://zenodo.org/records/3518038) | `data/lulc/` | global (download once) |
| OCO-2 XCO2 | [NASA GES DISC](https://disc.gsfc.nasa.gov/) | `data/oco2/` | global (2018–2022) |
| GLC30 land-cover tiles | [Zenodo 8239305](https://zenodo.org/records/8239305) | `data/glc30/` | global (download tiles for your area) |
| PAR half-hourly (PARin) | [CM SAF SAFIRA](https://wui.cmsaf.eu) or [EUMETSAT Data Store](https://data.eumetsat.int/) | `data/parh/` | per scene |

### Processing domain

ERA5 is only retrieved over a bounding box, and that box defines the pipeline's
**processing domain**. The default is Europe and Africa:

| | latitude | longitude |
|---|---|---|
| default domain | −35° … 75° | −20° … 60° |

An ECOSTRESS orbit is a full circuit of the globe, so most of its granules fall
outside this box. Such a granule has no atmospheric profile to be corrected with:
the ERA5 crop collapses to a single grid cell, the bilinear interpolation weights
sum to zero, and RTTOV rejects the resulting profile with

```
fatal error in module rttov_check_profiles.F90: invalid 2m air temperature
```

Granules outside the domain are therefore **excluded at download time**, and TES
skips any that are already on disk (synced from S3 or copied by hand), reporting
them in its end-of-run summary:

```
TES summary: 4 processed, 12 outside the processing domain, 0 failed (out of 16 granule(s))
```

Because CMR indexes only the `L1B_GEO` collection spatially — a `bounding_box`
search against `L1B_RAD` matches nothing at all — the download filter resolves the
in-domain scenes from GEO and then matches RAD and CLOUD on the `{orbit}_{scene}`
key they share with it.

Every CMR query paginates (`CMR-Search-After`). ECOSTRESS acquires more than 270
granules a day and CMR returns them ordered by start time, so reading a single
page truncated each day at its own page size and made the late-evening orbits
invisible: on 2026-04-20 all 19 granules of orbit 44173 (21:57–23:16 UTC) were
dropped and the downloader reported `No L1B_RAD granules found`. `--list-orbits`
counted the same day as 200 granules in 12 orbits instead of 283 in 16.

To process a different region, set `EEH2_DOMAIN` to `"N,W,S,E"`; ERA5 is then
requested over that box and the granule filter follows it:

```bash
EEH2_DOMAIN="40,60,-10,150" docker compose run eeh-pipeline --start-date 2026-08-01 --end-date 2026-08-01
```

ERA5 filenames record the timestamp only, never the box the file was requested
for, so a `data/era5/` filled for one domain is reused as-is by a run over
another. Change `EEH2_DOMAIN` without clearing those files and every pixel
outside the old box takes the value of the nearest grid edge. Both
`ReadERA5Data.py` warn when a granule exceeds the grid they opened:

```
granule lat [47.65, 52.23] lon [4.08, 10.94] exceeds the ERA5 grid
lat [49.00, 52.00] lon [2.00, 7.00] — outside pixels use edge values
```

Delete `data/era5/` when changing the domain, or download to a separate
`--output-dir`.

`tools/download_sample.py --whole-orbit` fetches every granule of an orbit
regardless of the domain. The extra granules still cannot be processed unless
`EEH2_DOMAIN` covers them, so this is mainly useful for archiving.

### Re-download guards: what counts as "already have it"

Every ancillary download is skipped when the data is already on disk, but "some
file of this kind is present" is never a sufficient test — the sample data shipped
for the reference orbit (2025-01-01, orbit 36798) would then suppress every later
download, and the step that needs the data fails instead with a missing-key error.
Each product is therefore checked against what the processing step will actually
ask for:

| product | skipped when |
|---|---|
| FCOVER / LAI | **every** 10-day composite the date needs is present — one for STIC, two for GPP (see below) |
| PAR | a **readable** `PARin<YYYYMMDD>*.nc` for that date is present (see below) |
| GLC30 | every 5° longitude strip the footprint needs is present |
| OCO-2 | a file for that **month and day** is present (see below) |
| ERA5 | the hourly file for each granule hour is present |

**PAR is stored unpacked.** The EUMETSAT Data Store answers with a ZIP holding the
NetCDF plus `EOPMetadata.xml` and `manifest.xml`, and the NetCDF member carries the
same name as the archive. The downloader extracts that member and keeps only it, so
`data/parh/` holds real NetCDF files. A ZIP saved under a `.nc` name is not usable:
GPP selects the file by name and then fails with
`OSError: [Errno -51] NetCDF: Unknown file format`, one granule at a time. Any such
file already on disk -- from an earlier version, a manual copy or the S3 bucket -- is
unpacked in place on the next PAR run, and the skip guard now opens the file rather
than trusting its extension.

Only `PARin` (half-hourly instantaneous) is downloaded. GPP takes the overpass time
slice and the daily mean from that one file, so the `PARdm` and `PARmm` products of
the same collection are not needed.

**CLMS 10-day composites** are labelled by the *end* of their period: days 1-10 →
`YYYYMM10`, 11-20 → `YYYYMM20`, 21+ → `YYYYMM<last day of month>`.

**STIC and GPP do not ask for the same composites.** STIC's `Read_FVC` indexes its
product map with the composite that *contains* the date, so a composite from a
neighbouring period is of no use to it — hence `KeyError: '20260810'` for
2026-08-01 if the wrong dekad is on disk. GPP instead *interpolates* LAI and FCOVER
between the two composites that **bracket** the date, so it needs a second, earlier
one: on 2026-08-01 it reads `20260731` and `20260810`. Downloading only STIC's
composite leaves GPP raising `StopIteration` in `interpolate_LAI_value` on every
granule. On a composite boundary the two needs do not even overlap — 2026-08-10
sends STIC to `20260810` and GPP to `20260731` and `20260820`, three files for one
date. The downloader therefore fetches the union, and skips per composite rather
than per directory.

This is invisible on the reference orbit: 2025-01-01 is day-of-year 1, which falls
in GPP's `doy <= 10` special case, where it reads a single composite and never
interpolates.

When CLMS has not yet published a composite the date needs, the downloader warns and
fetches the others rather than substituting a neighbouring period — the consumers
select by filename, so a stand-in would only move the failure further from its
cause. The run still fails on that date until the real composite appears.

**OCO-2** covers 2018-08-09 to 2022-02-28 only. For any later date GPP's `read_oco2`
falls back to the same month and day in the nearest year it has, so the downloader
resolves that substitute date itself and fetches it: 2026-08-01 is served by
2021-08-01. The download window used to be a fixed fortnight around the reference
orbit's 2025-01-01, which meant every other month and day found no file at all and
GPP failed inside `wrap_to_ECOSTRESS` on a `None` dataset rather than reporting a
missing input.

**GLC30** tiles are selected from the footprint of the granules being processed,
so a different orbit needs different tiles. The footprint is taken from each GEO
granule's `StandardMetadata` bounding coordinates (no pixel data is read), unioned
over the granules of the requested orbit and date, and clipped **per granule** to
the processing domain — an orbit circling the globe would otherwise yield the whole
domain and request land cover for all of Europe and Africa. Downloaded ZIPs are
recorded in `glc30/.glc30_fetched_zips`, because GLC30 has no tile over open water:
a longitude strip the swath crosses at sea would otherwise look permanently missing
and its multi-GB ZIP would be re-fetched on every run.

On a first run there is no GEO granule on disk yet — the S3 sync runs before the
ECOSTRESS download — so the footprint is resolved from **CMR granule metadata**
instead, which carries the same bounding coordinates. Without that fallback the
tile filter was silently disabled and the sync fell back to the *global* tile set:
43 GB of central and eastern Asian land cover for an orbit over Europe, and still
running when it was stopped. The reference orbit never showed it, because its GEO
granules were already on disk. Both sources agree to 5e-8 deg and select the same
tiles. If the footprint cannot be resolved at all, GLC30 is **skipped** rather than
synced globally — an unknown footprint is not a reason to fetch hundreds of GB.

### Failure reporting

TES, STIC and GPP catch per-granule failures so one bad granule does not abort the
batch, which means **exit code 0 does not mean every granule succeeded**. Two things
make the outcome visible:

- each step writes a JSON log of per-granule tracebacks to its own output
  directory (`--error-log-dir`, wired by `run_pipeline.py` to `/output/<step>`; the
  default `.` is `/app`, which is not a mounted volume, so those logs used to be
  lost when the container exited);
- `run_pipeline.py` reports how many files a step actually wrote, and warns
  explicitly when a step exits cleanly having written none.

### ERA5 file naming convention

ERA5 files must be split into **hourly** GRIB files following this naming pattern:

```
era5_single_levels-YYYY_MM_DD_HH:00.grib
era5_37levels-YYYY_MM_DD_HH:00.grib
```

`tools/download_sample.py` handles this automatically: it issues one bulk CDS request
per day and splits the result into hourly files locally (each message is routed
on its own GRIB validity time, so no external tool is needed).

**Which hours are downloaded** is derived from the ECOSTRESS granules on disk
for the requested date — no hardcoded hours:

| Dataset | Hours fetched |
|---|---|
| `era5_37levels` (pressure-levels) | acquisition hour **H** and **H+1** of each granule (TES interpolates between the two) |
| `era5_single_levels` | the same **H**/**H+1**, plus **11:00, 12:00, 13:00** UTC for STIC's daily max Ta |

So an orbit acquired at 13:15 UTC pulls pressure-levels at 13h/14h and
single-levels at 11h/12h/13h/14h. Day rollover is handled (an acquisition at
23:40 also fetches 00:00 of the next day), and `--orbit` restricts the derived
hours to that orbit's granules. If no granule is on disk to derive from, the
script warns and falls back to the reference orbit's hours.

> ECOSTRESS must therefore be downloaded before ERA5 — `tools/download_sample.py`
> already orders it that way.

## Project structure

```
eeh-public/
├── README.md                 # this file
├── Dockerfile                # Docker image (Python + GDAL + RTTOV)
├── docker-compose.yml        # container orchestration
├── pyproject.toml            # Python dependencies (for uv/pip)
├── uv.lock                   # locked dependency versions
├── .env.example              # configuration template (credentials + paths)
├── run_pipeline.py           # pipeline orchestrator (TES → STIC → GPP [→ COG])
├── tools/                    # supporting scripts (also runnable standalone)
│   ├── download_sample.py    # sample data downloader
│   └── h5_to_cog.py          # H5 → Cloud Optimized GeoTIFF converter
├── TES/                      # Land Surface Temperature
│   ├── TES_main.py           # entry point
│   ├── AtmCorrection.py      # RTTOV atmospheric correction
│   ├── ReadECOL1BData.py     # ECOSTRESS L1B reader
│   ├── ReadERA5Data.py       # ERA5 reader
│   ├── TES_vec.py            # TES algorithm (vectorized)
│   └── README.md
├── STIC/                     # Evapotranspiration
│   ├── STIC_main.py          # entry point
│   ├── PySTIC.py             # STIC model core
│   ├── Radiation.py          # radiation calculations
│   ├── SMWetness.py          # soil moisture
│   ├── TOA_Radiance.py       # TOA radiance
│   ├── LUT.py                # ET daily upscaling look-up
│   ├── ReadECOSTRESSData.py  # ECOSTRESS L2 reader
│   ├── ReadAncillaryData.py  # FCOVER, MOTA, LULC reader
│   ├── ReadERA5Data.py       # ERA5 reader
│   └── README.md
├── GPP/                      # Gross Primary Productivity
│   ├── GPP_main.py                           # entry point
│   ├── ReadData_EEH_GPP_final_TH.py         # data reader
│   ├── LookUpTable_LUE_HH_gsFULL-globe.csv  # LUE coefficients
│   └── README.md
```

## Output files

| Step | Output pattern | Contents |
|------|---------------|----------|
| TES | `EEH2TES_L2_LSTE_<orbit>_<datetime>_0000_00.h5` | LST, emissivity (bands 2/4/5), broadband emissivity, QA |
| STIC | `EEH2STIC_L3_ET_<orbit>_<datetime>_0000_00.h5` | LE, H, G, Rn, gah, gsc, Ms, Mrz, ETD |
| GPP | `EEH2STIC_L3_GPP_<key>_0000_00.h5` | GPPd, WUEd, GPP_gs, GPP_LUE, GPPi, LUCC |

### COG output (with `--cog`)

When `--cog` is passed, H5 outputs are resampled onto a regular lat/lon grid (~67 m) and converted to Cloud Optimized GeoTIFF:

| Product | Output dir | Source | Bands |
|---------|-----------|--------|-------|
| TES | `output/cog/tes/` | TES H5 | LST, Emis2, Emis4, Emis5, BBE, qa, land_fraction, cloud_mask |
| STIC | `output/cog/stic/` | STIC H5 | LE, H, G, Rn, gah, gsc, Ms, Mrz, ETD |
| GPP | `output/cog/gpp/` | GPP H5 | GPPd, GPPi, GPP_gs, GPP_LUE, LUCC |
| WUE | `output/cog/wue/` | GPP H5 | WUEd |

COG files use EPSG:4326, DEFLATE compression, 512×512 tiles, and include overviews. The standalone converter can also be used directly:

```bash
python tools/h5_to_cog.py --eeh-dir ./output/tes --geo-dir ./data/geo --cloud-dir ./data/cloud --out-dir ./output/cog/tes
python tools/h5_to_cog.py --eeh-dir ./output/gpp --geo-dir ./data/geo --out-dir ./output/cog/wue --force-product WUE
```

## References

- **TES**: Gillespie et al. (1998). *A temperature and emissivity separation algorithm for ASTER images.* IEEE TGRS, 36(4), 1113-1126.
- **STIC**: Mallick et al. (2014, 2016, 2018); Hu & Mallick et al. (2023). *Evaluating European ECOSTRESS Hub Evapotranspiration Products.* WRR, 59(8).
- **GPP/WUE**: Lin & Mallick et al.; Hu et al. See ATBD on https://ecostresshub.eu/

## License

MIT — see individual source files for author details.

Authors: Tian Hu, Yoanne Didry — Luxembourg Institute of Science and Technology (LIST)
