# Scope
This code is to estimate gross primary productivity (GPP) and water use efficiency (WUE) at the daily scale based on the estimated evapotranspiration and conductances from the ECOSTRESS and ancillary data.

# Model
The hybrid model was initially developed by Lin and Mallick et al . Hu et al. refined the hybrid model for further performance improvement.

More information about the model can be found in the ATBD of the developed model, available on the EEH landing page: https://ecostresshub.eu/

# Inputs
The command-line arguments are listed as below:
1) `--stic-files` / `--stic-dir` : EEH2STIC L3 ET/STIC file(s) to process (STIC output)
2) `--stic-pattern` : glob pattern to select STIC files (only used with `--stic-dir`, default: `*EEH2STIC*.h5`)
3) `--output-dir`  : directory for storing output GPP/WUE HDF5 data
4) `--geo-dir`     : directory for ECOSTRESS L1B_GEO
5) `--cloud-dir`   : directory for ECOSTRESS cloud mask
6) `--parh-dir`    : directory for hourly PAR data (PARin)
7) `--lai-dir`     : directory for LAI data (300m, 10-day)
8) `--fvc-dir`     : directory for FCOVER data (300m, 10-day)
9) `--era5-dir`    : directory for ERA5 data (single level)
10) `--oco2-dir`   : directory for OCO-2 XCO2 data
11) `--glc30-dir`  : directory for GLC30 land cover tiles
12) `--lut-file`   : path to the LookUpTable_LUE_HH_gsFULL-globe.csv file
13) `--geo-pattern` / `--cloud-pattern` : regex patterns (one capturing group, anchored with `$`) for GEO/CLOUD filenames
14) `--parh-pattern` / `--lai-pattern` / `--fvc-pattern` / `--era5-pattern` / `--oco2-pattern` / `--glc30-pattern` : glob patterns for the corresponding ancillary files (default: `*`)
15) `--cache-dir`  : directory to persist file-name lookup maps (JSON cache), avoids rescanning slow/networked storage (default: `<output-dir>/.gpp_filename_cache`; pass `""` to disable)
16) `--force-rebuild-cache` : ignore any existing cache and rescan all directories
17) `--error-log-dir` : directory for the execution error JSON log, if any (default: current directory)

# Usage
To run the code, type the following in command line:

python GPP_main.py \
    --stic-files /data/output_stic/EEH2STIC_L3_ET_..._0000_00.h5 \
    --geo-dir /data/geo \
    --cloud-dir /data/cloud \
    --parh-dir /data/parh \
    --lai-dir /data/lai \
    --fvc-dir /data/fvc \
    --era5-dir /data/era5 \
    --oco2-dir /data/oco2 \
    --glc30-dir /data/glc30 \
    --lut-file /data/LookUpTable_LUE_HH_gsFULL-globe.csv \
    --output-dir /data/output_gpp

To process a whole batch of files instead of a single file, replace
`--stic-files` with:

    --stic-dir /data/output_stic \
    --stic-pattern "*EEH2STIC_L3_ET*.h5"

If your local GEO/CLOUD filenames do not follow the standard ECOSTRESS
naming convention, override the matching patterns, e.g.:

    --geo-pattern '(ECOv002_L1B_GEO.*)_\d{4}_\d{2}\.h5$' \
    --cloud-pattern '(ECOv002_L2_CLOUD.*)_\d{4}_\d{2}\.h5$'

# Optional
--force-rebuild-cache : ignore the file-name cache and force a rescan of
                         all 8 ancillary data directories


--output-dir /data/output_gpp
