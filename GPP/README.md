# Scope
This code is to estimate gross primary productivity (GPP) and water use efficiency (WUE) at the daily scale based on the estimated evapotranspiration and conductances from the ECOSTRESS and ancillary data.

# Model
The hybrid model was initially developed by Lin and Mallick et al . Hu et al. refined the hybrid model for further performance improvement.

More information about the model can be found in the ATBD of the developed model, available on the EEH landing page: https://ecostresshub.eu/

#Inputs
The command-line arguments are listed as below:

--stic-files / --stic-dir : EEH2STIC L3 ET/STIC file(s) to process (STIC output, e.g. EEH2STIC_L3_ET_..._0000_00.h5)
--output-dir : directory for storing output GPP/WUE HDF5 data
--geo-dir : directory for ECOSTRESS L1B_GEO
--cloud-dir : directory for ECOSTRESS cloud mask
--parh-dir : directory for hourly PAR data (PARin)
--lai-dir : directory for LAI data (300m, 10-day)
--fvc-dir : directory for FCOVER data (300m, 10-day)
--era5-dir : directory for ERA5 data (single level)
--oco2-dir : directory for OCO-2 XCO2 data
--glc30-dir : directory for GLC30 land cover tiles
--lut-file : path to the LookUpTable_LUE_HH_gsFULL-globe.csv file
--geo-pattern (optional) : regex pattern (one capturing group) used to match ECOSTRESS L1B_GEO filenames. Default: (ECOv002_L1B_GEO.*)_\d{4}_\d{2}\.h5
--cloud-pattern (optional) : regex pattern (one capturing group) used to match ECOSTRESS L2_CLOUD filenames. Default: (ECOv002_L2_CLOUD.*)_\d{4}_\d{2}\.h5
--error-log-dir (optional) : directory to write the execution error JSON log to, if any. Default: current directory

#Usage
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


To process a whole batch of files (instead of a single file), replace --stic-files with:
--stic-dir /data/output_stic \
--stic-pattern "*EEH2STIC_L3_ET*.h5"


--output-dir /data/output_gpp
