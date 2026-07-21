# Scope
This code is to estimate gross primary productivity (GPP) and water use efficiency (WUE) at the daily scale based on the estimated evapotranspiration and conductances from the ECOSTRESS and ancillary data.

# Model
The hybrid model was initially developed by Lin and Mallick et al . Hu et al. refined the hybrid model for further performance improvement.

More information about the model can be found in the ATBD of the developed model, available on the EEH landing page: https://ecostresshub.eu/

# Inputs
Input parameters (Order must be the same as below) are:
1) filename of ECOSTRESS ET data 
2) directory for storing output HDF5 data
3) map_cld
4) map_geo
5) directory for PAR data 
6) directory for LAI data
7) directory for FVC data
8) directory for ERA5 data
9) directory for OCO2 data 
10) directory for LULC data
11) directory for cloud mask data 
12) directory for ECOSTRESS L1B_GEO

# Usage
To run the code, type the following in command line:

python myProcessor.py input.txt

For processing each image, all the 12 parameters should be on the same line. Different lines represent the processings for different images.




