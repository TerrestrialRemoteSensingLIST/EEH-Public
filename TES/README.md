# Scope
This code is to estimate land surface temperature and emissivity from the ECOSTRESS top-of-atmosphere radiance data using the Temperature and Emissivity Separation (TES) algorithm. 

# Critical dependency
The RTTOV lib (version 13.2, https://nwp-saf.eumetsat.int/site/software/rttov/) is used for atmospheric correction. Users need to install the lib before using this code.

# Algorithm
The original TES algorithm was developed by Gillespie et al. (1998). More information about the TES algorithm can be found in the paper below:

A. Gillespie, S. Rokugawa, T. Matsunaga, J. S. Cothern, S. Hook and A. B. Kahle, "A temperature and emissivity separation algorithm for Advanced Spaceborne Thermal Emission and Reflection Radiometer (ASTER) images," in IEEE Transactions on Geoscience and Remote Sensing, vol. 36, no. 4, pp. 1113-1126, July 1998, doi: 10.1109/36.700995.

# Inputs
This step used to require a fixed 4-column input.txt config file (RAD
filename / GEO directory / ERA5 directory / output directory). This has
been replaced by command-line arguments (see Usage below):

1) `--input-files` / `--input-dir` : ECOSTRESS L1B_RAD (V002) file(s) to process
2) `--geo-dir`     : directory for ECOSTRESS L1B_GEO (V002)
3) `--era5-dir`    : directory for ERA5 data
4) `--output-dir`  : directory for storing output HDF5 data
5) `--rttov-installdir` (+ optional `--rttov-wrapper-dir`, `--rttov-coef-file`) : local RTTOV installation

# Usage
To run the code, type the following in command line:

python TES_main.py \
    --input-files /data/ECOv002_L1B_RAD_..._0000_00.h5 /data/ECOv002_L1B_RAD_..._0000_01.h5 \
    --geo-dir /data/geo \
    --era5-dir /data/era5 \
    --output-dir /data/output \
    --rttov-installdir /root/rttov \
    --rttov-wrapper-dir /root/rttov/wrapper


# Citation
The following paper should be cited when using this code:

Tian Hu, Kaniska Mallick, Glynn C. Hulley, Lluís Perez Planells, Frank M. Göttsche, Martin Schlerf, Patrik Hitzelberger, Yoanne Didry, Zoltan Szantoi, Itziar Alonso, José A. Sobrino, Dražen Skoković, Jean-Louis Roujean, Gilles Boulet, Philippe Gamet, Simon Hook, Continental-scale evaluation of three ECOSTRESS land surface temperature products over Europe and Africa: Temperature-based validation and cross-satellite comparison, Remote Sensing of Environment, Volume 282, 2022, 113296, https://doi.org/10.1016/j.rse.2022.113296

