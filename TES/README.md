# Scope
This code is to estimate land surface temperature and emissivity from the ECOSTRESS top-of-atmosphere radiance data using the Temperature and Emissivity Separation (TES) algorithm. 

# Critical dependency
The RTTOV lib (version 13.2, https://nwp-saf.eumetsat.int/site/software/rttov/) is used for atmospheric correction. Users need to install the lib before using this code.

# Algorithm
The original TES algorithm was developed by Gillespie et al. (1998). More information about the TES algorithm can be found in the paper below:

A. Gillespie, S. Rokugawa, T. Matsunaga, J. S. Cothern, S. Hook and A. B. Kahle, "A temperature and emissivity separation algorithm for Advanced Spaceborne Thermal Emission and Reflection Radiometer (ASTER) images," in IEEE Transactions on Geoscience and Remote Sensing, vol. 36, no. 4, pp. 1113-1126, July 1998, doi: 10.1109/36.700995.

# Inputs
The command-line arguments are listed as below:
The command-line arguments are listed as below:
1) `--input-files` / `--input-dir` : ECOSTRESS L1B_RAD (V002) file(s) to process
2) `--pattern`     : glob pattern to select RAD files (only used with `--input-dir`, default: `*L1B_RAD*.h5`)
3) `--geo-dir`     : directory for ECOSTRESS L1B_GEO (V002)
4) `--era5-dir`    : directory for ERA5 data
5) `--output-dir`  : directory for storing output HDF5 data
6) `--rttov-installdir` (+ optional `--rttov-wrapper-dir`, `--rttov-coef-file`) : local RTTOV installation
7) `--cache-dir`   : directory to persist the GEO file-name lookup map (JSON cache), avoids rescanning slow/networked storage (default: `<output-dir>/.tes_filename_cache`; pass `""` to disable caching)
8) `--force-rebuild-cache` : ignore any existing GEO file-name cache and rescan the directory
9) `--error-log-dir` : directory for the execution error JSON log, if any (default: current directory)

# Usage
To run the code, type the following in command line:

python TES_main.py \
    --input-files /data/ECOv002_L1B_RAD_..._0000_00.h5 /data/ECOv002_L1B_RAD_..._0000_01.h5 \
    --geo-dir /data/geo \
    --era5-dir /data/era5 \
    --output-dir /data/output \
    --rttov-installdir /root/rttov \
    --rttov-wrapper-dir /root/rttov/wrapper

To process a whole batch of files instead of listing them individually,
replace `--input-files` with:

    --input-dir /data/rad \
    --pattern "ECOv002_L1B_RAD_*.h5"

# Optional
--alpha1 / --alpha2 / --alpha3   : override the default TES coefficients
                                    (0.9895, 0.7994, 0.8572)
--rttov-lib-preload               : advanced workaround only, preloads a
                                    compiled RTTOV f2py wrapper .so via
                                    ctypes before importing pyrttov
--force-rebuild-cache             : ignore the GEO file-name cache and
                                    force a directory rescan


# Citation
The following paper should be cited when using this code:

Tian Hu, Kaniska Mallick, Glynn C. Hulley, Lluís Perez Planells, Frank M. Göttsche, Martin Schlerf, Patrik Hitzelberger, Yoanne Didry, Zoltan Szantoi, Itziar Alonso, José A. Sobrino, Dražen Skoković, Jean-Louis Roujean, Gilles Boulet, Philippe Gamet, Simon Hook, Continental-scale evaluation of three ECOSTRESS land surface temperature products over Europe and Africa: Temperature-based validation and cross-satellite comparison, Remote Sensing of Environment, Volume 282, 2022, 113296, https://doi.org/10.1016/j.rse.2022.113296

