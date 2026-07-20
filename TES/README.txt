# This code is to estimate land surface temperature and emissivity from the ECOSTRESS top-of-atmosphere radiance data using the Temperature and Emissivity Separation (TES) algorithm with the support of the RTTOV lib for atmospheric correction
# The original TES algorithm was developed by Gillespie et al. (1998). More information about the TES algorithm can be found in the paper below:
# A. Gillespie, S. Rokugawa, T. Matsunaga, J. S. Cothern, S. Hook and A. B. Kahle, "A temperature and emissivity separation algorithm for Advanced Spaceborne Thermal Emission and Reflection Radiometer (ASTER) images," in IEEE Transactions on Geoscience and Remote Sensing, vol. 36, no. 4, pp. 1113-1126, July 1998, doi: 10.1109/36.700995.


# Input parameters (Order must be the same as below. Example is given in input.txt)
1) filename of ECOSTRESS L1B_RAD
2) directory for ECOSTRESS L1B_GEO
3) directory for ERA5 data
4) directory for storing output HDF5 data


# For processing each image, all the required parameters should be on the same line
# Different Lines represent the processings for different images


# To run the code, type the following in command line
python TES_main.py input.txt


# More detailed information about this software can be found on the software description document on the EEH landing page: https://ecostresshub.eu/


# The following paper should be cited when using this code:
Tian Hu, Kaniska Mallick, Glynn C. Hulley, Lluís Perez Planells, Frank M. Göttsche, Martin Schlerf, Patrik Hitzelberger, Yoanne Didry, Zoltan Szantoi, Itziar Alonso, José A. Sobrino, Dražen Skoković, Jean-Louis Roujean, Gilles Boulet, Philippe Gamet, Simon Hook,
Continental-scale evaluation of three ECOSTRESS land surface temperature products over Europe and Africa: Temperature-based validation and cross-satellite comparison,
Remote Sensing of Environment, Volume 282, 2022, 113296, https://doi.org/10.1016/j.rse.2022.113296

