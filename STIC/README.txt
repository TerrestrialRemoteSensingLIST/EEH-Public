# Input parameters (Order must be the same as below. Example is given in input.txt)
1) filename of ECOSTRESS LST data
2) directory for ECOSTRESS L1B_GEO
3) directory for ECOSTRESS cloud mask
4) directory ofr FVC data
5) directory for blacksky albedo
6) directory for whitesky albedo
7) directory for LULC data
8) directory for ERA5 data
9) directory for storing output HDF5 data


# For processing each image, all the 9 parameters should be on the same line
# Different Lines represent the processings for different images


# To run the code, type the following in command line
python STIC_main.py input.txt
