#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Script to run the STIC model for ET estimation

Created on March 18 2021
@author: Tian Hu at LIST

© 2022 – Luxembourg Institute of Science and Technology
Authors : Tian Hu (tian.hu@list.lu), Kaniska Mallick
Code licensed under MIT
SPDX-License-Identifier: MIT
"""

#Main function for the STIC model
import sys
import os
import numpy as np
import h5py
#from osgeo import gdal, osr
from ReadECOSTRESSData import *
from ReadAncillaryData import *
from ReadERA5Data import * 
from PySTIC import *

#import cache function
sys.path.insert(0,'../utils')
#sys.path.insert(0,'/root/algos/eeh/utils')
from S3_cache import cache_S3_with_pattern
from S3_cache import cache_S3_with_indexes

import traceback
import json 

#Parse the config file to obtain all the directories
def parse_input_config(config_file):
    data = np.loadtxt(config_file,dtype=bytes).astype(str)
    
    if data.ndim == 1:
        row = 1
        col = data.shape[0]
        data = data.reshape((row,col))
    elif data.ndim == 2:
        row = data.shape[0]
        col = data.shape[1]
    
    if col == 9:
        ready = 1
    else:
        ready = 0
        
    return (data,row,col,ready)  

def generate_raster_file(driver, filename, data, x_dim, y_dim,
                         geo_transform, proj_wkt,
                         no_data_value, data_type):
    '''
    Description:
        Creates a raster file on disk for the specified data, using the
        specified driver.

    Note: It is assumed that the driver supports setting of the no data
          value.
          It is the callers responsibility to fix it if it does not.

    Note: It is assumed that the caller specified the correct file
          extension in the filename parameter for the specfied driver.
    '''

    try:
        raster = driver.Create(filename, x_dim, y_dim, 1, data_type)

        raster.SetGeoTransform(geo_transform)
        raster.SetProjection(proj_wkt)
        raster.GetRasterBand(1).WriteArray(data)
        raster.GetRasterBand(1).SetNoDataValue(no_data_value)
        raster.FlushCache()

        # Cleanup memory
        del raster

    except Exception:
        print('Failed to generate file {}'.format(filename))
#def generate_raster_file

#Output the estimates to HDF5 file
def generate_hdf5_file(directory_output,ET,H,G,Rn,orbit_str,date_str,hour_str,min_str,sec_str):     
    base_filename = 'EEHSTIC_L3_ET_' + orbit_str + '_' + date_str + 'T' + hour_str + min_str + sec_str + '_0000_00.h5'
    filename = os.path.join(directory_output,base_filename)
    f_h5 = h5py.File(filename,'w')
    
    dset = f_h5.create_dataset('LE',data = np.float32(ET))
    dset.attrs['long_name'] = 'Latent heat flux'
    dset.attrs['units'] = 'W.m-2'
    dset.attrs['format'] = 'scaled'
    dset.attrs['coordsys'] = 'cartesian'
    dset.attrs['valid_range'] = np.array([0,2000])
    dset.attrs['fill_value'] = -9999
    dset.attrs['scale_factor'] = 1
    dset.attrs['add_offset'] = 0

    dset = f_h5.create_dataset('H',data = np.float32(H))
    dset.attrs['long_name'] = 'Sensible heat flux'
    dset.attrs['units'] = 'W.m-2'
    dset.attrs['format'] = 'scaled'
    dset.attrs['coordsys'] = 'cartesian'
    dset.attrs['valid_range'] = np.array([0,2000])
    dset.attrs['fill_value'] = -9999
    dset.attrs['scale_factor'] = 1
    dset.attrs['add_offset'] = 0

    dset = f_h5.create_dataset('G',data = np.float32(G))
    dset.attrs['long_name'] = 'Soil heat flux'
    dset.attrs['units'] = 'W.m-2'
    dset.attrs['format'] = 'scaled'
    dset.attrs['coordsys'] = 'cartesian'
    dset.attrs['valid_range'] = np.array([0,2000])
    dset.attrs['fill_value'] = -9999
    dset.attrs['scale_factor'] = 1
    dset.attrs['add_offset'] = 0

    dset = f_h5.create_dataset('Rn',data = np.float32(Rn))
    dset.attrs['long_name'] = 'Net radiation'
    dset.attrs['units'] = 'W.m-2'
    dset.attrs['format'] = 'scaled'
    dset.attrs['coordsys'] = 'cartesian'
    dset.attrs['valid_range'] = np.array([0,2000])
    dset.attrs['fill_value'] = -9999
    dset.attrs['scale_factor'] = 1
    dset.attrs['add_offset'] = 0

    f_h5.close() 

    print("STIC Model Running Completed")
    return
    
def run_STIC_from_config_file(config_file):
    (config_data,row,col,ready) = parse_input_config(config_file)

    #Cache the GEO and CLOUD file name
    #map_geo = cache_S3_with_indexes('L1B_GEO',0,43)
    map_geo = cache_S3_with_pattern('L1B_GEO_V002','/ECOSTRESS_RW/',r'(ECOv002_L1B_GEO.*)_\d{4}_\d{2}.h5')
    #map_fvc = cache_S3_with_pattern('FCOVER')
    map_fvc = cache_S3_with_pattern('FCOVER','/ECOSTRESS_RW/',r'c_gls_FCOVER300(?:-RT\d+)?_(\d{8})0000_GLOBE_.*\.nc')
    #map_albdir1 = cache_S3_with_indexes('Albedo_Directional',0,18)
    map_albdir2 = cache_S3_with_pattern('MOTA','/ECOSTRESS_RW/',r'MCD43C3\.A(\d{7}).*.hdf')
    #map_albdir2 = cache_S3_with_indexes('MOTA',9,16)
   # map_albhem1 = cache_S3_with_indexes('Albedo_Hemispherical',0,18)
    #map_albhem2 = cache_S3_with_indexes('MOTA',9,16)
    map_albhem2 = map_albdir2
    
    print('Done caching ')
    #done cashing

    map_error = dict()

    if ready == 1:
        for i in range(row):
            filename_lste = config_data[i,0]
            directory_geo = config_data[i,1]
            directory_cld = config_data[i,2]
            directory_fvc = config_data[i,3]
            directory_alb_dir = config_data[i,4]
            directory_alb_hem = config_data[i,5]
            directory_lulc = config_data[i,6]
            directory_era5 = config_data[i,7]
            directory_output = config_data[i,8]

            print('Run STIC model on ECOSTRESS data ' + filename_lste)           
            
            # Read the ECOSTRESS data
            (lst_eco,lse_eco,lat_eco,lon_eco,_,watermask,cloudmask,year,month,day,hour,minute,second,
             date_str,hour_str,min_str,sec_str,orbit_str,_) = Read_ECOSTRESS(filename_lste,directory_geo,directory_cld,map_geo)
            
            # Converting from K to Celsus degree
            lst_eco = lst_eco - 273.15
            
            # Converting time from UTM to local solar time (in seconds)
            time_delta = lon_eco/15.
            time_ls = float(hour) + float(minute)/60. + float(second)/3600. + time_delta
            time_ls[time_ls < 0] += 24
            time_ls[time_ls >= 24] = time_ls[time_ls >= 24] % 24
            time_ls = time_ls*3600            
            
            # Mask out pixels covered by cloud
            lst_eco[cloudmask == 1] = np.nan
            lse_eco[cloudmask == 1] = np.nan
            
            # Mask out pixels covered by water
            lst_eco[watermask == 1] = np.nan
            lse_eco[watermask == 1] = np.nan
            
            print('Reading ECOSTRESS data finished!')
            
            # Read and interpolate the ancillary data based on ECOSTRESS geolocation and observation time
            (fvc,alb_dir,alb_hem,lulc) = Read_Ancillary(directory_fvc,directory_alb_dir,directory_alb_hem,directory_lulc,
                                           lat_eco,lon_eco,year,month,day,map_fvc,map_albdir1,map_albdir2,map_albhem1,map_albhem2)
            
            print('Reading ancillary CGLS data finished!')
            
            (t_s,rh_s,_,sr_dir,sr_dif) = Read_ERA5(directory_era5,lat_eco,lon_eco,year,month,day,hour,minute,second)
            
            t_s = t_s - 273.15 # Converting from K to Celsus degree
            rh_s = rh_s*100 # Converting from 0-1 to percentage            
               
            print('Reading ERA5 data finished!')
            
            try:
                (RN_STIC,G_STIC,H_STIC,LE_STIC,_,gsc_STIC,
                        _,Ms_STIC,Mrz_STIC,converged) = STIC(lst_eco,lse_eco,t_s,rh_s,sr_dir,
                                                                sr_dif,alb_dir,alb_hem,fvc,time_ls)

                print('Running the STIC model finished!')
                
    #             #Extract the samples for testing the MATLAB version of STIC
    #             lst_eco1 = lst_eco[converged]
    #             lse_eco1 = lse_eco[converged]
    #             t_s1 = t_s[converged]
    #             rh_s1 = rh_s[converged]
    #             sr_dir1 = sr_dir[converged]
    #             sr_dif1 = sr_dif[converged]
    #             alb_dir1 = alb_dir[converged]
    #             alb_hem1 = alb_hem[converged]
    #             fvc1 = fvc[converged]
    #             time_ls1 = time_ls[converged]
    #             RN_STIC1 = RN_STIC[converged]
    #             G_STIC1 = G_STIC[converged]
    #             H_STIC1 = H_STIC[converged]
    #             LE_STIC1 = LE_STIC[converged]
                
    #             index = np.random.randint(lst_eco1.shape[0],size=(500))
    #             data = np.transpose(np.array([lst_eco1[index],lse_eco1[index],t_s1[index],rh_s1[index],sr_dir1[index],
    #                                 sr_dif1[index],alb_dir1[index],alb_hem1[index],fvc1[index],time_ls1[index],
    #                                          RN_STIC1[index],G_STIC1[index],H_STIC1[index],LE_STIC1[index]]))
    #             outputname = '/home/tian/ecostress/data/temp/STIC/data.txt'
    #             np.savetxt(outputname,data,fmt='%12.6f')          
                
                #Mask out the invalid pixels
                RN_STIC[np.isnan(RN_STIC)] = -9999
                G_STIC[np.isnan(G_STIC)] = -9999
                H_STIC[np.isnan(H_STIC)] = -9999
                LE_STIC[np.isnan(LE_STIC)] = -9999
                
                #Mask out the non-converged pixels
                RN_STIC[~converged] = -9999
                G_STIC[~converged] = -9999
                H_STIC[~converged] = -9999
                LE_STIC[~converged] = -9999
                
                #Mask out the night-time pixels
                G_STIC[RN_STIC < 0] = -9999
                H_STIC[RN_STIC < 0] = -9999
                LE_STIC[RN_STIC < 0] = -9999
                RN_STIC[RN_STIC < 0] = -9999
                
                #Mask out ET for certain land surface types
                LE_STIC[(lulc == 50) | (lulc == 60) | (lulc == 70)] = -9999
                H_STIC[(lulc == 50) | (lulc == 60) | (lulc == 70)] = -9999
                G_STIC[(lulc == 50) | (lulc == 60) | (lulc == 70)] = -9999
                RN_STIC[(lulc == 50) | (lulc == 60) | (lulc == 70)] = -9999
                
                # Output the variables in HDF5 file
                generate_hdf5_file(directory_output,LE_STIC,H_STIC,G_STIC,RN_STIC,orbit_str,
                                date_str,hour_str,min_str,sec_str)
                
                print('Outputting to HDF5 file finished')
            except Exception as err:
                print(str(err))
                traceback.print_exc()
                print('Error -> stopping there')
                map_error[filename_lste] = traceback.format_exc()
    else:
        message = 'ERROR!! Check the input directories!'
        raise Exception(message)
    
    if map_error:
        pid = str(os.getpid())
        with open('execution_'+pid+'_errors.txt', 'w') as file:
            file.write(json.dumps(map_error))
 
if __name__ == '__main__':
    args = sys.argv
    if len(args) > 1:
        config_file = args[1]
        run_STIC_from_config_file(config_file)
    else:
        message = 'ERROR!! Configuration data required!'
        raise Exception(message) 
