#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Script to run the STIC model for ET estimation

Created on March 18 2021
@author: Tian Hu at LIST

"""
import sys
import os
import numpy as np
import h5py
from datetime import datetime
#from osgeo import gdal, osr
from ReadECOSTRESSData import *
from ReadAncillaryData import *
from ReadERA5Data import * 
from PySTIC import *
from TOA_Radiance import *
from LUT import *

from S3_cache import cache_S3_with_pattern

import traceback
import json 

import wasdi

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
        wasdi.wasdiLog('Failed to generate file {}'.format(filename))
#def generate_raster_file

def adjust_lulc(lulc):
    lulc[(lulc == 80) | (lulc == 200)] = 0 # water
    lulc[(lulc == 111) | (lulc == 121)] = 1 # evergreen conifer
    lulc[(lulc == 112) | (lulc == 122)] = 1 # evergreen broadleaf
    lulc[(lulc == 113) | (lulc == 123)] = 1 # deciduous conifer
    lulc[(lulc == 114) | (lulc == 124)] = 1 # deciduous broadleaf
    lulc[(lulc == 115) | (lulc == 116) | (lulc == 125) | (lulc == 126)] = 1 # mixed forest
    lulc[lulc == 20] = 6 # woody savanna
    lulc[(lulc == 30) | (lulc == 100)] = 3 # grassland
    lulc[lulc == 90] = 5 # wetland
    lulc[lulc == 40] = 2 # cropland
    lulc[lulc == 50] = 8 # urban
    lulc[lulc == 70] = 9 # snow
    lulc[lulc == 60] = 7 # barren surfaces
    lulc[lulc == 0] = 0 # exclusion of invalid pixels

    return lulc

#Output the estimates to HDF5 file
def generate_hdf5_file(directory_output,ET,H,G,Rn,gah,gsc,Ms,Mrz,ETD,orbit_str,date_str,hour_str,min_str,sec_str):     
    base_filename = 'EEH2STIC_L3_ET_' + orbit_str + '_' + date_str + 'T' + hour_str + min_str + sec_str + '_0000_00.h5'
    filename = os.path.join(directory_output,base_filename)
    wasdi.wasdiLog("Writing output to : "+filename)
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

    dset = f_h5.create_dataset('gah',data = np.float32(gah))
    dset.attrs['long_name'] = 'Aerodynamic conductance'
    dset.attrs['units'] = 'm.s-1'
    dset.attrs['format'] = 'scaled'
    dset.attrs['coordsys'] = 'cartesian'
    dset.attrs['valid_range'] = np.array([0,1])
    dset.attrs['fill_value'] = -9999
    dset.attrs['scale_factor'] = 1
    dset.attrs['add_offset'] = 0

    dset = f_h5.create_dataset('gsc',data = np.float32(gsc))
    dset.attrs['long_name'] = 'Surface conductance'
    dset.attrs['units'] = 'm.s-1'
    dset.attrs['format'] = 'scaled'
    dset.attrs['coordsys'] = 'cartesian'
    dset.attrs['valid_range'] = np.array([0,1])
    dset.attrs['fill_value'] = -9999
    dset.attrs['scale_factor'] = 1
    dset.attrs['add_offset'] = 0

    dset = f_h5.create_dataset('Ms',data = np.float32(Ms))
    dset.attrs['long_name'] = 'Water stress'
    dset.attrs['units'] = 'unitless'
    dset.attrs['format'] = 'scaled'
    dset.attrs['coordsys'] = 'cartesian'
    dset.attrs['valid_range'] = np.array([0,1])
    dset.attrs['fill_value'] = -9999
    dset.attrs['scale_factor'] = 1
    dset.attrs['add_offset'] = 0

    dset = f_h5.create_dataset('Mrz',data = np.float32(Mrz))
    dset.attrs['long_name'] = 'Water stress root zone'
    dset.attrs['units'] = 'unitless'
    dset.attrs['format'] = 'scaled'
    dset.attrs['coordsys'] = 'cartesian'
    dset.attrs['valid_range'] = np.array([0,1])
    dset.attrs['fill_value'] = -9999
    dset.attrs['scale_factor'] = 1
    dset.attrs['add_offset'] = 0

    dset = f_h5.create_dataset('ETD',data = np.float32(ETD))
    dset.attrs['long_name'] = 'Daily ET'
    dset.attrs['units'] = 'mm.day-1'
    dset.attrs['format'] = 'scaled'
    dset.attrs['coordsys'] = 'cartesian'
    dset.attrs['valid_range'] = np.array([0,2000])
    dset.attrs['fill_value'] = -9999
    dset.attrs['scale_factor'] = 1
    dset.attrs['add_offset'] = 0

    f_h5.close() 

    wasdi.wasdiLog("STIC Model Running Completed")
    
    return base_filename
    
def run_STIC_from_wasdi_config_file():
 
    #For local testing only
    #wasdi.getPath('ecostress_fast/L1B_GEO_V002/ECOv002_L1B_GEO_00478_013_20180806T133126_0712_04.h5')
    #wasdi.getPath('ecostress_fast/L2_CLOUD_V002/ECOv002_L2_CLOUD_00478_013_20180806T133126_0712_04.h5')
    #wasdi.getPath('ecostress_fast/FCOVER/c_gls_FCOVER300_201800600000_GLOBE_PROBAV_V1.0.1.nc')
    #wasdi.getPath('ecostress_fast/Albedo_Directional/c_gls_ALDH_201808060000_GLOBE_PROBAV_V1.5.1.nc')
    #wasdi.getPath('ecostress_fast/Albedo_Hemispherical/c_gls_ALBH_201808060000_GLOBE_PROBAV_V1.5.1.nc')
    #wasdi.getPath('ecostress_fast/LULC/PROBAV_LC100_global_v3.0.1_2018-conso_Discrete-Classification-map_EPSG-4326.tif')
    #wasdi.getPath('ecostress_fast/ERA5/era5_single_levels-2018_08_06_14:00.grib')
    #wasdi.getPath('ecostress_fast/ERA5/era5_single_levels-2018_08_06_15:00.grib')
    #wasdi.getPath('ecostress_fast/ERA5/era5_single_levels-2018_08_06_11:00.grib')
    ###wasdi.getPath('ecostress_fast/ERA5/era5_single_levels-2018_08_06_12:00.grib')
    #wasdi.getPath('ecostress_fast/MOTA/MCD43C3.A2018218.061.2021347183125.hdf')


    #Cache the GEO and CLOUD file name from the working directory
    #ECOv002_L1B_GEO_00424_001_20180803T011700_0712_04.h5
    map_geo = cache_S3_with_pattern('L1B_GEO',wasdi.getSavePath()+'ecostress_fast/',r'(ECOv002_L1B_GEO.*)_\d{4}_\d{2}.h5')
    
    #TODO : better handling
    map_cloud   = cache_S3_with_pattern('CLOUD',wasdi.getSavePath()+'ecostress_fast/',r'(ECOv002_L2_CLOUD.*)_\d{4}_\d{2}.h5')
    map_fvc = cache_S3_with_pattern('FCOVER',wasdi.getSavePath()+'ecostress_fast/',r'c_gls_FCOVER300(?:-RT\d+)?_(\d{8})0000_GLOBE_.*\.nc')
    #c_gls_FCOVER300_201908200000_GLOBE_PROBAV_V1.0.1.nc or c_gls_FCOVER300-RT6_202407310000_GLOBE_OLCI_V1.1.2.nc
    map_albdir1 = {}
    map_albhem1 = {}
    map_albdir2 = cache_S3_with_pattern('MOTA',wasdi.getSavePath()+'ecostress_fast/',r'MCD43C3\.A(\d{7}).*.hdf')
    map_albhem2 = map_albdir2

    wasdi.wasdiLog('Done caching ')
    #done cashing

    map_error = dict()

    lste_files = wasdi.getParameter('L2_LSTE_FILES') #EEHTES folder i.e. ecostress_fast/EEH2/EEHTES/
    row = len(lste_files)

    wasdi.updateStatus("RUNNING", 0)
    for i in range(row):

        filename_lste = wasdi.getPath(lste_files[i])

        directory_geo = wasdi.getSavePath()+wasdi.getParameter('L1B_GEO_PATH')
        directory_cld = wasdi.getSavePath()+wasdi.getParameter('L2_CLOUD_PATH')
        directory_fvc = wasdi.getSavePath()+wasdi.getParameter('FCOVER_PATH')
        directory_alb_dir = wasdi.getSavePath()+wasdi.getParameter('ALBEDODIR_PATH')
        directory_alb_hem = wasdi.getSavePath()+wasdi.getParameter('ALBEDOHEM_PATH')
        directory_lulc = wasdi.getSavePath()+wasdi.getParameter('LULC_PATH')
        directory_era5 = wasdi.getSavePath()+wasdi.getParameter('ERA5_PATH')
        directory_output = wasdi.getSavePath()+wasdi.getParameter('STIC_OUTPUT')

        wasdi.wasdiLog('Run STIC model on ECOSTRESS data ' + filename_lste)           
        
        # Read the ECOSTRESS data
        (lst_eco,lse_eco,lat_eco,lon_eco,_,watermask,cloudmask,year,month,day,hour,minute,second,
         date_str,hour_str,min_str,sec_str,orbit_str,_) = Read_ECOSTRESS(filename_lste,directory_geo,directory_cld,map_geo,map_cloud)
        
        # Converting from K to Celsus degree
        lst_eco = lst_eco - 273.15

        # Calculating day of year and decimal time 
        #date = datetime(float(year),float(month),float(day))
        date = datetime(int(year),int(month),int(day))
        doy = date.timetuple().tm_yday
        time_decimal = float(hour) + float(minute)/60. + float(second)/3600
        
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
        
        wasdi.wasdiLog('Reading ECOSTRESS data finished!')
        
        # Read and interpolate the ancillary data based on ECOSTRESS geolocation and observation time
        (fvc,alb_dir,alb_hem,lulc) = Read_Ancillary(directory_fvc,directory_alb_dir,directory_alb_hem,directory_lulc,
                                       lat_eco,lon_eco,year,month,day,map_fvc,map_albdir1,map_albdir2,map_albhem1,map_albhem2)
        lulc_adjusted = adjust_lulc(lulc)
        wasdi.wasdiLog('Reading ancillary CGLS data finished!')
        
        (t_s,rh_s,_,sr_dir,sr_dif,ta_max) = Read_ERA5(directory_era5,lat_eco,lon_eco,year,month,day,hour,minute,second)
        
        t_s = t_s - 273.15 # Converting from K to Celsus degree
        rh_s = rh_s*100 # Converting from 0-1 to percentage            
           
        wasdi.wasdiLog('Reading ERA5 data finished!')
        
        try:
            (RN_STIC,G_STIC,H_STIC,LE_STIC,gah_STIC,gsc_STIC,
                    T0_STIC,Ms_STIC,Mrz_STIC,converged,Lin,Lout,Lnet,Sin,Sout) = STIC(lst_eco,lse_eco,t_s,rh_s,sr_dir,
                                                            sr_dif,alb_dir,alb_hem,fvc,time_ls)

            wasdi.wasdiLog('Running the STIC model finished!')

            #Mask out the invalid pixels
            RN_STIC[np.isnan(RN_STIC)] = -9999
            G_STIC[np.isnan(G_STIC)] = -9999
            H_STIC[np.isnan(H_STIC)] = -9999
            LE_STIC[np.isnan(LE_STIC)] = -9999
            gah_STIC[np.isnan(gah_STIC)] = -9999
            gsc_STIC[np.isnan(gsc_STIC)] = -9999
            Ms_STIC[np.isnan(Ms_STIC)] = -9999
            Mrz_STIC[np.isnan(Mrz_STIC)] = -9999
            
            #Mask out the non-converged pixels
            RN_STIC[~converged] = -9999
            G_STIC[~converged] = -9999
            H_STIC[~converged] = -9999
            LE_STIC[~converged] = -9999
            gah_STIC[~converged] = -9999
            gsc_STIC[~converged] = -9999
            Ms_STIC[~converged] = -9999
            Mrz_STIC[~converged] = -9999
            
            #Mask out the night-time pixels
            G_STIC[RN_STIC < 0] = -9999
            H_STIC[RN_STIC < 0] = -9999
            LE_STIC[RN_STIC < 0] = -9999
            RN_STIC[RN_STIC < 0] = -9999
            gah_STIC[RN_STIC < 0] = -9999
            gsc_STIC[RN_STIC < 0] = -9999
            Ms_STIC[RN_STIC < 0] = -9999
            Mrz_STIC[RN_STIC < 0] = -9999
            
            #Mask out ET for non-vegetated land surface types
            LE_STIC[(lulc_adjusted < 1) | (lulc_adjusted > 6)] = -9999
            H_STIC[(lulc_adjusted < 1) | (lulc_adjusted > 6)] = -9999
            G_STIC[(lulc_adjusted < 1) | (lulc_adjusted > 6)] = -9999
            RN_STIC[(lulc_adjusted < 1) | (lulc_adjusted > 6)] = -9999
            gah_STIC[(lulc_adjusted < 1) | (lulc_adjusted > 6)] = -9999
            gsc_STIC[(lulc_adjusted < 1) | (lulc_adjusted > 6)] = -9999
            Ms_STIC[(lulc_adjusted < 1) | (lulc_adjusted > 6)] = -9999
            Mrz_STIC[(lulc_adjusted < 1) | (lulc_adjusted > 6)] = -9999

            #Calculate daily ET from the instantaneous estimtes
            (RgTOAiMJ,RgTOAInt) = f_TOARadiance(doy,time_decimal,lat_eco,lon_eco)
            LE_STIC1 = LE_STIC * 30 * 60 / (2264.76*1e3) #Unit conversion: mm
            ET_DLY = f_ETDaily(LE_STIC1,float(hour),t_s,ta_max,Lin,Lout,Sin,Sout,RgTOAiMJ,RgTOAInt,lulc_adjusted)
            ET_DLY[LE_STIC == -9999] = -9999
            
            # Output the variables in HDF5 file
            generate_hdf5_file(directory_output,LE_STIC,H_STIC,G_STIC,RN_STIC,gah_STIC,gsc_STIC,Ms_STIC,Mrz_STIC,ET_DLY,
                               orbit_str,date_str,hour_str,min_str,sec_str)
            
            wasdi.wasdiLog('Outputting to HDF5 file finished')
        except Exception as err:
            traceback.print_exc()
            wasdi.wasdiLog('Error -> stopping there')
            map_error[filename_lste] = traceback.format_exc()
    
    if map_error:
        pid = str(os.getpid())
        with open('execution_'+pid+'_errors.txt', 'w') as file:
            file.write(json.dumps(map_error))
 
def run():
    wasdi.wasdiLog("Starting STIC algorithm")
    run_STIC_from_wasdi_config_file()
    wasdi.updateStatus("DONE",100)
    wasdi.wasdiLog("STIC algorithm done")

if __name__ == '__main__':
    wasdi.init("./config.json")
    run()
