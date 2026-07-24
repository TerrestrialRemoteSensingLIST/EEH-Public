#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Script to run the Temperature Emissivity Separation (TES) algorithm for LST estimation

Created on April 15 2021
@author: Tian Hu at LIST

Modified to run fully locally (no S3 access required):
 - Removed cache_S3_with_indexes / S3_cache dependency.
 - Added build_local_geo_map(), which scans a local directory of GEO (.h5)
   files and builds the same {canonical_name: actual_filename} mapping that
   Read_L1B_Data() expects, so no code changes were needed downstream.
"""
import sys
import os
import numpy as np
import h5py
import cv2
from ReadECOL1BData import Read_L1B_Data
from ReadERA5Data import Read_ERA5
from AtmCorrection import runRTTOV
from TES_vec import LST_Estimate
from scipy.interpolate import griddata

#rttov manual insertion from path
#import ctypes
#lib = ctypes.CDLL('/volume/rttov/lib/rttov_wrapper_f2py.so')

import traceback
import json


def build_local_geo_map(directory_geo, key_index1=0, key_index2=41,
                         exclude_substrings=('.xml', '.dmrpp')):
    """
    Local, pure-Python replacement for cache_S3_with_indexes().

    Original S3 version (S3_cache.py) worked by:
      1. `ls -U <mount_folder><folder>` (excluding .xml / .dmrpp files),
      2. keying each resulting filename by a fixed character slice
         thing[key_index1:key_index2],
      3. mapping that key -> the full filename.

    Read_L1B_Data() builds its lookup key as:
        'ECOv002_L1B_GEO_' (16 chars) + orbit_str (9) + '_' (1)
        + date_str (8) + 'T' (1) + hour_str+min_str+sec_str (6)
      = 41 characters total, i.e. filename[0:41].

    This function reproduces that exact keying, but scans a local
    directory directly instead of shelling out to `ls` against an S3
    mount, and takes directory_geo as an argument instead of a
    hardcoded mount_folder/folder pair.

    Parameters
    ----------
    directory_geo : str
        Local folder containing the GEO .h5 files.
    key_index1, key_index2 : int
        Character slice [key_index1:key_index2] of the filename used as
        the lookup key. Defaults (0, 41) match the existing GEO naming
        convention and Read_L1B_Data()'s key construction. Only change
        these if your local filenames differ from that convention.
    exclude_substrings : tuple of str
        Files containing any of these substrings are skipped (mirrors
        the original `grep -v xml | grep -v dmrpp` filtering).

    Returns
    -------
    dict : {filename[key_index1:key_index2]: actual_filename}
    """
    if not os.path.isdir(directory_geo):
        raise FileNotFoundError(f"GEO directory not found: {directory_geo}")

    geo_map = {}

    for fname in sorted(os.listdir(directory_geo)):
        if any(sub in fname for sub in exclude_substrings):
            continue
        if not os.path.isfile(os.path.join(directory_geo, fname)):
            continue

        key = fname[key_index1:key_index2]

        if key in geo_map:
            print(f'Warning: duplicate GEO key "{key}" -> keeping "{geo_map[key]}", '
                  f'ignoring "{fname}"')
            continue

        geo_map[key] = fname

    if not geo_map:
        print(f'Warning: no GEO files found in {directory_geo}')

    return geo_map


def parse_input_config(config_file):
    data = np.loadtxt(config_file,dtype=bytes).astype(str)
    
    if data.ndim == 1:
        row = 1
        col = data.shape[0]
        data = data.reshape((row,col))
    elif data.ndim == 2:
        row = data.shape[0]
        col = data.shape[1]
    
    if col == 4:
        ready = 1
    else:
        ready = 0
        
    return (data,row,col,ready)
    

def generate_hdf5_file_V2(directory_output, lst, emib2, emib4, emib5, bbe, mask, qa, orbit_str, date_str, hour_str, min_str, sec_str, version):
    lst_scaled = np.uint16(lst/0.02)
    emi_b2_scaled = np.uint8((emib2 - 0.489999)/0.002)
    emi_b4_scaled = np.uint8((emib4 - 0.489999)/0.002)
    emi_b5_scaled = np.uint8((emib5 - 0.489999)/0.002)
    bbe_scaled = np.uint8((bbe - 0.489999)/0.002)

    lst_scaled[mask] = 0
    emi_b2_scaled[mask] = 0
    emi_b4_scaled[mask] = 0
    emi_b5_scaled[mask] = 0
    bbe_scaled[mask] = 0

    base_filename = 'EEHTES_L2_LSTE_' + version + '_'+ orbit_str + '_' + date_str + 'T' + hour_str + min_str + sec_str + '_0000_00.h5'
    filename = os.path.join(directory_output,base_filename)

    with h5py.File(filename, 'w') as f_lst:
        # Adding GZIP compression with compression level 4
        dset_names = ['Emis2', 'Emis4', 'Emis5', 'BBE', 'LST', 'qa']
        datasets = {
            'Emis2': emi_b2_scaled,
            'Emis4': emi_b4_scaled,
            'Emis5': emi_b5_scaled,
            'BBE': bbe_scaled,
            'LST': lst_scaled,
            'qa': np.int16(qa)
        }
        attributes = {
            'Emis2': {'long_name': 'Band 2 Emissivity', 'units': 'n/a', 'format': 'scaled', 'coordsys': 'cartesian', 'valid_range': [1, 255], 'fill_value': 0, 'scale_factor': 0.002, 'add_offset': 0.489999},
            'Emis4': {'long_name': 'Band 4 Emissivity', 'units': 'n/a', 'format': 'scaled', 'coordsys': 'cartesian', 'valid_range': [1, 255], 'fill_value': 0, 'scale_factor': 0.002, 'add_offset': 0.489999},
            'Emis5': {'long_name': 'Band 5 Emissivity', 'units': 'n/a', 'format': 'scaled', 'coordsys': 'cartesian', 'valid_range': [1, 255], 'fill_value': 0, 'scale_factor': 0.002, 'add_offset': 0.489999},
            'BBE': {'long_name': 'Broad Band Emissivity', 'units': 'n/a', 'format': 'scaled', 'coordsys': 'cartesian', 'valid_range': [1, 255], 'fill_value': 0, 'scale_factor': 0.002, 'add_offset': 0.489999},
            'LST': {'long_name': 'Land Surface Temperature', 'units': 'K', 'format': 'scaled', 'coordsys': 'cartesian', 'valid_range': [7500, 65535], 'fill_value': 0, 'scale_factor': 0.02, 'add_offset': 0},
            'qa': {'long_name': 'Quality Assurance', 'units': 'N/A', 'format': 'scaled', 'coordsys': 'cartesian', 'valid_range': [-5, 5], 'fill_value': -9999, 'scale_factor': 1, 'add_offset': 0}
        }

        for name in dset_names:
            dset = f_lst.create_dataset(name, data=datasets[name], compression='gzip', compression_opts=4)
            for key, value in attributes[name].items():
                dset.attrs[key] = value

def generate_hdf5_file(directory_output, lst, emib2, emib4, emib5, bbe, mask, qa, orbit_str, date_str, hour_str, min_str, sec_str, version):  
    
    lst_scaled = np.uint16(lst/0.02)       
    emi_b2_scaled = np.uint8((emib2 - 0.489999)/0.002)
    emi_b4_scaled = np.uint8((emib4 - 0.489999)/0.002)
    emi_b5_scaled = np.uint8((emib5 - 0.489999)/0.002)
    bbe_scaled = np.uint8((bbe - 0.489999)/0.002)
    
    lst_scaled[mask] = 0
    emi_b2_scaled[mask] = 0
    emi_b4_scaled[mask] = 0
    emi_b5_scaled[mask] = 0
    bbe_scaled[mask] = 0
    
    base_filename = 'EEHTES_L2_LSTE_' + version + '_'+ orbit_str + '_' + date_str + 'T' + hour_str + min_str + sec_str + '_0000_00.h5'
    filename = os.path.join(directory_output,base_filename)
    f_lst = h5py.File(filename,'w')
    
    dset = f_lst.create_dataset('Emis2',data = emi_b2_scaled)
    dset.attrs['long_name'] = 'Band 2 Emissivity'
    dset.attrs['units'] = 'n/a'
    dset.attrs['format'] = 'scaled'
    dset.attrs['coordsys'] = 'cartesian'
    dset.attrs['valid_range'] = np.array([1,255])
    dset.attrs['fill_value'] = 0
    dset.attrs['scale_factor'] = 0.002
    dset.attrs['add_offset'] = 0.489999
    
    dset = f_lst.create_dataset('Emis4',data = emi_b4_scaled)
    dset.attrs['long_name'] = 'Band 4 Emissivity'
    dset.attrs['units'] = 'n/a'
    dset.attrs['format'] = 'scaled'
    dset.attrs['coordsys'] = 'cartesian'
    dset.attrs['valid_range'] = np.array([1,255])
    dset.attrs['fill_value'] = 0
    dset.attrs['scale_factor'] = 0.002
    dset.attrs['add_offset'] = 0.489999
    
    dset = f_lst.create_dataset('Emis5',data = emi_b5_scaled)
    dset.attrs['long_name'] = 'Band 5 Emissivity'
    dset.attrs['units'] = 'n/a'
    dset.attrs['format'] = 'scaled'
    dset.attrs['coordsys'] = 'cartesian'
    dset.attrs['valid_range'] = np.array([1,255])
    dset.attrs['fill_value'] = 0
    dset.attrs['scale_factor'] = 0.002
    dset.attrs['add_offset'] = 0.489999
    
    dset = f_lst.create_dataset('BBE',data = bbe_scaled)
    dset.attrs['long_name'] = 'Broad Band Emissivity'
    dset.attrs['units'] = 'n/a'
    dset.attrs['format'] = 'scaled'
    dset.attrs['coordsys'] = 'cartesian'
    dset.attrs['valid_range'] = np.array([1,255])
    dset.attrs['fill_value'] = 0
    dset.attrs['scale_factor'] = 0.002
    dset.attrs['add_offset'] = 0.489999
    
    dset = f_lst.create_dataset('LST',data = lst_scaled)
    dset.attrs['long_name'] = 'Land Surface Temperature'
    dset.attrs['units'] = 'K'
    dset.attrs['format'] = 'scaled'
    dset.attrs['coordsys'] = 'cartesian'
    dset.attrs['valid_range'] = np.array([7500,65535])
    dset.attrs['fill_value'] = 0
    dset.attrs['scale_factor'] = 0.02
    dset.attrs['add_offset'] = 0
    
    dset = f_lst.create_dataset('qa',data = np.int16(qa))
    dset.attrs['long_name'] = 'Quality Assurance'
    dset.attrs['units'] = 'N/A'
    dset.attrs['format'] = 'scaled'
    dset.attrs['coordsys'] = 'cartesian'
    dset.attrs['valid_range'] = np.array([-5,5])
    dset.attrs['fill_value'] = -9999
    dset.attrs['scale_factor'] = 1
    dset.attrs['add_offset'] = 0
    
    f_lst.close()   
    

def run_TES_from_config_file(config_file):
    (config_data,row,_,ready) = parse_input_config(config_file)

    map_error = dict()
    geo_map_cache = {}  # directory_geo -> map_geo, so we don't rescan the same folder repeatedly

    if ready == 1:
        for i in range(row):           
            filename_rad = config_data[i,0]
            directory_geo = config_data[i,1]
            directory_era5 = config_data[i,2]
            directory_output = config_data[i,3]

            if directory_geo not in geo_map_cache:
                print(f'Building local GEO file index for {directory_geo}')
                geo_map_cache[directory_geo] = build_local_geo_map(directory_geo)
            map_geo = geo_map_cache[directory_geo]

            print('Run Temperature Emissivity Separation algorithm on ECOSTRESS data ' + filename_rad)

            try:
                (lat_eco,lon_eco,alt_eco,lf,_,vza,sza,year,month,day,hour,minute,second,
                r2,rqa2,r4,rqa4,r5,rqa5,
                 date_str,hour_str,min_str,sec_str,orbit_str) = Read_L1B_Data(filename_rad,directory_geo,map_geo)           
            except Exception as e:
                print(str(e))
                traceback.print_exc()
                print('Error while opening L1B Data -> stopping there')
                map_error[filename_rad] = traceback.format_exc()
                continue

            print('Reading ECOSTRESS L1B data completed')
            
            # Resize the ECOSTRESS data with a 10-by-10 window
            nsubset = 10
            nrow = np.ceil(lat_eco.shape[0]/nsubset).astype(int)
            ncol = np.ceil(lat_eco.shape[1]/nsubset).astype(int)
            
            lat_eco_agg = cv2.resize(lat_eco, dsize=(ncol, nrow), interpolation=cv2.INTER_LINEAR)
            lon_eco_agg = cv2.resize(lon_eco, dsize=(ncol, nrow), interpolation=cv2.INTER_LINEAR)
            alt_eco_agg = cv2.resize(alt_eco, dsize=(ncol, nrow), interpolation=cv2.INTER_LINEAR)
            lf_agg = cv2.resize(lf, dsize=(ncol, nrow), interpolation=cv2.INTER_LINEAR)           
            vza_agg = cv2.resize(vza, dsize=(ncol, nrow), interpolation=cv2.INTER_LINEAR)
            sza_agg = cv2.resize(sza, dsize=(ncol, nrow), interpolation=cv2.INTER_LINEAR)
#             time_agg = cv2.resize(time, dsize=(ncol, nrow), interpolation=cv2.INTER_LINEAR)
            watermask_agg = np.ones((nrow, ncol))
            watermask_agg[lf_agg > 0.5] = 0 # Land 
            
            print('Aggregating ECOSTRESS L1B data completed')
            
            try:
                #print('Directory ERA5='+directory_era5)
                #print('LAT_ECO='+str(lat_eco_agg))
                #print('LON_ECO='+str(lon_eco_agg))
                #print('year='+str(year))
                #print('month='+str(month))
                #print('day='+str(day))
                #print('hour='+str(hour))
                #print('minute='+str(minute))
                #print('second='+str(second))
                (p_era5, new_t, new_q, new_q2m, new_t2m, new_skt,_)=Read_ERA5(directory_era5,lat_eco_agg,lon_eco_agg,year,month,day,hour,minute,second)
                new_q[new_q < 0.1e-10] = 0.1e-10
            except Exception as e:
                print(str(e))
                traceback.print_exc()
                print('Error while opening ERA5 Data -> stopping there')
                map_error[filename_rad] = traceback.format_exc()
                continue
            
            print('Reading and interpolating ERA5 data completed')
            
            # Calculate surface pressure based on the ECOSTRESS altitude data
            sp = 1013.25 * (1.0 - 2.225577e-5 * alt_eco_agg)**5.25588 # Altitude unit: m
            new_sp = sp.reshape(sp.shape[0]*sp.shape[1])

            new_t2m

            try:
                (trans0, upclear0, dnclear0) = runRTTOV(p_era5, new_t, new_q, new_sp, new_q2m, new_t2m, 
                                     new_skt, vza_agg, sza_agg, watermask_agg,
                                     lat_eco_agg, lon_eco_agg, alt_eco_agg, year, month, day)
            except Exception as e:
                print(str(e))
                traceback.print_exc()
                print('Error while running RTTOV -> stopping there')
                map_error[filename_rad] = traceback.format_exc()
                continue

            
            print('Running RTTOV completed')
            
            # Convert from 1-dimension to 3-dimension
            upclear = upclear0.reshape(lat_eco_agg.shape[0],lat_eco_agg.shape[1],3)
            dnclear = dnclear0.reshape(lat_eco_agg.shape[0],lat_eco_agg.shape[1],3)
            trans = trans0.reshape(lat_eco_agg.shape[0],lat_eco_agg.shape[1],3)
            
            # Interpolate the atmospheric parameters
            upclear_f = np.zeros((lat_eco.shape[0],lat_eco.shape[1],3))
            dnclear_f = np.zeros((lat_eco.shape[0],lat_eco.shape[1],3))
            trans_f = np.zeros((lat_eco.shape[0],lat_eco.shape[1],3))
            
            xy = np.vstack((np.ravel(lon_eco_agg), np.ravel(lat_eco_agg))).T
            
            z = np.ravel(upclear[:,:,0])
            upclear_f[:,:,0] = griddata(xy,z,(lon_eco,lat_eco),method='nearest')
            z = np.ravel(upclear[:,:,1])
            upclear_f[:,:,1] = griddata(xy,z,(lon_eco,lat_eco),method='nearest')
            z = np.ravel(upclear[:,:,2])
            upclear_f[:,:,2] = griddata(xy,z,(lon_eco,lat_eco),method='nearest')
            
            print('Interpolating RTTOV upL output completed')
            
            z = np.ravel(dnclear[:,:,0])
            dnclear_f[:,:,0] = griddata(xy,z,(lon_eco,lat_eco),method='nearest')
            z = np.ravel(dnclear[:,:,1])
            dnclear_f[:,:,1] = griddata(xy,z,(lon_eco,lat_eco),method='nearest')
            z = np.ravel(dnclear[:,:,2])
            dnclear_f[:,:,2] = griddata(xy,z,(lon_eco,lat_eco),method='nearest')
            
            print('Interpolating RTTOV dnL output completed')
            
            z = np.ravel(trans[:,:,0])
            trans_f[:,:,0] = griddata(xy,z,(lon_eco,lat_eco),method='nearest')
            z = np.ravel(trans[:,:,1])
            trans_f[:,:,1] = griddata(xy,z,(lon_eco,lat_eco),method='nearest')
            z = np.ravel(trans[:,:,2])
            trans_f[:,:,2] = griddata(xy,z,(lon_eco,lat_eco),method='nearest')
            
            print('Interpolating RTTOV trans output completed')

            coefficients_sets = [
                (0.9895, 0.7994, 0.8572),  # Coefficients for EEH TES
                (0.9692, 0.8117, 0.9957),  # Coefficients for input samples for SAIL271
                (0.9824, 0.8931, 0.9757)   # Coefficients for SAIL271
            ]

            version_name = ['A1','A2','A3']
            v = 0
            for alpha1, alpha2, alpha3 in coefficients_sets:
                try:
                    version = version_name[v]
                    (lst, emib2, emib4, emib5, _, qa) = LST_Estimate(r2, r4, r5, upclear_f, dnclear_f, trans_f, alpha1, alpha2, alpha3)
                    
                    bbe = np.zeros(emib2.shape)
                    mask = np.logical_and.reduce((emib2 >0, emib4 >0, emib5 >0))
                    bbe[mask] = 0.3287*emib2[mask] + 0.3783*emib4[mask] + 0.3158*emib5[mask] - 0.0255
                    
                    # Creating mask for the invaild values
                    mask = np.logical_or.reduce((rqa2 > 1,rqa4 > 1,rqa5 > 1, lst == 0, bbe == 0))
                    
                    print('Estimating LST completed')

                    generate_hdf5_file_V2(directory_output, lst, emib2, emib4, emib5, bbe, mask, qa, 
                                    orbit_str, date_str, hour_str, min_str, sec_str, version)
                    
                    print('Outputting to HDF5 completed')
                    v = v+1
                except Exception as err:
                    print(str(err))
                    traceback.print_exc()
                    print('Error -> stopping there')
                    map_error[filename_rad] = traceback.format_exc()
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
        run_TES_from_config_file(config_file)
    else:
        message = 'ERROR!! Configuration Data Required!'
        raise Exception(message)