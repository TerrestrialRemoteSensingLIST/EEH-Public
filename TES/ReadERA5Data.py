#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Script to run the Temperature Emissivity Separation (TES) algorithm for LST estimation

Created on April 15 2021
@author: Tian Hu at LIST

© 2026 – Luxembourg Institute of Science and Technology
Authors : Tian Hu (tian.hu@list.lu)
SPDX-License-Identifier: MIT
"""

# Read the ERA5 data

import cfgrib
import numpy as np
import os
from datetime import datetime,timedelta

# Linear interpolation function
def linearinterp(time, variables):
    y1 = variables[0]
    y2 = variables[1]
    outp = time/60*(y2-y1) + y1
    return outp

# Bilinear interpolation function
def bilinearinterp(x,y,variable,dimension):
    if dimension == 2:     
        x0 = np.floor(x).astype(int)
        x1 = x0 + 1
        y0 = np.floor(y).astype(int)
        y1 = y0 + 1

        x0 = np.clip(x0, 0, variable.shape[1]-1)
        x1 = np.clip(x1, 0, variable.shape[1]-1)
        y0 = np.clip(y0, 0, variable.shape[0]-1)
        y1 = np.clip(y1, 0, variable.shape[0]-1)
        
        Ia = variable[y0,x0]
        Ib = variable[y1,x0]
        Ic = variable[y0,x1]
        Id = variable[y1,x1]
        
        wa = (x1-x) * (y1-y)
        wb = (x1-x) * (y-y0)
        wc = (x-x0) * (y1-y)
        wd = (x-x0) * (y-y0)

        return wa*Ia + wb*Ib + wc*Ic + wd*Id
    
    if dimension == 3:
        
        x0 = np.floor(x).astype(int)
        x1 = x0 + 1
        y0 = np.floor(y).astype(int)
        y1 = y0 + 1

        x0 = np.clip(x0, 0, variable.shape[2]-1)
        x1 = np.clip(x1, 0, variable.shape[2]-1)
        y0 = np.clip(y0, 0, variable.shape[1]-1)
        y1 = np.clip(y1, 0, variable.shape[1]-1)
        
        levels = variable.shape[0]
        outp = np.empty((levels,x.shape[0],x.shape[1]),dtype=np.float64)

        for i in range(levels):
            Ia = variable[i,y0,x0]
            Ib = variable[i,y1,x0]
            Ic = variable[i,y0,x1]
            Id = variable[i,y1,x1]

            wa = (x1-x) * (y1-y)
            wb = (x1-x) * (y-y0)
            wc = (x-x0) * (y1-y)
            wd = (x-x0) * (y-y0)
            
            outp[i,:,:] = wa*Ia + wb*Ib + wc*Ic + wd*Id

        return outp
    return None


def Read_ERA5(directory_era5,lat_eco,lon_eco,year,month,day,hour,minute,second):
    #Construct ERA5 filenames accordingly
    
    time = float(minute) + float(second)/60
    time1_str = year + month + day + hour  
#     time1_str = date_str + hour_str
    
    time2 = datetime.strptime(time1_str,'%Y%m%d%H') + timedelta(hours=1)
    time2_str = time2.strftime('%Y%m%d%H')

    filename1 = time1_str[0:4] + '_' + time1_str[4:6] + '_' + time1_str[6:8] + '_' + time1_str[8:]
    filename2 = time2_str[0:4] + '_' + time2_str[4:6] + '_' + time2_str[6:8] + '_' + time2_str[8:]

    era5_filename = []
    era5_filename.append(filename1)
    era5_filename.append(filename2)

    #Crop the global ERA5 data to the spatial range of the whole ECOSTRESS image
    #ERA5 data is from lat: 75 to -35, from lon: -20 to 60
    max_lat = np.amax(lat_eco)
    min_lat = np.amin(lat_eco)
    max_lon = np.amax(lon_eco)
    min_lon = np.amin(lon_eco)

    max_lat = np.clip(max_lat, -35, 75)   # Latitudes entre -35 et 75
    min_lat = np.clip(min_lat, -35, 75)
    max_lon = np.clip(max_lon, -20, 60)   # Longitudes entre -20 et 60
    min_lon = np.clip(min_lon, -20, 60)

    row1 = np.floor((75 - max_lat)/0.25).astype(int)
    row2 = np.ceil((75 - min_lat)/0.25).astype(int)
    col1 = np.floor((min_lon + 20)/0.25).astype(int)
    col2 = np.ceil((max_lon + 20)/0.25).astype(int)

    # p_era5
    # lat_era5
    # lon_era5

    t_array = []
    q_array = []
    
    sp_array = []
    q2m_array = []
    t2m_array = []
    skt_array = []
    for i in range(2):
        #Atmospheric profiles at 37 levels
        base_filename = 'era5_37levels-' + era5_filename[i] + ':00.grib'
        print('Reading ERA5 file '+base_filename)
        filename = os.path.join(directory_era5,base_filename)

        ds_t = cfgrib.open_dataset(filename,engine='cfgrib',backend_kwargs={'filter_by_keys': {'shortName': 't'}, 'indexpath': ''})
        t = ds_t.t.data[:,row1:row2+1,col1:col2+1]
        t_array.append(t)
        
        p_era5 = ds_t.isobaricInhPa.data
        lat_era5 = ds_t.latitude.data[row1:row2+1]
        lon_era5 = ds_t.longitude.data[col1:col2+1]  
        ds_t.close()
        del ds_t,t    

        ds_q = cfgrib.open_dataset(filename,engine='cfgrib',backend_kwargs={'filter_by_keys': {'shortName': 'q'}, 'indexpath': ''})
        q = ds_q.q.data[:,row1:row2+1,col1:col2+1]
        q_array.append(q)
        ds_q.close()
        del ds_q,q
        
        #Atmospheric profiles at the surface level
        base_filename1 = 'era5_single_levels-' + era5_filename[i] + ':00.grib'
        print('Reading ERA5 file '+base_filename1)
        filename1 = os.path.join(directory_era5,base_filename1)

        #Data from the single level profile
        ds_sp = cfgrib.open_dataset(filename1,engine='cfgrib',backend_kwargs={'filter_by_keys': {'shortName': 'sp'}, 'indexpath': ''})
        sp = ds_sp.sp.data[row1:row2+1,col1:col2+1]/100.0 # converting from Pa to hPa
        sp_array.append(sp) 
    
        ds_d2m = cfgrib.open_dataset(filename1,engine='cfgrib',backend_kwargs={'filter_by_keys': {'shortName': '2d'}, 'indexpath': ''})
        d2m = ds_d2m.d2m.data[row1:row2+1,col1:col2+1]
        e0 = 6.113
        c_water = 5423
        t0 = 273.15
        q2m = (0.622 * e0 * np.exp(c_water * (d2m - t0)/(d2m * t0)))/sp # dewpoint temperature to specific humidity
        q2m_array.append(q2m)
        ds_sp.close()
        ds_d2m.close()
        del q2m,d2m,ds_d2m, ds_sp

        ds_t2m = cfgrib.open_dataset(filename1,engine='cfgrib',backend_kwargs={'filter_by_keys': {'shortName': '2t'}, 'indexpath': ''})
        t2m = ds_t2m.t2m.data[row1:row2+1,col1:col2+1]
        t2m_array.append(t2m)
        ds_t2m.close()
        del t2m,ds_t2m

        ds_skt = cfgrib.open_dataset(filename1,engine='cfgrib',backend_kwargs={'filter_by_keys': {'shortName': 'skt'}, 'indexpath': ''})
        skt = ds_skt.skt.data[row1:row2+1,col1:col2+1]
        skt_array.append(skt)
        ds_skt.close()
        del skt,ds_skt       
        
    tempinterp_t = linearinterp(time,t_array)
    tempinterp_q = linearinterp(time,q_array)
    tempinterp_sp = linearinterp(time,sp_array)
    tempinterp_q2m = linearinterp(time,q2m_array)
    tempinterp_t2m = linearinterp(time,t2m_array)
    tempinterp_skt = linearinterp(time,skt_array)

    del t_array
    del q_array
    del sp_array
    del q2m_array
    del t2m_array
    del skt_array

    #FIXME uncomment later
    lon_eco = np.clip(lon_eco, np.min(lon_era5)+1e-6, np.max(lon_era5)-1e-6)
    lat_eco = np.clip(lat_eco, np.min(lat_era5)+1e-6, np.max(lat_era5)-1e-6)

    x = (lon_eco - lon_era5[0]) / 0.25
    y = (lat_era5[0] - lat_eco) / 0.25

    spatiotempinterp_t = bilinearinterp(x,y,tempinterp_t,3)
    spatiotempinterp_q = bilinearinterp(x,y,tempinterp_q,3)
    spatiotempinterp_sp = bilinearinterp(x,y,tempinterp_sp,2)
    spatiotempinterp_q2m = bilinearinterp(x,y,tempinterp_q2m,2)
    spatiotempinterp_t2m = bilinearinterp(x,y,tempinterp_t2m,2)
    spatiotempinterp_skt = bilinearinterp(x,y,tempinterp_skt,2)

    del tempinterp_t
    del tempinterp_q
    del tempinterp_sp
    del tempinterp_q2m
    del tempinterp_t2m
    del tempinterp_skt

    ################Flatten the profiles###############
    new_t = spatiotempinterp_t.reshape(spatiotempinterp_t.shape[0],(spatiotempinterp_t.shape[1]*spatiotempinterp_t.shape[2]))
    new_q = spatiotempinterp_q.reshape(spatiotempinterp_q.shape[0],(spatiotempinterp_q.shape[1]*spatiotempinterp_q.shape[2]))
    new_q2m = spatiotempinterp_q2m.reshape(spatiotempinterp_q2m.shape[0]*spatiotempinterp_q2m.shape[1])
    new_t2m = spatiotempinterp_t2m.reshape(spatiotempinterp_t2m.shape[0]*spatiotempinterp_t2m.shape[1])
    new_skt = spatiotempinterp_skt.reshape(spatiotempinterp_skt.shape[0]*spatiotempinterp_skt.shape[1])

    del spatiotempinterp_t
    del spatiotempinterp_q
    del spatiotempinterp_q2m
    del spatiotempinterp_t2m
    del spatiotempinterp_skt

    return (p_era5, new_t, new_q, new_q2m, new_t2m, new_skt,spatiotempinterp_sp)
