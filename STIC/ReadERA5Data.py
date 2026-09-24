#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Script to run the STIC model for ET estimation

Created on September 1 2024
@author: Tian Hu, Kaniska Mallick, Yoanne Didry at LIST

© 2026 – Luxembourg Institute of Science and Technology
Authors : Tian Hu (tian.hu@list.lu), Yoanne Didry (yoanne.didry@list.lu)
Code licensed under MIT
SPDX-License-Identifier: MIT
"""

#Functions for reading the ERA5 data

import cfgrib
import logging
import numpy as np
import os
from datetime import datetime,timedelta

logger = logging.getLogger(__name__)

#Linear interpolation function
def _grid_span(axis, lo, hi):
    """Index range of `axis` covering [lo, hi], clamped to the axis itself.

    The ERA5 grid origin follows EEH2_DOMAIN, so it has to be read from the
    file. Deriving indices from the default domain's corner (75, -20) instead
    made every narrowed domain slice out of range: the slices came back empty
    and STIC died on `np.min` of a zero-size array.
    """
    step = abs(float(axis[1] - axis[0]))
    if float(axis[0]) > float(axis[-1]):        # descending, as ERA5 latitude
        i1 = int(np.floor((float(axis[0]) - hi) / step))
        i2 = int(np.ceil((float(axis[0]) - lo) / step))
    else:
        i1 = int(np.floor((lo - float(axis[0])) / step))
        i2 = int(np.ceil((hi - float(axis[0])) / step))
    i1 = max(i1, 0)
    i2 = min(max(i2, i1), len(axis) - 1)
    return i1, i2


def linearinterp(time, variables):
    y1 = variables[0]
    y2 = variables[1]
    outp = time/60*(y2-y1) + y1
    
    return outp

#Bilinear interpolation function
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

#Calculate the maximum air temperature
def calculatetamax(ta_array):
    ta1 = ta_array[0]
    ta2 = ta_array[1]
    ta3 = ta_array[2]

    return np.maximum(np.maximum(ta1, ta2), ta3)

def Read_ERA5(directory_era5,lat_eco,lon_eco,year,month,day,hour,minute,second):
    #Construct ERA5 filenames accordingly
    
    time = float(minute) + float(second)/60
    time1_str = year + month + day + hour  
    
    time2 = datetime.strptime(time1_str,'%Y%m%d%H') + timedelta(hours=1)
    time2_str = time2.strftime('%Y%m%d%H')

    filename1 = time1_str[0:4] + '_' + time1_str[4:6] + '_' + time1_str[6:8] + '_' + time1_str[8:]
    filename2 = time2_str[0:4] + '_' + time2_str[4:6] + '_' + time2_str[6:8] + '_' + time2_str[8:]

    era5_filename = []
    era5_filename.append(filename1)
    era5_filename.append(filename2)

    #Construct ERA5 filenames for calcualting maximum Ta

    time3_str = year + month + day + '11'
    time4_str = year + month + day + '12'
    time5_str = year + month + day + '13'

    filename3 = time3_str[0:4] + '_' + time3_str[4:6] + '_' + time3_str[6:8] + '_' + time3_str[8:]
    filename4 = time4_str[0:4] + '_' + time4_str[4:6] + '_' + time4_str[6:8] + '_' + time4_str[8:]
    filename5 = time5_str[0:4] + '_' + time5_str[4:6] + '_' + time5_str[6:8] + '_' + time5_str[8:]

    era5_ta_filename = []
    era5_ta_filename.append(filename3)
    era5_ta_filename.append(filename4)
    era5_ta_filename.append(filename5)

    #Crop the ERA5 data to the spatial range of the whole ECOSTRESS image.
    #The grid extent follows EEH2_DOMAIN, so it is read from the file below
    #rather than assumed to be the default lat 75..-35, lon -20..60.
    max_lat = np.amax(lat_eco)
    min_lat = np.amin(lat_eco)
    max_lon = np.amax(lon_eco)
    min_lon = np.amin(lon_eco)

    row1 = row2 = col1 = col2 = None

    sp_array = []
    t2m_array = []
    rh_array = []
    uv_array = []
    ssrd_array = []
    fdir_array = []
    for i in range(2):
        base_filename = 'era5_single_levels-' + era5_filename[i] + ':00.grib'
        filename = os.path.join(directory_era5,base_filename)
        
        #Surface pressure (Pa)
        ds_sp = cfgrib.open_dataset(filename,engine='cfgrib',backend_kwargs={'filter_by_keys': {'shortName': 'sp'}, 'indexpath': ''})
        if row1 is None:
            lat_axis = ds_sp.latitude.data
            lon_axis = ds_sp.longitude.data
            row1, row2 = _grid_span(lat_axis, min_lat, max_lat)
            col1, col2 = _grid_span(lon_axis, min_lon, max_lon)
            if (min_lat < np.amin(lat_axis) or max_lat > np.amax(lat_axis)
                    or min_lon < np.amin(lon_axis) or max_lon > np.amax(lon_axis)):
                # Not fatal: the interpolation clamps to the grid edge below.
                # Worth saying, because the granule then carries edge values for
                # every pixel outside the grid, and ERA5 filenames record no
                # extent, so a cached narrow-domain file looks identical.
                logger.warning(
                    "granule lat [%.2f, %.2f] lon [%.2f, %.2f] exceeds the ERA5 "
                    "grid lat [%.2f, %.2f] lon [%.2f, %.2f] — outside pixels "
                    "use edge values", min_lat, max_lat, min_lon, max_lon,
                    float(np.amin(lat_axis)), float(np.amax(lat_axis)),
                    float(np.amin(lon_axis)), float(np.amax(lon_axis)))
        sp = ds_sp.sp.data[row1:row2+1,col1:col2+1]
        sp_array.append(sp) 
        
        lat_era5 = ds_sp.latitude.data[row1:row2+1]
        lon_era5 = ds_sp.longitude.data[col1:col2+1]  
        
        #Temperature at 2m (K)
        ds_t2m = cfgrib.open_dataset(filename,engine='cfgrib',backend_kwargs={'filter_by_keys': {'shortName': '2t'}, 'indexpath': ''})
        t2m = ds_t2m.t2m.data[row1:row2+1,col1:col2+1] 
        t2m_array.append(t2m)
        
        #Dewpoint temperature at 2m (K)
        ds_d2m = cfgrib.open_dataset(filename,engine='cfgrib',backend_kwargs={'filter_by_keys': {'shortName': '2d'}, 'indexpath': ''})
        d2m = ds_d2m.d2m.data[row1:row2+1,col1:col2+1]
        
        b = 17.625
        c = 243.04
        t0 = 273.15
        t_c = t2m - t0
        dew_c = d2m - t0      
        rh = np.exp(b*c*(dew_c-t_c)/((c+dew_c)*(c+t_c))) #converting from dewpoint temperature to rh (0-1)     
        rh_array.append(rh)
        ds_t2m.close()
        ds_d2m.close()
        ds_sp.close()
        del rh,sp,ds_sp,d2m,ds_d2m,t2m,ds_t2m
             
        #U at 10m (m s-1)
        ds_u10m = cfgrib.open_dataset(filename,engine='cfgrib',backend_kwargs={'filter_by_keys': {'shortName': '10u'}, 'indexpath': ''})
        u10m = ds_u10m.u10.data[row1:row2+1,col1:col2+1]
               
        #V at 10m (m s-1)
        ds_v10m = cfgrib.open_dataset(filename,engine='cfgrib',backend_kwargs={'filter_by_keys': {'shortName': '10v'}, 'indexpath': ''})
        v10m = ds_v10m.v10.data[row1:row2+1,col1:col2+1]
        
        uv = np.sqrt(u10m**2 + v10m**2)    
        uv_array.append(uv)
        ds_u10m.close()
        ds_v10m.close()
        del uv,u10m,ds_u10m,v10m,ds_v10m
        
        #Total solar radiation at surface (J m-2)
        ds_ssrd = cfgrib.open_dataset(filename,engine='cfgrib',backend_kwargs={'filter_by_keys': {'shortName': 'ssrd'}, 'indexpath': ''})
        ssrd = ds_ssrd.ssrd.data[row1:row2+1,col1:col2+1]/3600 #converting from J m-2 to W m-2
        ssrd_array.append(ssrd)
        
        ds_ssrd.close()
        del ssrd,ds_ssrd
        
        #Direct solar radiation at surface (J m-2)
        ds_fdir = cfgrib.open_dataset(filename,engine='cfgrib',backend_kwargs={'filter_by_keys': {'shortName': 'fdir'}, 'indexpath': ''})
        fdir = ds_fdir.fdir.data[row1:row2+1,col1:col2+1]/3600 #converting from J m-2 to W m-2
        fdir_array.append(fdir)
        
        ds_fdir.close()
        del fdir,ds_fdir

    ta_array = []
    for i in range(3):
        base_filename = 'era5_single_levels-' + era5_ta_filename[i] + ':00.grib'
        filename = os.path.join(directory_era5,base_filename)

        #Temperature at 2m (K)
        ds_t2m = cfgrib.open_dataset(filename,engine='cfgrib',backend_kwargs={'filter_by_keys': {'shortName': '2t'}, 'indexpath': ''})
        t2m = ds_t2m.t2m.data[row1:row2+1,col1:col2+1] 
        ta_array.append(t2m)

        ds_t2m.close()
        del ds_t2m,t2m
    
    #Interpolate variables temporally and spatially
    #tempinterp_sp = linearinterp(time,sp_array)
    #spatiotempinterp_sp =  bilinearinterp((lon_eco - lon_era5[0])/0.25,(lat_era5[0] - lat_eco)/0.25,tempinterp_sp,2)
    
    lon_eco = np.clip(lon_eco, np.min(lon_era5)+1e-6, np.max(lon_era5)-1e-6)
    lat_eco = np.clip(lat_eco, np.min(lat_era5)+1e-6, np.max(lat_era5)-1e-6)

    x = (lon_eco - lon_era5[0]) / 0.25
    y = (lat_era5[0] - lat_eco) / 0.25
    
    tempinterp_t2m = linearinterp(time,t2m_array)
    spatiotempinterp_t2m = bilinearinterp(x,y,tempinterp_t2m,2)
    
    tempinterp_rh = linearinterp(time,rh_array)
    spatiotempinterp_rh = bilinearinterp(x,y,tempinterp_rh,2)
    
    tempinterp_uv = linearinterp(time,uv_array)
    spatiotempinterp_uv = bilinearinterp(x,y,tempinterp_uv,2)
    
    tempinterp_ssrd = linearinterp(time,ssrd_array)
    spatiotempinterp_ssrd = bilinearinterp(x,y,tempinterp_ssrd,2)
    
    tempinterp_fdir = linearinterp(time,fdir_array)
    spatiotempinterp_fdir = bilinearinterp(x,y,tempinterp_fdir,2)
    
    spatiotempinterp_fdif = spatiotempinterp_ssrd - spatiotempinterp_fdir

    #Calculate the maximum air temperature
    ta_max = calculatetamax(ta_array)
    spatiotempinterp_tamax = bilinearinterp(x,y,ta_max,2)
    
    return (spatiotempinterp_t2m,spatiotempinterp_rh,spatiotempinterp_uv,spatiotempinterp_fdir,spatiotempinterp_fdif,spatiotempinterp_tamax)

  
