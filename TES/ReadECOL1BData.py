#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Script to run the Temperature Emissivity Separation (TES) algorithm for LST estimation

Created on April 15 2021
@author: Tian Hu at LIST

© 2022 – Luxembourg Institute of Science and Technology
Authors : Tian Hu, Kaniska Mallick
Code licensed under MIT
SPDX-License-Identifier: MIT
"""

# Read the ECOSTRESS L1B_GEO and L1B_RAD data
import h5py
import numpy as np
import glob
import os
from datetime import datetime, timedelta

def timetransform(time):
    t = datetime(2000, 1, 1, 12, 0) + timedelta(seconds=time)
    
    year = t.strftime('%Y')
    month = t.strftime('%m')
    day = t.strftime('%d')
    hour = t.strftime('%H')
    minute = t.strftime('%M')
    second = t.strftime('%S')   
#     time = float(minute) + float(second)/60

    return (year,month,day,hour,minute,second)

def Read_L1B_Data(filename_rad,directory_geo,map_geo):
    #Read radiances
    f_rad = h5py.File(filename_rad, 'r')
    
    r2 = np.array(f_rad['Radiance']['radiance_2']) #Unit: W.m-2.sr-1.um
    r2[r2 < 0] = np.nan
    rqa2 = np.array(f_rad['Radiance']['data_quality_2'])

    r4 = np.array(f_rad['Radiance']['radiance_4'])
    r4[r4 < 0] = np.nan
    rqa4 = np.array(f_rad['Radiance']['data_quality_4'])

    r5 = np.array(f_rad['Radiance']['radiance_5'])
    r5[r5 < 0] = np.nan
    rqa5 = np.array(f_rad['Radiance']['data_quality_5'])
    
    f_rad.close()
    del f_rad

    #Substract the temporal information
    date_str = filename_rad[-26:-18]
    hour_str = filename_rad[-17:-15]
    min_str = filename_rad[-15:-13]
    sec_str = filename_rad[-13:-11]
    
    #Substract the orbit and track information
    orbit_str = filename_rad[-36:-27]
    
    #Compose the filename of geolocation file
    filename_geo = 'ECOSTRESS_L1B_GEO_' + orbit_str + '_' + date_str + 'T' + hour_str + min_str + sec_str
    geo_match = map_geo[filename_geo]
    filename = os.path.join(directory_geo,geo_match)
    f_geo = h5py.File(filename, 'r')
    
    lat_eco = np.array(f_geo['Geolocation']['latitude']) #[-90,90]
    lon_eco = np.array(f_geo['Geolocation']['longitude']) #[-180,180]
    alt_eco = np.array(f_geo['Geolocation']['height']) #Unit: m
    
    lf = np.array(f_geo['Geolocation']['land_fraction'])/100 #[0,1]
    watermask = np.ones(lf.shape,dtype=np.int32)
    watermask[lf > 0.5] = 0 #land

    vza = np.abs(np.array(f_geo['Geolocation']['view_zenith'])) #[0,90]
    sza = np.abs(np.array(f_geo['Geolocation']['solar_zenith'])) #[0,90]
    linetime = np.array(f_geo['Geolocation']['line_start_time_j2000'])
    #geomatching = np.array(f_geo['L1GEOMetadata']['OrbitCorrectionPerformed'])
    
#     timetransform_v = np.vectorize(timetransform)
#     linetime1 = timetransform_v(linetime)
#     time = np.ones(lat_eco.shape)*linetime1.reshape(linetime1.shape[0],1)

    linetime_mean = np.mean(linetime)
    (year,month,day,hour,minute,second) = timetransform(linetime_mean)
    
    f_geo.close()
    del f_geo
 
    return (lat_eco,lon_eco,alt_eco,lf,watermask,vza,sza,year,month,day,hour,minute,second,
            r2,rqa2,r4,rqa4,r5,rqa5,date_str,hour_str,min_str,sec_str,orbit_str)
