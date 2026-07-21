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

# Functions for read the ECOSTRESS L2_LSTE, L1B_GEO and L2_CLOUD data 
import h5py
import numpy as np
import os
from datetime import datetime,timedelta
import wasdi

def timetransform(time):
    t = datetime(2000, 1, 1, 12, 0) + timedelta(seconds=time)
    
    year = t.strftime('%Y')
    month = t.strftime('%m')
    day = t.strftime('%d')
    hour = t.strftime('%H')
    minute = t.strftime('%M')
    second = t.strftime('%S')  
    
    return (year,month,day,hour,minute,second)

def substring(string,n):  
    return string[n]

def Read_ECOSTRESS(filename_lste,directory_geo,directory_cld,map_geo_files,map_cloud_files):
    ########Read LST and LSE data###########
    f_lste = h5py.File(filename_lste, 'r')
    wasdi.wasdiLog('Opening ECOSTRESS LSTE data finished')
    
    lst = np.array(f_lste['LST']).astype(float)
    #lst = np.array(f_lste['SDS']['LST']).astype(float)
    lst = lst*0.02
    lst[(lst < 150) | (lst > 400)] = np.nan
    
    lse = np.array(f_lste['BBE']).astype(float)
    #lse = np.array(f_lste['SDS']['EmisWB']).astype(float)
    lse = lse*0.002 + 0.489999 
    lse[(lse < 0.491999) | (lse > 0.999999)] = np.nan
    
    f_lste.close()
    del f_lste

    #Substract the temporal information
    date_str = filename_lste[-26:-18]
    hour_str = filename_lste[-17:-15]
    min_str = filename_lste[-15:-13]
    sec_str = filename_lste[-13:-11]
    
    #Substract the orbit and track information
    orbit_str = filename_lste[-36:-27]
    
    wasdi.wasdiLog('Reading ECOSTRESS LSTE data finished!')
    
    #########Read geolocation############
    #Compose the filename of geolocation file
    #ECOv002_L1B_GEO_00424_001_20180803T011700_0712_04.h5
    filename_geo = 'ECOv002_L1B_GEO_' + orbit_str + '_' + date_str + 'T' + hour_str + min_str + sec_str
    geo_match = map_geo_files[filename_geo]
    filename = os.path.join(directory_geo,geo_match)
    print(f"Trying to read GEO file:{filename}")
    f_geo = h5py.File(filename, 'r')
    
    wasdi.wasdiLog('Opening ECOSTRESS Geolocation data finished')
    lat_eco = np.array(f_geo['Geolocation']['latitude']) #[-90,90]
    lon_eco = np.array(f_geo['Geolocation']['longitude']) #[-180,180]
    alt_eco = np.array(f_geo['Geolocation']['height']) #Unit: m

    lf = np.array(f_geo['Geolocation']['land_fraction'])/100 #[0,1]
    watermask = np.ones(lf.shape,dtype=np.int32)
    watermask[lf > 0.5] = 0 #land
    
    linetime = np.array(f_geo['Geolocation']['line_start_time_j2000'])
    linetime_mean = np.mean(linetime)
    (year,month,day,hour,minute,second) = timetransform(linetime_mean)
    
    prj_wkt = f_geo['Geolocation'].attrs['Projection_WKT']
    
    f_geo.close()
    del f_geo,lf
    
    wasdi.wasdiLog('Reading ECOSTRESS geolocation data finished!')
    
    #########Read cloud mask#########
    #Compose the filename of cloudmask file
    #ECOv002_L2_CLOUD_00424_001_20180803T011700_0712_04.h5
    filename_cld = 'ECOv002_L2_CLOUD_' + orbit_str + '_' + date_str + 'T' + hour_str + min_str + sec_str
    cloud_match = map_cloud_files[filename_cld]
    filename = os.path.join(directory_cld,cloud_match)
    f_cld = h5py.File(filename, 'r')
    
    wasdi.wasdiLog('Opening ECOSTRESS cloud mask finished')
    
    cld_eco = np.array(f_cld['SDS']['Cloud_final'])
    
    #Convert from int to 8-bit string
    binary_repr_v = np.vectorize(np.binary_repr)
    cld_eco_bin = binary_repr_v(cld_eco,8)
    
    #Extract the bit information
    substring_v = np.vectorize(substring)
    cmask_temp = substring_v(cld_eco_bin,6)
    cmask = np.ones(cmask_temp.shape,dtype=np.int8)
    cmask[cmask_temp == '0'] = 0
    
    f_cld.close()
    del f_cld
    
    wasdi.wasdiLog('Reading ECOSTRESS cloud mask data finished!')
    
    return (lst,lse,lat_eco,lon_eco,alt_eco,watermask,cmask,year,month,day,hour,minute,second,
         date_str,hour_str,min_str,sec_str,orbit_str,prj_wkt)
