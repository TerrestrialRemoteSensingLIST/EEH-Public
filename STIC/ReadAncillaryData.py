# Read the ancillary CGLS NDVI, albedo_dir, albedo_hem and LULC data

import netCDF4 as nc
import numpy as np
import os
import calendar
from pyhdf.SD import SD, SDC
from datetime import datetime
from osgeo import gdal

def warp_raster(x,y,variable):
    #This function aims at finding the matching data for ecostress swath
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
#def warp_raster

def warp_lulc(x,y,lulc):
    #This function aims at finding the matching data for ecostress swath
    xx = np.round(x).astype(int)
    yy = np.round(y).astype(int)
    
    xx = np.clip(xx, 0, lulc.shape[1]-1)
    yy = np.clip(yy, 0, lulc.shape[0]-1)
    
    eco_lulc = lulc[yy,xx]
    
    return eco_lulc

def Read_NDVI(directory_ndvi,latitude,longitude,year_str,month_str,day_str,map_ndvi):
    max_lat = np.amax(latitude)
    min_lat = np.amin(latitude)
    max_lon = np.amax(longitude) 
    min_lon = np.amin(longitude) 
   
    #Calculate the spatial range of the ecostress image in fvc data
    res = 1./336
    row1 = np.floor((80 - max_lat)/res).astype(int)
    row2 = np.ceil((80 - min_lat)/res).astype(int)
    col1 = np.floor((min_lon + 180)/res).astype(int)
    col2 = np.ceil((max_lon + 180)/res).astype(int)
    
    day = int(day_str)
    if day <= 10:
        sday = '01'
    elif day <= 20:
        sday = '11'
    else:
        sday = '21'
        
    key_ndvi = year_str + month_str + sday 
    filename_ndvi = map_ndvi[key_ndvi]
    filename = os.path.join(directory_ndvi,filename_ndvi)

    f_ndvi = nc.Dataset(filename)
    #Extract the corresponding fvc 
    var_ndvi = f_ndvi.variables['NDVI']
    
    if float(year_str) < 2021:
        ndvi = np.array(var_ndvi[row1:row2+1,col1:col2+1]) #the data is scaled automatically in nc data
    else:
        ndvi = np.array(var_ndvi[0,row1:row2+1,col1:col2+1]) #the data is scaled automatically in nc data
    
    ndvi[(ndvi < -0.08) | (ndvi > 0.92)] = np.nan
    
    var_lat = f_ndvi.variables['lat']
    lat_ndvi = np.array(var_lat[row1:row2+1])
    
    var_lon = f_ndvi.variables['lon']
    lon_ndvi = np.array(var_lon[col1:col2+1])
    
    x = (longitude - np.amin(lon_ndvi))/res
    y = (np.amax(lat_ndvi) - latitude)/res
    
    eco_ndvi = warp_raster(x,y,ndvi)
    
#     inds = np.isnan(eco_ndvi)
#     eco_ndvi[inds] = -9999
    
    f_ndvi.close()
    return eco_ndvi
    
def Read_LAI(directory_lai,latitude,longitude,year_str,month_str,day_str,map_lai):
    max_lat = np.amax(latitude)
    min_lat = np.amin(latitude)
    max_lon = np.amax(longitude) 
    min_lon = np.amin(longitude) 
   
    #Calculate the spatial range of the ecostress image in fvc data
    res = 1./336
    row1 = np.floor((80 - max_lat)/res).astype(int)
    row2 = np.ceil((80 - min_lat)/res).astype(int)
    col1 = np.floor((min_lon + 180)/res).astype(int)
    col2 = np.ceil((max_lon + 180)/res).astype(int)
    
    day = int(day_str)
    if day <= 10:
        sday = '10'
    elif day <= 20:
        sday = '20'
    else:
        lastday = calendar.monthrange(int(year_str),int(month_str))[1]
        sday = str(lastday)

    key_lai = year_str + month_str + sday    
    filename_lai = map_lai[key_lai]
    filename = os.path.join(directory_lai,filename_lai)

    f_lai = nc.Dataset(filename)   
    #Extract the corresponding fvc 
    var_lai = f_lai.variables['LAI']
    
    #Calculate the accumulated year
    accum_year = float(year_str) + float(month_str)/12.0
    if accum_year < 2020.7:
        lai = np.array(var_lai[row1:row2+1,col1:col2+1]) #the data is scaled automatically in nc data
    else:
        lai = np.array(var_lai[0,row1:row2+1,col1:col2+1]) #the data is scaled automatically in nc data
    
    lai[(lai < 0) | (lai > 7)] = np.nan
    
    var_lat = f_lai.variables['lat']
    lat_lai = np.array(var_lat[row1:row2+1])
    
    var_lon = f_lai.variables['lon']
    lon_lai = np.array(var_lon[col1:col2+1])
    
    x = (longitude - np.amin(lon_lai))/res
    y = (np.amax(lat_lai) - latitude)/res
    
    eco_lai = warp_raster(x,y,lai)

#     inds = np.isnan(eco_lai)
#     eco_lai[inds] = -9999
    
    f_lai.close()
    return eco_lai
    
def Read_FVC(directory_fvc,latitude,longitude,year_str,month_str,day_str,map_fvc):
    max_lat = np.amax(latitude)
    min_lat = np.amin(latitude)
    max_lon = np.amax(longitude) 
    min_lon = np.amin(longitude) 
   
    #Calculate the spatial range of the ecostress image in fvc data
    res = 1./336
    row1 = np.floor((80 - max_lat)/res).astype(int)
    row2 = np.ceil((80 - min_lat)/res).astype(int)
    col1 = np.floor((min_lon + 180)/res).astype(int)
    col2 = np.ceil((max_lon + 180)/res).astype(int)
    
    day = int(day_str)
    if day <= 10:
        sday = '10'
    elif day <= 20:
        sday = '20'
    else:
        lastday = calendar.monthrange(int(year_str),int(month_str))[1]
        sday = str(lastday)

    key_fvc = year_str + month_str + sday     
    filename_fvc = map_fvc[key_fvc]
    filename = os.path.join(directory_fvc,filename_fvc)
    print('Opening FVC file='+filename)

    f_fvc = nc.Dataset(filename)   
    #Extract the corresponding fvc 
    var_fvc = f_fvc.variables['FCOVER']
    
    #Calculate the accumulated year
    accum_year = float(year_str) + float(month_str)/12.0
    if accum_year < 2020.7:
        fvc = np.array(var_fvc[row1:row2+1,col1:col2+1]) #the data is scaled automatically in nc data
    else:
        fvc = np.array(var_fvc[0,row1:row2+1,col1:col2+1]) #the data is scaled automatically in nc data
    
    fvc[(fvc < 0) | (fvc > 1)] = np.nan
    fvc[fvc == 1] = 0.99
    
    var_lat = f_fvc.variables['lat']
    lat_fvc = np.array(var_lat[row1:row2+1])
    
    var_lon = f_fvc.variables['lon']
    lon_fvc = np.array(var_lon[col1:col2+1])
    
    x = (longitude - np.amin(lon_fvc))/res
    y = (np.amax(lat_fvc) - latitude)/res
    
    eco_fvc = warp_raster(x,y,fvc)

#     inds = np.isnan(eco_fvc)
#     eco_fvc[inds] = -9999
    
    f_fvc.close()
    return eco_fvc
    
def Read_ALB_DIR(directory_alb_dir,latitude,longitude,year_str,month_str,day_str,map_albdir1,map_albdir2):
    max_lat = np.amax(latitude)
    min_lat = np.amin(latitude)
    max_lon = np.amax(longitude) 
    min_lon = np.amin(longitude)
    
    time = datetime(int(year_str),int(month_str),int(day_str))
    doy = (time - datetime(int(year_str), 1, 1)).days + 1
    doy_str = '{0:03d}'.format(doy)
    #Calculate the accumulated year
    accum_year = float(year_str) + float(month_str)/12.0
    #if accum_year < 2020.5:
    if False:
        #Calculate the spatial range of the ecostress image in fvc data
        res = 1./112.0
        row1 = int(np.floor((80 - max_lat)/res))
        row2 = int(np.ceil((80 - min_lat)/res))
        col1 = int(np.floor((min_lon + 180)/res))
        col2 = int(np.ceil((max_lon + 180)/res))

        day = int(day_str)
        if day <= 10:
            sday = '0'
        elif day <= 20:
            sday = '1'
        else:
            sday = '2'

        key_albdir = 'c_gls_ALDH_' + year_str + month_str + sday
        filename_alb_dir = map_albdir1[key_albdir]
        filename = os.path.join(directory_alb_dir,filename_alb_dir)

        f_alb = nc.Dataset(filename)   
        #Extract the corresponding albedo 
        var_alb = f_alb.variables['AL_DH_BB']

        alb = np.array(var_alb[0,row1:row2+1,col1:col2+1]) #the data is scaled automatically in nc data
        alb[(alb < 0) | (alb > 1)] = np.nan

        var_lat = f_alb.variables['lat']
        lat_alb = np.array(var_lat[row1:row2+1])

        var_lon = f_alb.variables['lon']
        lon_alb = np.array(var_lon[col1:col2+1])

        x = (longitude - np.amin(lon_alb))/res
        y = (np.amax(lat_alb) - latitude)/res

        eco_alb = warp_raster(x,y,alb)
        
#         inds = np.isnan(eco_alb)
#         eco_alb[inds] = -9999
        
        f_alb.close()       
        return eco_alb
    else:
        #Calculate the spatial range of the ecostress image in fvc data
        res = 0.05
        row1 = int(np.floor((89.975 - max_lat)/res))
        row2 = int(np.ceil((89.975 - min_lat)/res))
        col1 = int(np.floor((min_lon + 179.975)/res))
        col2 = int(np.ceil((max_lon + 179.975)/res))

        key_albdir = year_str + doy_str
        filename_alb_dir = map_albdir2[key_albdir]
        filename = os.path.join(directory_alb_dir,filename_alb_dir)
    
        print('Reading file '+filename)
        f_alb = SD(filename, SDC.READ)
        #Extract the corresponding albedo
        var_alb = f_alb.select('Albedo_BSA_shortwave')
        
        alb = var_alb[row1:row2+1,col1:col2+1].astype(float)*0.001 
        alb[(alb < 0) | (alb > 1)] = np.nan
        
        var_lat = np.linspace(89.975,-89.975,3600)
        lat_alb = var_lat[row1:row2+1]
        
        var_lon = np.linspace(-179.975,179.975,7200)
        lon_alb = var_lon[col1:col2+1]
        
        x = (longitude - np.amin(lon_alb))/res
        y = (np.amax(lat_alb) - latitude)/res
        
        eco_alb = warp_raster(x,y,alb)
        
        f_alb.end()
        return eco_alb 
    
def Read_ALB_HEM(directory_alb_hem,latitude,longitude,year_str,month_str,day_str,map_albhem1,map_albhem2):
    max_lat = np.amax(latitude)
    min_lat = np.amin(latitude)
    max_lon = np.amax(longitude) 
    min_lon = np.amin(longitude)
    
    time = datetime(int(year_str),int(month_str),int(day_str))
    doy = (time - datetime(int(year_str), 1, 1)).days + 1
    doy_str = '{0:03d}'.format(doy)
    
    #Calculate the accumulated year
    accum_year = float(year_str) + float(month_str)/12.0

    if False: #accum_year < 2020.5:
        #Calculate the spatial range of the ecostress image in fvc data
        res = 1./112.0
        row1 = np.floor((80 - max_lat)/res).astype(int)
        row2 = np.ceil((80 - min_lat)/res).astype(int)
        col1 = np.floor((min_lon + 180)/res).astype(int)
        col2 = np.ceil((max_lon + 180)/res).astype(int)

        day = int(day_str)
        if day <= 10:
            sday = '0'
        elif day <= 20:
            sday = '1'
        else:
            sday = '2'

        key_albhem = 'c_gls_ALBH_' + year_str + month_str + sday
        filename_alb_hem = map_albhem1[key_albhem]
        filename = os.path.join(directory_alb_hem,filename_alb_hem)

        f_alb = nc.Dataset(filename)   
        #Extract the corresponding fvc 
        var_alb = f_alb.variables['AL_BH_BB']

        alb = np.array(var_alb[0,row1:row2+1,col1:col2+1]) #the data is scaled automatically in nc data
        alb[(alb < 0) | (alb > 1)] = np.nan

        var_lat = f_alb.variables['lat']
        lat_alb = np.array(var_lat[row1:row2+1])

        var_lon = f_alb.variables['lon']
        lon_alb = np.array(var_lon[col1:col2+1])

        x = (longitude - np.amin(lon_alb))/res
        y = (np.amax(lat_alb) - latitude)/res

        eco_alb = warp_raster(x,y,alb)

    #     inds = np.isnan(eco_alb)
    #     eco_alb[inds] = -9999

        f_alb.close()
        return eco_alb
    else:
        #Calculate the spatial range of the ecostress image in fvc data
        res = 0.05
        row1 = int(np.floor((89.975 - max_lat)/res))
        row2 = int(np.ceil((89.975 - min_lat)/res))
        col1 = int(np.floor((min_lon + 179.975)/res))
        col2 = int(np.ceil((max_lon + 179.975)/res))

        key_albhem = year_str + doy_str
        filename_alb_hem = map_albhem2[key_albhem]
        filename = os.path.join(directory_alb_hem,filename_alb_hem)
    
        f_alb = SD(filename, SDC.READ)
        #Extract the corresponding albedo
        var_alb = f_alb.select('Albedo_WSA_shortwave')
        
        alb = var_alb[row1:row2+1,col1:col2+1].astype(float)*0.001 
        alb[(alb < 0) | (alb > 1)] = np.nan
        
        var_lat = np.linspace(89.975,-89.975,3600)
        lat_alb = var_lat[row1:row2+1]
        
        var_lon = np.linspace(-179.975,179.975,7200)
        lon_alb = var_lon[col1:col2+1]
        
        x = (longitude - np.amin(lon_alb))/res
        y = (np.amax(lat_alb) - latitude)/res
        
        eco_alb = warp_raster(x,y,alb)
        
        f_alb.end()
        return eco_alb
    
def Read_LULC(directory_lulc,latitude,longitude,year_str):  
    max_lat = np.amax(latitude)
    min_lat = np.amin(latitude)
    max_lon = np.amax(longitude) 
    min_lon = np.amin(longitude) 
   
    #Calculate the spatial range of the ecostress image in fvc data
    res = 1./1008
    lat_up = 80 - 0.5*res
    lon_left = -180 + 0.5*res
    row1 = np.floor((lat_up - max_lat)/res)
    row2 = np.ceil((lat_up - min_lat)/res)
    col1 = np.floor((min_lon - lon_left)/res)
    col2 = np.ceil((max_lon - lon_left)/res)
    
    filename_lulc = 'PROBAV_LC100_global_v3.0.1_2018-conso_Discrete-Classification-map_EPSG-4326.tif'
    filename = os.path.join(directory_lulc,filename_lulc)
    dataset = gdal.Open(filename,gdal.GA_ReadOnly)
    band = dataset.GetRasterBand(1)
    xsize = int(col2-col1+1)
    ysize = int(row2-row1+1)
    lulc = band.ReadAsArray(int(col1),int(row1),xsize,ysize)
    
    lat_lulc_max = lat_up - row1*res
    lon_lulc_min = lon_left + col1*res
    
    x = (longitude - lon_lulc_min)/res
    y = (lat_lulc_max - latitude)/res
    
    eco_lulc = warp_lulc(x,y,lulc)
    
    band = None
    dataset = None
    
    return eco_lulc
    
def Read_Ancillary(directory_fvc,directory_alb_dir,directory_alb_hem,
                   directory_lulc,latitude,longitude,year_str,month_str,day_str,map_fvc,map_albdir1,map_albdir2,map_albhem1,map_albhem2):   
    
    #Read FVC data
    fvc = Read_FVC(directory_fvc,latitude,longitude,year_str,month_str,day_str,map_fvc)
    
    print('Reading FVC data finished!')
    
    #Read directional-hemispherical albedo data
    alb_dir = Read_ALB_DIR(directory_alb_dir,latitude,longitude,year_str,month_str,day_str,map_albdir1,map_albdir2)
    
    print('Reading ALB_DIR data finished!')
    
    #Read bi-hemispherical albedo data
    alb_hem = Read_ALB_HEM(directory_alb_hem,latitude,longitude,year_str,month_str,day_str,map_albhem1,map_albhem2)
    
    print('Reading ALB_HEM data finished!')
    
    #Read LULC data 
    lulc = Read_LULC(directory_lulc,latitude,longitude,year_str)
    
    print('Reading LULC data finished!')
 
    return (fvc,alb_dir,alb_hem,lulc)
