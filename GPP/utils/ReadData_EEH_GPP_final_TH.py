#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Script to run the hybrid model for GPP and WUE estimation

Created on May 1 2025
@author: Ziyu Lin at LIST

© 2025 – Luxembourg Institute of Science and Technology
Authors : Ziyu Lin, Kaniska Mallick, Tian Hu (tian.hu@list.lu)
Code licensed under MIT
SPDX-License-Identifier: MIT
"""

import sys
sys.path.insert(0,'./utils')

import os
import numpy as np
import h5py
import glob
import calendar
import pandas as pd
from datetime import datetime,timedelta
import xarray as xr
import rioxarray
from rioxarray.merge import merge_arrays
from rasterio.enums import Resampling
import pandas as pd
import re
import time
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from shapely.geometry import Point, Polygon, box
import xarray as xr
import geopandas as gpd
import re
import sys
from pyhdf.SD import SD, SDC
from scipy.interpolate import griddata

import wasdi
 
def load_S3_paths(path_txt): ##'S3_file_paths_V3.txt'
    # Initialize dictionary to hold categorized file paths
    categories = [
        'geo', 'lste', 'cld', 'fvc', 'fapar', 'lai', 'albdir', 'albhem',
        'era5', 'oco2', 'mota', 'CI', 'GLC30', 'LUCC', 'PAR', 'PARH', 'STIC'
    ]
    file_dict = {cat: [] for cat in categories}

    # Read and categorize paths
    with open(path_txt, 'r') as f:
        for line in f:
            if ': ' in line:
                category, path = line.strip().split(': ', 1)
                if category in file_dict:
                    file_dict[category].append(path)

    return file_dict


def bluesky_albedo_calculator(path_MOTA, eco_bound):
    """
    Calculate blue sky albedo from black and white sky albedo using Long & Ackerman (2000).
    
    Args:
        path_MOTA (str): Path to MODIS MCD43C3 HDF file.
        eco_bound:  (lat_max, lat_min, lon_max, lon_min)

    Returns:
        xarray.DataArray: Blue sky albedo with quality mask and coordinates

    """
    lat_max, lat_min, lon_max, lon_min = eco_bound
    # Load HDF4 datasets
    hdf = SD(path_MOTA, SDC.READ)
    BSA = hdf.select('Albedo_BSA_shortwave').get()
    WSA = hdf.select('Albedo_WSA_shortwave').get()
    SZA = hdf.select('Local_Solar_Noon').get()
    mask = hdf.select('Albedo_Quality').get()

    # Create lat/lon grids
    rows, cols = BSA.shape
    lat = np.linspace(90, -90, rows)
    lon = np.linspace(-180, 180, cols)
    

    # Stack into DataArray with band dimension
    data_stack = np.stack([BSA, WSA, SZA, mask], axis=0)
    da_stack = xr.DataArray(
        data_stack,
        dims=['band', 'y', 'x'],
        coords={
            'band': ['BSA', 'WSA', 'SZA','QA'],
            'y': lat,
            'x': lon
        }
    ).rio.write_crs("EPSG:4326")

    # Clip to bounding box
    da_clipped = da_stack.rio.clip_box(
        minx=lon_min, miny=lat_min, maxx=lon_max, maxy=lat_max,
        auto_expand=False
    )

    # Extract clipped bands
    BSA_clip = da_clipped.sel(band='BSA').values
    WSA_clip = da_clipped.sel(band='WSA').values
    SZA_clip = da_clipped.sel(band='SZA').values
    QA_clip = da_clipped.sel(band='QA').values
 

    # Compute blue sky albedo
    SZA_rad = np.radians(SZA_clip)
    cos_sza = np.clip(np.cos(SZA_rad), 1e-6, None)# mask non-positive values
    f_dif = (cos_sza ** -0.8425) * 0.123
    scale_factor = 1000
    bluesky = (BSA_clip * (1 - f_dif) + WSA_clip * f_dif) / scale_factor
    bluesky[QA_clip != 0] = np.nan

    # Wrap result in DataArray
    lon_1d, lat_1d = da_clipped['x'].values, da_clipped['y'].values
    lon_grid, lat_grid = np.meshgrid(lon_1d, lat_1d)
    da_bluesky = xr.DataArray(
        bluesky,
         coords={
            'y': (['y', 'x'], lat_grid),
            'x': (['y', 'x'], lon_grid)
        },
        dims=['y', 'x']
    )


    return da_bluesky



def wrap_to_ECOSTRESS(dataset, lat_eco, lon_eco, method='linear'):
    """
    Regrid MODIS data to ECOSTRESS resolution using 2D bilinear interpolation.

    Args:
        dataset (xarray.DataArray): Source DataArray with 2D lat/lon coordinates.
        lat_eco, lon_eco (np.ndarray): Target lat/lon grids (2D).
        method (str): Interpolation method ('nearest', 'linear', 'cubic').

    Returns:
        xarray.DataArray: Interpolated DataArray on ECOSTRESS grid.
    """
    # Extract source lat/lon and data
    lat_src_2d, lon_src_2d = get_lat_lon(dataset)
    data_src = dataset.values
    
    # Flatten source grid
    points = np.column_stack((lat_src_2d.ravel(), lon_src_2d.ravel()))
    values = data_src.ravel()

    # Flatten target grid
    target_points = np.column_stack((lat_eco.ravel(), lon_eco.ravel()))

    # Interpolate using bilinear method
    interpolated = griddata(points, values, target_points, method='linear')

    # Reshape to target grid shape
    interpolated_2d = interpolated.reshape(lat_eco.shape)
  

    # Return as xarray.DataArray
    da_interp = xr.DataArray(
        interpolated_2d,
        dims=('y', 'x'),
        coords={'lat': (('y', 'x'), lat_eco), 'lon': (('y', 'x'), lon_eco)},
        name='interpolated'
    )

    return da_interp


def filter_tile_bounds(path, eco_bound):
    """
    Extract bounding box from filename, and check if tile overlaps with eco_bound.
    eco_bound = (lat_max, lat_min, lon_max, lon_min)
    """
    match = re.search(r'([EW])(\d+)([NS])(\d+)', path)
    if not match:
        return None
    
    # coordinate at upper left corner
    lon_dir = match.group(1)
    lon = int(match.group(2))
    lat_dir = match.group(3)
    lat = int(match.group(4))
    tile_size=5 #5*5 degree per tile

    # Adjust latitude sign
    lat = lat if lat_dir == 'N' else -lat
    lon = lon if lon_dir == 'E' else -lon

    # Define tile bounding box 
    t_min_lon = lon
    t_max_lon = lon + tile_size
    t_max_lat = lat
    t_min_lat = lat - tile_size 
    
    e_max_lat, e_min_lat, e_max_lon, e_min_lon = eco_bound
    
  
    # Check for intersection
    intersects = not (
        t_max_lon < e_min_lon or
        t_min_lon > e_max_lon or
        t_max_lat < e_min_lat or
        t_min_lat > e_max_lat
    )
    
    return intersects


def read_and_wrap_GLC30(list_tif, lat_eco, lon_eco, eco_bound, year):
    """
    Reads GLC30 raster files, clips to bounding box, aligns to ECOSTRESS grid,
    and updates ECOSTRESS lat/lon grid with GLC30 values.

    Parameters:
    - list_tif: List of file paths to georeferenced raster images
    - lat_eco, lon_eco (np.ndarray): Target lat/lon grids (2D)
    - eco_bound: Tuple of (lat_max, lat_min, lon_max, lon_min)
    - year: Year of interest (2016–2022). If outside range, defaults to last band

    Returns:
    - xarray.DataArray with GLC30 values aligned to ECOSTRESS grid
    """
    lat_max, lat_min, lon_max, lon_min = eco_bound
    #band_ind = year - 2016 if year <= 2022 else 7
    band_ind = year - 2016 if year <= 2022 else 6

    glc30_grid = xr.DataArray(
    data=np.ones(lat_eco.shape, dtype=np.int16) * -9999, 
    coords={
        "y": (["y", "x"], lat_eco),
        "x": (["y", "x"], lon_eco)
        },
        dims=["y", "x"],
        name="GLC30"
    )

 
    for path in list_tif:
        try:
            # Load raster
            da = rioxarray.open_rasterio(path, masked=True)

            # Select band
            da_band = da.isel(band=band_ind)
            
            # Create mask based on bounds
#             lat = da_band['y']
#             lon = da_band['x']
#             mask_bound = (lat >= lat_min) & (lat <= lat_max) & (lon >= lon_min) & (lon <= lon_max)     
#             da_clipped = da_band.where(mask_bound, drop=True)  # drop=False keeps shape, but fills outside with NaN
        
            # Clip to bounding box
            da_clipped = da_band.rio.clip_box(
                minx=lon_min, miny=lat_min, maxx=lon_max, maxy=lat_max,
                auto_expand=False
            )
  
            scaled = (da_clipped.data / 10).astype(np.int16)
            nonveg_mask = np.isin(scaled, [14, 15, 19, 20, 21, 22])
            scaled[nonveg_mask] = 0
            # Replace 1 with 2
            scaled[scaled == 1] = 2 #CRO
            da_clipped.data = scaled
            
            # plt.figure()
            # im= plt.imshow(da_clipped, cmap='viridis')#, extent=[lon_min, lon_max, lat_min, lat_max])
            # plt.colorbar(im, label='lucc')
            # plt.tight_layout()
            # plt.show()
        
            # Align to ECOSTRESS grid using nearest-neighbor interpolation
            da_aligned = da_clipped.interp(
                x=glc30_grid['x'],
                y=glc30_grid['y'],
                method="nearest"
            )

            # Update ECOSTRESS grid where nodata
            valid = (glc30_grid >= 0)  # or whatever your valid range is,including -9999 and nan
            glc30_grid = xr.where(valid, glc30_grid, da_aligned)
            

     
        except Exception as e:
            print(f"Error processing {path}: {e}")
            
        # SAV: Grass and Shrub found in Africa between latitude 15° North and 30 degrees South
        scaled = glc30_grid.data
        scaled[(glc30_grid['y'].values<15)&(glc30_grid['y'].values>-30)&(scaled == 13)]=14 #grass SAV 
        scaled[(glc30_grid['y'].values<15)&(glc30_grid['y'].values>-30)&(scaled == 12)]=15 #woody WSA   
        glc30_grid.data = scaled
        
    return glc30_grid


def Create_cloud_mask(path_cld):
    """Create a mask to select non-cloudy pixels based on provided data."""
    
    with h5py.File(path_cld, "r") as f_cld:
        # Extract cloud
        eco_cld =  np.array(f_cld['SDS']['Cloud_final']) 
        # Create cloud mask (binary) 
        # FOR v001 data 0=clear 1= cloud -1 or 255=missing
        cmask = (eco_cld == 0).astype(np.int8)
        
    return cmask


def Extract_ECOSTRESS_from_mask(file_path, mask_cld_ffp, variables=["LST"]):
    """Extract values from ECOSTRESS footprint using provided footprint and cloud masks."""
    
    # Open the HDF5 file using xarray with the appropriate engine
    ds = xr.open_dataset(file_path, engine="h5netcdf", phony_dims="sort")

    # Validate variables
    missing_vars = [var for var in variables if var not in ds.variables]
    if missing_vars:
        raise ValueError(f"Variables not found in dataset: {missing_vars}")

    # Extract values
    # .where(mask_cld_ffp)
    return [ds[var].values.astype(float) for var in variables]

def Extract_STIC(filename_lste, mask_cld_ffp):
    """Extract STIC outputs values of EEH2"""
    variables = ["ETD", "gah", "gsc"]
    if np.sum(mask_cld_ffp) > 0:  # Ensuring valid data within footprint&cloud mask
        ETD, gah, gsc = Extract_ECOSTRESS_from_mask(filename_lste, mask_cld_ffp,variables)
        gt = gsc*gah/(gsc+gah)
        # Apply mask and filter positive values
        ETD_filtered = np.where(mask_cld_ffp&(ETD<20)&(ETD>0), ETD, np.nan) #remove abnormal pixels
        gt_filtered = np.where(mask_cld_ffp&(gah>0)&(gsc>0), gt, np.nan).clip(0.0001,0.2) #range adaped from STIC
        return [ETD_filtered, gt_filtered]
    else:
        nan_array = np.full_like(mask_cld_ffp, np.nan, dtype=float)
        return [nan_array for i in range(len(variables))]

    
def read_geo(path_geo):
    """read geolocation data."""
    
    with h5py.File(path_geo, "r") as f_geo:
        # Extract latitude and longitude information
        eco_lat =  np.array(f_geo['Geolocation']['latitude'])
        eco_lon =  np.array(f_geo['Geolocation']['longitude'])
        sza_lon =  np.array(f_geo['Geolocation']['solar_zenith']) #in degree
   
    return eco_lat,eco_lon,sza_lon

def timetransform(time):
    t = datetime(2000, 1, 1, 12, 0) + timedelta(seconds=time)
    
    year = t.strftime('%Y')
    month = t.strftime('%m')
    day = t.strftime('%d')
    hour = t.strftime('%H')
    minute = t.strftime('%M')
    second = t.strftime('%S')   

    return (year,month,day,hour,minute,second)

def regex_extract_ecostress(filename):
    # Create a regex pattern to extract relevant details
    pattern = r".*_(\d+)_(\d+)_(\d{4})(\d{2})(\d{2})T(\d{2})(\d{2})(\d{2}).*"
    
    # Apply regex search
    match = re.search(pattern, filename)
    print(f"Searching {pattern} in {filename}")

    if match:
        orbit_str = match.group(1)      # Extract  Orbit number
        scene_str = match.group(2)    # Extract Scene ID
        year_str = match.group(3)      # Extract Year
        month_str = match.group(4)     # Extract Month
        day_str = match.group(5)      # Extract Date
        hour_str = match.group(6)    # Extract Timestamp
        min_str = match.group(7)    
        sec_str = match.group(8)    
        return (orbit_str,scene_str,year_str,month_str,day_str,hour_str,min_str,sec_str)
    else:
        print("No match found in file name.")
        return None

def extract_from_global_data(path_globe, eco_bound, variable='LAI'):
    """
    Extract and interpolate geospatial data from a global NetCDF file within a specified bounding box.

    Args:
        path_globe (str): Path to NetCDF file, only supports 2d data without time dim
        eco_bound (tuple): (lat_max, lat_min, lon_max, lon_min)
        hour (int): Overpass hour.
        minute (int): Overpass minute.

    Returns:
       tuple: (DataArray at overpass time, DataArray summed over time)
    """
  
    # read data without preloading
    if path_globe.endswith(('.nc', '.nc4', '.h5', '.hdf5', '.cdf')):
        # Handle NetCDF/HDF file
        ds = xr.open_dataset(path_globe, decode_timedelta=False, engine='netcdf4')
     
    elif path_globe.endswith('.grib'):
        ds = xr.open_dataset(path_globe,decode_timedelta=False,engine='cfgrib',backend_kwargs={ "errors": "ignore" })
        # datasets = [xr.open_dataset(
        #             path_globe,
        #             engine="cfgrib",
        #             backend_kwargs={"filter_by_keys": {"shortName": var}}
        #         ) for var in variable]
        # datasets = [ds for ds in datasets if ds is not None]
        # ds = xr.merge(datasets)
    elif path_globe.endswith('.tif'):
        ds = xr.open_dataset(path_globe, decode_timedelta=False,engine="rasterio")


    else:
        raise ValueError("Only .nc,.nc4,.h5,.hdf5,.cdf,.grib,.tif files are supported.")
     # standarize the dims name
    # print(ds)
    
    
    # Select the variable(s) as a DataArray
    ds_data = ds[variable]
    # Get the actual spatial dimension names, !! important!!
    lat_name = [name for name in ds_data.coords if ('lat' in name)|(name=='y')][0]
    lon_name = [name for name in ds_data.coords if ('lon' in name)|(name=='x')][0]
    dims_2d = (lat_name, lon_name) #2d

    #Set spatial dimensions for rioxarray
    ds_data.rio.set_spatial_dims(x_dim=lon_name, y_dim=lat_name, inplace=True)
    # Ensure CRS is set
    ds_data.rio.write_crs("EPSG:4326", inplace=True)
   
    # print(ds_data)
    
    # Bounding box
    lat_max, lat_min, lon_max, lon_min = eco_bound
    # mask_bound = (lat_2d >= lat_min) & (lat_2d <= lat_max) & \
    #              (lon_2d >= lon_min) & (lon_2d <= lon_max)

    # Read from bound to save time !!!!
    # Clip to bounding box
    da_clipped = ds_data.rio.clip_box(
        minx=lon_min, miny=lat_min, maxx=lon_max, maxy=lat_max,
        auto_expand=True  # ensures it trims to the box
    )
    # print(da_clipped)
    
    #prepare 2-d lat lon for matching with ECOSTRESS
    # Extract 1D coordinates
    lat_1d, lon_1d = da_clipped[lat_name].values,da_clipped[lon_name].values

    # Create 2D coordinate grids to match ecostress
    lon_2d, lat_2d = np.meshgrid(lon_1d, lat_1d)

    # Assign 2D coordinates
    da_clipped = da_clipped.assign_coords({
        lat_name: (dims_2d, lat_2d),
        lon_name: (dims_2d, lon_2d)
    })


    return da_clipped    

def extract_PAR_from_global_data(path_globe, eco_bound, hour, minute):
    """
    Extract and interpolate geospatial data from a global NetCDF file within a specified bounding box.

    Args:
        path_globe (str): Path to NetCDF file.
        eco_bound (tuple): (lat_max, lat_min, lon_max, lon_min)
        hour (int): Overpass hour.
        minute (int): Overpass minute.

    Returns:
       tuple: (DataArray at overpass time, DataArray avearaged over time)
    """
  
    
    if path_globe.endswith('.nc'):
        ds = xr.open_dataset(path_globe, decode_timedelta=False, engine='netcdf4')
        ds_data = ds['PAR']

        # Extract 1D coordinates
        lon_1d = ds.lon.values
        lat_1d = ds.lat.values

        # Create 2D coordinate grids
        lon_2d, lat_2d = np.meshgrid(lon_1d, lat_1d)

        # Assign 2D coordinates
        ds_data = ds_data.assign_coords(lat=(["lat", "lon"], lat_2d),
                                        lon=(["lat", "lon"], lon_2d))

    else:
        raise ValueError("Only .nc files are supported.")

    # Bounding box
    lat_max, lat_min, lon_max, lon_min = eco_bound
    mask_bound = (lat_2d >= lat_min) & (lat_2d <= lat_max) & \
                 (lon_2d >= lon_min) & (lon_2d <= lon_max)

    # Get time dimension size
    time_len = ds_data.sizes['time']

    # Expand mask to 3D: (time, lat, lon)
    mask_3d = np.broadcast_to(mask_bound, (time_len, *mask_bound.shape))

    # Wrap in DataArray with matching dims and coords
    mask_da = xr.DataArray(
        mask_3d,
        dims=ds_data.dims,
        coords={ 
            'time': ds_data['time'],
            'lat': ds_data['lat'],
            'lon': ds_data['lon']
        }
    )

    # Apply mask
    da_clipped = ds_data.where(mask_da, drop=True)
      
    # Estimate time index (assuming 30-min intervals)
    time_index = int(hour * 2 + np.round(minute / 30))

    # Extract time slice
    da_PARH = da_clipped.isel(time=time_index)

    # mean over time
    da_PARmean = da_clipped.mean(dim='time')

    return da_PARH, da_PARmean


def interpolate_ERA5_value(year, month, day, hour, minute, files_era5, eco_bound):

    # Step 1: find path before and after overpass time
    target_time = datetime(year, month, day, hour, minute)

    before = None
    after = None
    min_before_diff = float('inf')
    min_after_diff = float('inf')

    for path in files_era5:
        match = re.search(r'(\d{4})_(\d{2})_(\d{2})_(\d{2}):(\d{2})', path)
        if match:
            file_time = datetime(
                int(match.group(1)), int(match.group(2)),
                int(match.group(3)), int(match.group(4)), int(match.group(5))
            )
            diff = abs((file_time - target_time).total_seconds())

            if file_time < target_time and diff < min_before_diff:
                min_before_diff = diff
                before_path = path
                before_time = file_time
            elif file_time > target_time and diff < min_after_diff:
                min_after_diff = diff
                after_path = path
                after_time = file_time
    if before_path is None or after_path is None:
        raise ValueError("Target date is outside the range of observation ERA5.")

    # Step 2: Read data
    #t2m	2 m temperature	K	Air temperature at 2 meters above surface
    #d2m	2 m dew point temperature	K	Dew point at 2 meters, used to estimate humidity
    ds_before = extract_from_global_data(before_path, eco_bound, variable=['t2m','d2m'])
    ds_after  = extract_from_global_data(after_path, eco_bound, variable=['t2m','d2m'])

    # Step 3: Interpolate
    weight = (target_time - before_time) / (after_time - before_time)
    ds_interp = ds_before + weight * (ds_after - ds_before)

    # Step 4: Fill missing values
    ds_filled = ds_interp.fillna(ds_after).fillna(ds_before)
    
    # Step 5: calcululate vpd
    t2m = ds_filled['t2m'] #°C = K − 273.15
    t2m_C = t2m-273.15#in °C

    d2m = ds_filled['d2m']
    vpd_hPa = calculate_vpd_from_t2m_d2m(t2m, d2m)

    return t2m_C, vpd_hPa



def calculate_vpd_from_t2m_d2m(t2m, d2m):
    """
    Calculate VPD from 2-meter temperature and dew point temperature.

    Parameters:
    t2m : np.array
        Air temperature at 2 meters (Kelvin)
    d2m : np.array
        Dew point temperature at 2 meters (Kelvin)

    Returns:
    vpd : np.array
        Vapor Pressure Deficit in hPa
    """
    # Convert to Celsius
    T_C = t2m - 273.15
    TD_C = d2m - 273.15

    # Saturation vapor pressure
    e_s = 610.94 * np.exp((17.625 * T_C) / (T_C + 243.04))

    # Actual vapor pressure
    e_a = 610.94 * np.exp((17.625 * TD_C) / (TD_C + 243.04))

    # VPD in kPa
    vpd = (e_s - e_a) / 100.0
    return vpd


def get_lat_lon(dataset):
    # Try common latitude keys
    lat_keys = ['y','lat', 'latitude']
    lon_keys = ['x','lon', 'longitude']
    
    coord_names = [name for name, coord in dataset.coords.items()]
    # print(coord_names)

    lat = next((dataset[key].values for key in lat_keys if key in coord_names), None)
    lon = next((dataset[key].values for key in lon_keys if key in coord_names), None)

    if lat is None or lon is None:
        raise KeyError("Latitude or longitude key not found in dataset.")

    return lat, lon


def interpolate_LAI_value(year, month, day, files_lai, eco_bound):
    """read 10day data and interploate to target date"""
    def ensure_2d_or_time0(da):
        """Assure que le DataArray est 2D, prend time=0 si nécessaire"""
        if 'time' in da.dims:
            return da.isel(time=0)
        return da

    # Step 1: Build DOY list
    obs_doys = []
    for m in range(1, 13):
        last_day = calendar.monthrange(year, m)[1]
        obs_days = [10, 20, last_day]
        for d in obs_days:
            dt = datetime(year, m, d)
            doy = dt.timetuple().tm_yday
            obs_doys.append((doy, dt))

    # Step 2: Find surrounding DOYs
    target_date = datetime(year, month, day)
    target_doy = target_date.timetuple().tm_yday
    doy_list = [doy for doy, _ in obs_doys]

    # Special cases at edges
    if target_doy <= 10:
        sel_date = datetime(year, 1, 10)
        path = next(p for p in files_lai if f"_{year:04d}0110" in p)
        return ensure_2d_or_time0(extract_from_global_data(path, eco_bound, variable='LAI'))

    if target_doy >= 364:
        candidates = [p for p in files_lai if f"_{year:04d}12" in p]
        if not candidates:
            raise ValueError(f"Aucun fichier LAI trouvé pour décembre {year}")
        path = sorted(candidates)[-1]
        return ensure_2d_or_time0(extract_from_global_data(path, eco_bound, variable='LAI'))

    before_doy = max([d for d in doy_list if d < target_doy], default=None)
    after_doy = min([d for d in doy_list if d > target_doy], default=None)
    if before_doy is None or after_doy is None:
        raise ValueError("Target date is outside the range of observation DOYs.")

    # Get corresponding dates
    before_date = next(dt for doy, dt in obs_doys if doy == before_doy)
    after_date  = next(dt for doy, dt in obs_doys if doy == after_doy)
    year_before, month_before, day_before = before_date.year, before_date.month, before_date.day
    year_after,  month_after,  day_after  = after_date.year,  after_date.month,  after_date.day

    # Step 3: Get file paths
    path_LAI_before = next(p for p in files_lai if f"_{year_before:04d}{month_before:02d}{day_before:02d}" in p)
    path_LAI_after  = next(p for p in files_lai if f"_{year_after:04d}{month_after:02d}{day_after:02d}" in p)

    # Step 4: Read data
    ds_before = ensure_2d_or_time0(extract_from_global_data(path_LAI_before, eco_bound, variable='LAI'))
    ds_after  = ensure_2d_or_time0(extract_from_global_data(path_LAI_after, eco_bound, variable='LAI'))

    ds_before.load()
    ds_after.load()

    # Step 5: Interpolate
    weight = (target_doy - before_doy) / (after_doy - before_doy)
    ds_interp = ds_before + weight * (ds_after - ds_before)

    # Step 6: Fill missing values
    ds_filled = ds_interp.fillna(ds_after).fillna(ds_before)

    return ds_filled


def interpolate_FVC_value(year, month, day, files_fvc, eco_bound):
    """read 10day data and interploate to target date"""
    def ensure_2d_or_time0(da):
        """Assure que le DataArray est 2D, prend time=0 si nécessaire"""
        if 'time' in da.dims:
            return da.isel(time=0)
        return da

    # Step 1: Build DOY list
    obs_doys = []
    for m in range(1, 13):
        last_day = calendar.monthrange(year, m)[1]
        obs_days = [10, 20, last_day]
        for d in obs_days:
            dt = datetime(year, m, d)
            doy = dt.timetuple().tm_yday
            obs_doys.append((doy, dt))

    # Step 2: Find surrounding DOYs
    target_date = datetime(year, month, day)
    target_doy = target_date.timetuple().tm_yday
    doy_list = [doy for doy, _ in obs_doys]

    # Special cases at edges
    if target_doy <= 10:
        sel_date = datetime(year, 1, 10)
        path = next(p for p in files_fvc if f"_{year:04d}0110" in p)
        return ensure_2d_or_time0(extract_from_global_data(path, eco_bound, variable='FCOVER'))

    if target_doy >= 364:
        candidates = [p for p in files_fvc if f"_{year:04d}12" in p]
        if not candidates:
            raise ValueError(f"Aucun fichier FCOVER trouvé pour décembre {year}")
        path = sorted(candidates)[-1]
        return ensure_2d_or_time0(extract_from_global_data(path, eco_bound, variable='FCOVER'))

    before_doy = max([d for d in doy_list if d < target_doy], default=None)
    after_doy = min([d for d in doy_list if d > target_doy], default=None)
    if before_doy is None or after_doy is None:
        raise ValueError("Target date is outside the range of observation DOYs.")

    # Get corresponding dates
    before_date = next(dt for doy, dt in obs_doys if doy == before_doy)
    after_date  = next(dt for doy, dt in obs_doys if doy == after_doy)
    year_before, month_before, day_before = before_date.year, before_date.month, before_date.day
    year_after,  month_after,  day_after  = after_date.year,  after_date.month,  after_date.day

    # Step 3: Get file paths
    path_LAI_before = next(p for p in files_fvc if f"_{year_before:04d}{month_before:02d}{day_before:02d}" in p)
    path_LAI_after  = next(p for p in files_fvc if f"_{year_after:04d}{month_after:02d}{day_after:02d}" in p)

    # Step 4: Read data
    ds_before = ensure_2d_or_time0(extract_from_global_data(path_LAI_before, eco_bound, variable='FCOVER'))
    ds_after  = ensure_2d_or_time0(extract_from_global_data(path_LAI_after, eco_bound, variable='FCOVER'))

    ds_before.load()
    ds_after.load()

    # Step 5: Interpolate
    weight = (target_doy - before_doy) / (after_doy - before_doy)
    ds_interp = ds_before + weight * (ds_after - ds_before)

    # Step 6: Fill missing values
    ds_filled = ds_interp.fillna(ds_after).fillna(ds_before)

    return ds_filled


def find_lai_path(files_lai, year, month, day):
    date_str = f"_{year:04d}{month:02d}{day:02d}"
    matches = [p for p in files_lai if date_str in p]
    if matches:
        return matches[0]
    else:
        # Try previous year
        year_prev = year - 1
        date_str_prev = f"_{year_prev:04d}{month:02d}{day:02d}"
        matches_prev = [p for p in files_lai if date_str_prev in p]
        if matches_prev:
            print(f"Warning: Falling back to previous year {year_prev} for {month:02d}-{day:02d}")
            return matches_prev[0]
        else:
            raise FileNotFoundError(f"No LAI file found for {month:02d}-{day:02d} in {year} or {year_prev}")


def read_oco2(year, month, day, eco_bound, files_oco2):
    """
    Read OCO-2 XCO2 data for a specific date, or fallback to closest year with same DOY.
    EEH covering 2018-2025, while oco2 temporal coverage:2018-08-09 to 2022-02-28
    
    Args:
        year, month, day (int): Target date.
        eco_bound (tuple): (lat_max, lat_min, lon_max, lon_min)
        files_oco2 (list): List of available OCO2 file paths.
    
    Returns:
        xarray.DataArray or None: Clipped XCO2 data, or None if not found.
    """
    target_date = datetime(year, month, day)
    target_date_str = target_date.strftime("%Y%m%d")
    target_doy = target_date.timetuple().tm_yday

    # Try exact match first
    exact_matches = [f for f in files_oco2 if target_date_str in f]
    if exact_matches:
        oco2_path = exact_matches[0]
    else:
        # Fallback: match same date across available years
        # doy_str = f"{target_doy:03d}"
        month_day_str = f"{month:02d}{day:02d}"
        date_matches = [f for f in files_oco2 if f"{month_day_str}_" in f]

        if not date_matches:
            print(f"No OCO2 file found for month {month:02d} day {day:02d}")
            return None

        # Extract year from filename and find closest to target year
        def extract_year(f):
            try:
                date_part = f.split("_day_")[1][:8]
                return int(date_part[:4])
            except:
                return None

        date_matches = [(f, extract_year(f)) for f in date_matches if extract_year(f)]
        date_matches.sort(key=lambda x: abs(x[1] - year))

        oco2_path = date_matches[0][0]

    # Load and clip
    try:
        ds = extract_from_global_data(oco2_path, eco_bound, variable='XCO2')
        return ds
    except Exception as e:
        print(f"Error reading OCO2 file: {e}")
        return None


def read_ci(year,month,eco_bound, files_CI):
    matches = [f for f in files_CI if f'Month{month}' in f]
    ci_path = matches[0] if matches else None
    ds = extract_from_global_data(ci_path,eco_bound,variable='band_data')
    ds = ds.where(ds != 65535, np.nan) #remove non data
    ds.data= ds.data/1000.0 # scale to actual value
    return ds


def plot_heatmap(data, title, label, ax=None):
    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 5))

    im = ax.imshow(data, cmap='viridis', aspect='auto')
    ax.set_title(title)
    plt.colorbar(im, ax=ax, label=label)
    return ax


# Define mapper
GLC30_mapper = {
    'CRO': 2, 'EBF': 5, 'DBF': 6, 'ENF': 7, 'DNF': 8, 'MF': 9,
    'SHR': 12, 'GRA': 13, 'SAV': 14, 'WSA':15, 'WET': 18, 'nonV': 0
}


# Define discrete colormap
def plot_lucc(data, title="LUCC Classification", label="Class", ax=None):
    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 5))

    # Reverse mapper for legend
    value_to_label = {v: k for k, v in GLC30_mapper.items()}
    # Ensure data is a NumPy array
    arr = data.values if isinstance(data, xr.DataArray) else data
    # Remove NaNs and get unique integer values
    unique_vals = np.unique(arr[~np.isnan(arr)]).astype(int)


    # Discrete colormap
    cmap = plt.get_cmap('viridis', len(unique_vals))
    norm = mcolors.BoundaryNorm(boundaries=np.append(unique_vals, unique_vals[-1]+1), ncolors=len(unique_vals))

    # Plot
    im = ax.imshow(arr, cmap=cmap, norm=norm, aspect='auto')
    ax.set_title(title)

    # Custom legend
    handles = [
        plt.Line2D([0], [0], marker='s', color='w',
                   label=value_to_label.get(val, str(val)),
                   markerfacecolor=cmap(i), markersize=10)
        for i, val in enumerate(unique_vals)
    ]
    ax.legend(handles=handles, bbox_to_anchor=(1.05, 1), loc='upper left', title="LUCC Types", frameon=False)
    return ax


def big_leaf_calculator(ds, vector):
    """
    Returns:
    - gpp in gC m-2 s-1
    """
    # TODO:
    umolPhoton_to_mj = 1/(4.57 * 1e6) # Convert from µmolPhoton-1 to MJ-1
    LAI = ds['LAI']
    PAR = ds['PARH']*umolPhoton_to_mj
    temperature = ds['Tair']
    vpd = ds['VPD']

    LUEmax, Tmax, Tmin, Topt, k = vector

    #Temperature constraint
    numerator = (temperature - Tmin) * (temperature - Tmax)
    denominator = numerator - (temperature - Topt) ** 2
    denominator = np.where(denominator==0,np.nan,denominator)
    stressor_temp = numerator / denominator

    stressor_temp[(temperature <= Tmin) | (temperature >= Tmax)] = 0

    stressor_temp = np.clip(stressor_temp, 0, 1)

    #VPD constraint
    VPDmin = 9
    VPDmax = 40
    stressor_vpd = (VPDmax - vpd) / (VPDmax - VPDmin)
    
    stressor_vpd[vpd <= VPDmin] = 1
    stressor_vpd[vpd >= VPDmax] = 0

    stressor_vpd = np.clip(stressor_vpd, 0, 1)

    fPAR = 1 - np.exp(-k * LAI)
    APAR = PAR * fPAR

    lue = LUEmax * stressor_temp * stressor_vpd
    gpp = APAR * lue

    return gpp


def gpp_gs_calculator(gt,CO2, T_C, VPD_hPa):

    """
    GPP = (1-ratio_ci_ca)*gt*CO2_ppm
    Calculate the intercellular CO2 concentration ratio (c_ratio) based on 
    Medlyn et al. (2011) using temperature, VPD, and enzymatic kinetics.
    T. F. Keenan (2023) CO2 based GPP https://doi.org/10.1038/s41558-023-01867-2
    
    Parameters:
    - T_C: temperature in clecius
    - CO2: concertration in mol dry area, need to convert to ppm(mju mol)
    - VPD_hPa: Vapor Pressure Deficit in hPa
    
    Returns:
    - GPP in gC
    
    Note: 
    global long-term average value of ci is  0.7, typical for C3 plants (Prentice et al., 2014; Wang et al., 2017). 
    
    """

    # Universal mole gas constant
    R = 8.314  #J mol−1 K−1
    mol_to_umol = 1e6
    mp_to_molm2 = 40.088/1.6
    umolCO2_to_gC =  12 * 1e-6 # Convert from µmolCO2 to gC
    CO2_ppm =CO2 *mol_to_umol


    # Michaelis-Menten coefficients for Rubisco (temperature dependencies of the kinetic properties at 25C)
    # r25 is 39.97 kPa for Kc and 27.48 kPa for Ko
    # ΔH is 79.43 kJ mol−1 for Kc and 36.38 kJ mol−1 for Ko
    T = T_C+273.15
    Kc = 39.97  * np.exp(79.43 * (T - 298.15) / (298.15 * R * T))  # CO2 coefficient
    Ko = 27.48  * np.exp(36.38 * (T - 298.15) / (298.15 * R * T))  # O2 coefficient
    Po = 21  # Partial pressure of O2 in kPa


    # Michaelis-Menten coefficient for Rubisco
    # K is the Michaelis–Menten coefficient for Rubisco-limited photosynthesis at a pO2 of 21 kPa.
    K = Kc * (1 + Po / Ko)   #in kPa

    # Water viscosity correction
    # n_star: Relative viscosity of water at given temp at 25 °C
    # referenceuber, M. L. et al. New international formulation for the viscosity of H2O. J. Phys. Chem. Ref. Data 38, 101–125 (2009)
    n_star = 0.8903  
    # ξ parameter from Medlyn et al. (2011) - also called g1
    # b is the ratio of the cost of maintaining carboxylation relative to that of maintaining transpiration
    # Prentice, I. C., Dong, N., Gleason, S. M., Maire, V. & Wright, I. J. Balancing the costs of carbon gain and water transport: testing a new theoretical framework for plant functional ecology. Ecol. Lett. 17, 82–91 (2014).
    b = 0.35651 # ratio from Farquhar et al(1980) 0.27 at 25degree
    Ksi = (b * K / (1.6 * n_star)) ** 0.5  

    # Calculate intercellular CO2 ratio (c_ratio)
    VPD_kPa = VPD_hPa/10
    c_ratio = Ksi / (Ksi + np.sqrt(VPD_kPa)) 
    
    # photopiratory compensation point that depends on temperature
    # # CO2 compensation point in the  absence of dark respiration
    # r25 = 4.22 Pa, is the photorespiratory point at 25  ̊C
    # ΔH is the activation energy for Γ* (37.83 kJ mol−1)
    # Universal mole gas constant
    T_K = T_C +273.15
    R = 8.314  #J mol−1 K−1
    T_star = 4.22*np.exp(37830 * (T_K - 298.15) / (298.15 * R * T_K)) #in Pa
    ci = c_ratio*CO2_ppm #ppm
    ci_ratio_refined = ((ci-T_star)/(ci+2*T_star)*c_ratio).clip(0.6, 1) 
    
    GPP = (1-ci_ratio_refined)*(gt*mp_to_molm2)*CO2_ppm*umolCO2_to_gC

    return GPP


#Output the estimates to HDF5 file
def generate_hdf5_file_GPP(directory_output,dict_outputs,key_eco):     
    #base_filename = f'EEH2/EEHGPP/EEH2STIC_L3_GPP_{key_eco}_0000_00.h5'
    base_filename = f'EEH2STIC_L3_GPP_{key_eco}_0000_00.h5'
    filename = os.path.join(directory_output,base_filename)
    wasdi.wasdiLog("Writing output to : "+filename)
    with h5py.File(filename, 'w') as f_h5:
        dset = f_h5.create_dataset('GPPd',data = np.float32(dict_outputs['GPPd']))
        dset.attrs['long_name'] = 'Daily gross primary productivity'
        dset.attrs['units'] = 'gC.m-2.day-1'
        dset.attrs['format'] = 'scaled'
        dset.attrs['coordsys'] = 'cartesian'
        dset.attrs['valid_range'] = np.array([0,50])
        dset.attrs['fill_value'] = -9999
        dset.attrs['scale_factor'] = 1
        dset.attrs['add_offset'] = 0

        dset = f_h5.create_dataset('WUEd',data = np.float32(dict_outputs['WUEd']))
        dset.attrs['long_name'] = 'Daily water use efficiency'
        dset.attrs['units'] = 'gC.m-2/mm.H2O'
        dset.attrs['format'] = 'scaled'
        dset.attrs['coordsys'] = 'cartesian'
        dset.attrs['valid_range'] = np.array([0,100])
        dset.attrs['fill_value'] = -9999
        dset.attrs['scale_factor'] = 1
        dset.attrs['add_offset'] = 0

        dset = f_h5.create_dataset('GPP_gs',data = np.float32(dict_outputs['GPP_gs']))
        dset.attrs['long_name'] = 'instantaneous STIC-gs derived GPP'
        dset.attrs['units'] = 'gC.m-2.s-1'
        dset.attrs['format'] = 'scaled'
        dset.attrs['coordsys'] = 'cartesian'
        dset.attrs['valid_range'] = np.array([0,0.001])
        dset.attrs['fill_value'] = -9999
        dset.attrs['scale_factor'] = 1
        dset.attrs['add_offset'] = 0

        dset = f_h5.create_dataset('GPP_LUE',data = np.float32(dict_outputs['GPP_LUE']))
        dset.attrs['long_name'] = 'instantaneous big-leaf model derived GPP (not considering VPD stress)'
        dset.attrs['units'] =  'gC.m-2.s-1'
        dset.attrs['format'] = 'scaled'
        dset.attrs['coordsys'] = 'cartesian'
        dset.attrs['valid_range'] = np.array([0,0.001])
        dset.attrs['fill_value'] = -9999
        dset.attrs['scale_factor'] = 1
        dset.attrs['add_offset'] = 0

        dset = f_h5.create_dataset('GPPi',data = np.float32(dict_outputs['GPPmin']))
        dset.attrs['long_name'] = 'STIC-gs calibrated BL-LUE GPP'
        dset.attrs['notes'] = 'min(GPP_gs, GPP_BL) if VPD_hPa>9 else GPP_BL'
        dset.attrs['units'] =  'gC.m-2.s-1'
        dset.attrs['format'] = 'scaled'
        dset.attrs['coordsys'] = 'cartesian'
        dset.attrs['valid_range'] = np.array([0,0.001])
        dset.attrs['fill_value'] = -9999
        dset.attrs['scale_factor'] = 1
        dset.attrs['add_offset'] = 0

        dset = f_h5.create_dataset('LUCC',data = np.int16(dict_outputs['LUCC']))
        dset.attrs['long_name'] = 'GLC_FCS30D'
        dset.attrs['class_dict']:  "{'CRO': 2, 'EBF': 5, 'DBF': 6, 'ENF': 7, 'DNF': 8, 'MF': 9,'SHR': 12, 'GRA': 13, 'SAV': 14, 'WSA':15, 'WET': 18, 'nonV': 0}"
        dset.attrs['units'] =  'unitless'
        dset.attrs['format'] = 'scaled'
        dset.attrs['coordsys'] = 'cartesian'
        dset.attrs['valid_range'] = np.array([0,18])
        dset.attrs['fill_value'] = 0
        dset.attrs['scale_factor'] = 1
        dset.attrs['add_offset'] = 0

    wasdi.wasdiLog("GPP Output Completed")
    
    return base_filename


