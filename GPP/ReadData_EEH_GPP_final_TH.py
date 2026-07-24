#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Script to run the hybrid model for GPP and WUE estimation
Created on May 1 2025
@author: Ziyu Lin, Kaniska Mallick, Tian Hu, Yoanne Didry at LIST

© 2026 - Luxembourg Institute of Science and Technology
Authors : Ziyu Lin, Kaniska Mallick, Tian Hu (tian.hu@list.lu), and Yoanne Didry (yoanne.didry@list.lu)
Code licensed under MIT
SPDX-License-Identifier: MIT
"""

# Supporting functions

import os
import re
import calendar
from datetime import datetime, timedelta

import numpy as np
import h5py
import pandas as pd
import xarray as xr
import rioxarray
from rioxarray.merge import merge_arrays
from rasterio.enums import Resampling
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from shapely.geometry import Point, Polygon, box
import geopandas as gpd
from pyhdf.SD import SD, SDC
from scipy.interpolate import griddata

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
            'band': ['BSA', 'WSA', 'SZA', 'QA'],
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
    cos_sza = np.clip(np.cos(SZA_rad), 1e-6, None)  # mask non-positive values
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
    """
    lat_src_2d, lon_src_2d = get_lat_lon(dataset)
    data_src = dataset.values

    points = np.column_stack((lat_src_2d.ravel(), lon_src_2d.ravel()))
    values = data_src.ravel()
    target_points = np.column_stack((lat_eco.ravel(), lon_eco.ravel()))
    interpolated = griddata(points, values, target_points, method='linear')
    interpolated_2d = interpolated.reshape(lat_eco.shape)

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

    lon_dir = match.group(1)
    lon = int(match.group(2))
    lat_dir = match.group(3)
    lat = int(match.group(4))
    tile_size = 5  # 5*5 degree per tile
    lat = lat if lat_dir == 'N' else -lat
    lon = lon if lon_dir == 'E' else -lon
    t_min_lon = lon
    t_max_lon = lon + tile_size
    t_max_lat = lat
    t_min_lat = lat - tile_size

    e_max_lat, e_min_lat, e_max_lon, e_min_lon = eco_bound

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
    """
    lat_max, lat_min, lon_max, lon_min = eco_bound
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
            da = rioxarray.open_rasterio(path, masked=True)
            da_band = da.isel(band=band_ind)

            da_clipped = da_band.rio.clip_box(
                minx=lon_min, miny=lat_min, maxx=lon_max, maxy=lat_max,
                auto_expand=False
            )

            scaled = (da_clipped.data / 10).astype(np.int16)
            nonveg_mask = np.isin(scaled, [14, 15, 19, 20, 21, 22])
            scaled[nonveg_mask] = 0
            scaled[scaled == 1] = 2  # CRO
            da_clipped.data = scaled

            da_aligned = da_clipped.interp(
                x=glc30_grid['x'],
                y=glc30_grid['y'],
                method="nearest"
            )
            valid = (glc30_grid >= 0)
            glc30_grid = xr.where(valid, glc30_grid, da_aligned)

        except Exception as e:
            print(f"Error processing {path}: {e}")

        # SAV: Grass and Shrub found in Africa between latitude 15N and 30S
        scaled = glc30_grid.data
        scaled[(glc30_grid['y'].values < 15) & (glc30_grid['y'].values > -30) & (scaled == 13)] = 14  # grass SAV
        scaled[(glc30_grid['y'].values < 15) & (glc30_grid['y'].values > -30) & (scaled == 12)] = 15  # woody WSA
        glc30_grid.data = scaled

    da.close()

    return glc30_grid


def Create_cloud_mask(path_cld):
    """Create a mask to select non-cloudy pixels based on provided data."""
    with h5py.File(path_cld, "r") as f_cld:
        eco_cld = np.array(f_cld['SDS']['Cloud_final'])
        # FOR v001 data 0=clear 1=cloud -1 or 255=missing
        cmask = (eco_cld == 0).astype(np.int8)
    return cmask


def Extract_ECOSTRESS_from_mask(file_path, mask_cld_ffp, variables=["LST"]):
    """Extract values from ECOSTRESS footprint using provided footprint and cloud masks."""
    ds = xr.open_dataset(file_path, engine="h5netcdf", phony_dims="sort")
    missing_vars = [var for var in variables if var not in ds.variables]
    if missing_vars:
        raise ValueError(f"Variables not found in dataset: {missing_vars}")
    return [ds[var].values.astype(float) for var in variables]


def Extract_STIC(filename_lste, mask_cld_ffp):
    """Extract STIC outputs values of EEH2"""
    variables = ["ETD", "gah", "gsc"]
    if np.sum(mask_cld_ffp) > 0:
        ETD, gah, gsc = Extract_ECOSTRESS_from_mask(filename_lste, mask_cld_ffp, variables)
        gt = gsc * gah / (gsc + gah)
        ETD_filtered = np.where(mask_cld_ffp & (ETD < 20) & (ETD > 0), ETD, np.nan)
        gt_filtered = np.where(mask_cld_ffp & (gah > 0) & (gsc > 0), gt, np.nan).clip(0.0001, 0.2)
        return [ETD_filtered, gt_filtered]
    else:
        nan_array = np.full_like(mask_cld_ffp, np.nan, dtype=float)
        return [nan_array for _ in range(len(variables))]


def read_geo(path_geo):
    """read geolocation data."""
    with h5py.File(path_geo, "r") as f_geo:
        eco_lat = np.array(f_geo['Geolocation']['latitude'])
        eco_lon = np.array(f_geo['Geolocation']['longitude'])
        sza_lon = np.array(f_geo['Geolocation']['solar_zenith'])  # in degree
    return eco_lat, eco_lon, sza_lon


def timetransform(time):
    t = datetime(2000, 1, 1, 12, 0) + timedelta(seconds=time)
    year = t.strftime('%Y')
    month = t.strftime('%m')
    day = t.strftime('%d')
    hour = t.strftime('%H')
    minute = t.strftime('%M')
    second = t.strftime('%S')
    return (year, month, day, hour, minute, second)


def regex_extract_ecostress(filename):
    pattern = r".*_(\d+)_(\d+)_(\d{4})(\d{2})(\d{2})T(\d{2})(\d{2})(\d{2}).*"
    match = re.search(pattern, filename)
    print(f"Searching {pattern} in {filename}")
    if match:
        orbit_str = match.group(1)
        scene_str = match.group(2)
        year_str = match.group(3)
        month_str = match.group(4)
        day_str = match.group(5)
        hour_str = match.group(6)
        min_str = match.group(7)
        sec_str = match.group(8)
        return (orbit_str, scene_str, year_str, month_str, day_str, hour_str, min_str, sec_str)
    else:
        print("No match found in file name.")
        return None


def extract_from_global_data(path_globe, eco_bound, variable='LAI'):
    """
    Extract and interpolate geospatial data from a global raster/NetCDF file
    within a specified bounding box.
    """
    if path_globe.endswith(('.nc', '.nc4', '.h5', '.hdf5', '.cdf')):
        ds = xr.open_dataset(path_globe, decode_timedelta=False, engine='netcdf4')
    elif path_globe.endswith('.grib'):
        ds = xr.open_dataset(path_globe, decode_timedelta=False, engine='cfgrib',
                              backend_kwargs={"errors": "ignore"})
    elif path_globe.endswith('.tif'):
        ds = xr.open_dataset(path_globe, decode_timedelta=False, engine="rasterio")
    else:
        raise ValueError("Only .nc,.nc4,.h5,.hdf5,.cdf,.grib,.tif files are supported.")

    ds_data = ds[variable]
    lat_name = [name for name in ds_data.coords if ('lat' in name) | (name == 'y')][0]
    lon_name = [name for name in ds_data.coords if ('lon' in name) | (name == 'x')][0]
    dims_2d = (lat_name, lon_name)
    ds_data.rio.set_spatial_dims(x_dim=lon_name, y_dim=lat_name, inplace=True)
    ds_data.rio.write_crs("EPSG:4326", inplace=True)

    lat_max, lat_min, lon_max, lon_min = eco_bound
    da_clipped = ds_data.rio.clip_box(
        minx=lon_min, miny=lat_min, maxx=lon_max, maxy=lat_max,
        auto_expand=True
    )

    lat_1d, lon_1d = da_clipped[lat_name].values, da_clipped[lon_name].values
    lon_2d, lat_2d = np.meshgrid(lon_1d, lat_1d)
    da_clipped = da_clipped.assign_coords({
        lat_name: (dims_2d, lat_2d),
        lon_name: (dims_2d, lon_2d)
    })
    da_clipped = da_clipped.load()
    ds.close()
    return da_clipped


def extract_PAR_from_global_data(path_globe, eco_bound, hour, minute):
    """
    Extract and interpolate PAR data from a global NetCDF file within a
    specified bounding box, matching the overpass time.
    """
    if path_globe.endswith('.nc'):
        ds = xr.open_dataset(path_globe, decode_timedelta=False, engine='netcdf4')
        ds_data = ds['PAR']
        lon_1d = ds.lon.values
        lat_1d = ds.lat.values
        lon_2d, lat_2d = np.meshgrid(lon_1d, lat_1d)
        ds_data = ds_data.assign_coords(lat=(["lat", "lon"], lat_2d),
                                         lon=(["lat", "lon"], lon_2d))
    else:
        raise ValueError("Only .nc files are supported.")

    lat_max, lat_min, lon_max, lon_min = eco_bound
    mask_bound = (lat_2d >= lat_min) & (lat_2d <= lat_max) & \
                 (lon_2d >= lon_min) & (lon_2d <= lon_max)
    time_len = ds_data.sizes['time']
    mask_3d = np.broadcast_to(mask_bound, (time_len, *mask_bound.shape))
    mask_da = xr.DataArray(
        mask_3d,
        dims=ds_data.dims,
        coords={
            'time': ds_data['time'],
            'lat': ds_data['lat'],
            'lon': ds_data['lon']
        }
    )
    da_clipped = ds_data.where(mask_da, drop=True)

    time_index = int(hour * 2 + np.round(minute / 30))
    da_PARH = da_clipped.isel(time=time_index)
    da_PARmean = da_clipped.mean(dim='time')
    return da_PARH, da_PARmean


def interpolate_ERA5_value(year, month, day, hour, minute, files_era5, eco_bound):
    target_time = datetime(year, month, day, hour, minute)
    before_path = None
    after_path = None
    before_time = None
    after_time = None
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

    ds_before = extract_from_global_data(before_path, eco_bound, variable=['t2m', 'd2m'])
    ds_after = extract_from_global_data(after_path, eco_bound, variable=['t2m', 'd2m'])

    weight = (target_time - before_time) / (after_time - before_time)
    ds_interp = ds_before + weight * (ds_after - ds_before)
    ds_filled = ds_interp.fillna(ds_after).fillna(ds_before)

    t2m = ds_filled['t2m']
    t2m_C = t2m - 273.15
    d2m = ds_filled['d2m']
    vpd_hPa = calculate_vpd_from_t2m_d2m(t2m, d2m)
    return t2m_C, vpd_hPa


def calculate_vpd_from_t2m_d2m(t2m, d2m):
    """
    Calculate VPD from 2-meter temperature and dew point temperature.

    Parameters:
        t2m : Air temperature at 2 meters (Kelvin)
        d2m : Dew point temperature at 2 meters (Kelvin)
    Returns:
        vpd : Vapor Pressure Deficit in hPa
    """
    T_C = t2m - 273.15
    TD_C = d2m - 273.15
    e_s = 610.94 * np.exp((17.625 * T_C) / (T_C + 243.04))
    e_a = 610.94 * np.exp((17.625 * TD_C) / (TD_C + 243.04))
    vpd = (e_s - e_a) / 100.0
    return vpd


def get_lat_lon(dataset):
    lat_keys = ['y', 'lat', 'latitude']
    lon_keys = ['x', 'lon', 'longitude']
    coord_names = [name for name, coord in dataset.coords.items()]
    lat = next((dataset[key].values for key in lat_keys if key in coord_names), None)
    lon = next((dataset[key].values for key in lon_keys if key in coord_names), None)
    if lat is None or lon is None:
        raise KeyError("Latitude or longitude key not found in dataset.")
    return lat, lon


def interpolate_LAI_value(year, month, day, files_lai, eco_bound):
    """read 10day data and interpolate to target date"""
    def ensure_2d_or_time0(da):
        if 'time' in da.dims:
            return da.isel(time=0)
        return da

    obs_doys = []
    for m in range(1, 13):
        last_day = calendar.monthrange(year, m)[1]
        obs_days = [10, 20, last_day]
        for d in obs_days:
            dt = datetime(year, m, d)
            doy = dt.timetuple().tm_yday
            obs_doys.append((doy, dt))

    target_date = datetime(year, month, day)
    target_doy = target_date.timetuple().tm_yday
    doy_list = [doy for doy, _ in obs_doys]

    if target_doy <= 10:
        path = next(p for p in files_lai if f"_{year:04d}0110" in p)
        return ensure_2d_or_time0(extract_from_global_data(path, eco_bound, variable='LAI'))
    if target_doy >= 364:
        candidates = [p for p in files_lai if f"_{year:04d}12" in p]
        if not candidates:
            raise ValueError(f"Aucun fichier LAI trouve pour decembre {year}")
        path = sorted(candidates)[-1]
        return ensure_2d_or_time0(extract_from_global_data(path, eco_bound, variable='LAI'))
    before_doy = max([d for d in doy_list if d < target_doy], default=None)
    after_doy = min([d for d in doy_list if d > target_doy], default=None)
    if before_doy is None or after_doy is None:
        raise ValueError("Target date is outside the range of observation DOYs.")

    before_date = next(dt for doy, dt in obs_doys if doy == before_doy)
    after_date = next(dt for doy, dt in obs_doys if doy == after_doy)
    year_before, month_before, day_before = before_date.year, before_date.month, before_date.day
    year_after, month_after, day_after = after_date.year, after_date.month, after_date.day

    path_LAI_before = next(p for p in files_lai if f"_{year_before:04d}{month_before:02d}{day_before:02d}" in p)
    path_LAI_after = next(p for p in files_lai if f"_{year_after:04d}{month_after:02d}{day_after:02d}" in p)

    ds_before = ensure_2d_or_time0(extract_from_global_data(path_LAI_before, eco_bound, variable='LAI'))
    ds_after = ensure_2d_or_time0(extract_from_global_data(path_LAI_after, eco_bound, variable='LAI'))
    ds_before.load()
    ds_after.load()

    weight = (target_doy - before_doy) / (after_doy - before_doy)
    ds_interp = ds_before + weight * (ds_after - ds_before)
    ds_filled = ds_interp.fillna(ds_after).fillna(ds_before)
    return ds_filled


def interpolate_FVC_value(year, month, day, files_fvc, eco_bound):
    """read 10day data and interpolate to target date"""
    def ensure_2d_or_time0(da):
        if 'time' in da.dims:
            return da.isel(time=0)
        return da

    obs_doys = []
    for m in range(1, 13):
        last_day = calendar.monthrange(year, m)[1]
        obs_days = [10, 20, last_day]
        for d in obs_days:
            dt = datetime(year, m, d)
            doy = dt.timetuple().tm_yday
            obs_doys.append((doy, dt))

    target_date = datetime(year, month, day)
    target_doy = target_date.timetuple().tm_yday
    doy_list = [doy for doy, _ in obs_doys]

    if target_doy <= 10:
        path = next(p for p in files_fvc if f"_{year:04d}0110" in p)
        return ensure_2d_or_time0(extract_from_global_data(path, eco_bound, variable='FCOVER'))
    if target_doy >= 364:
        candidates = [p for p in files_fvc if f"_{year:04d}12" in p]
        if not candidates:
            raise ValueError(f"Aucun fichier FCOVER trouve pour decembre {year}")
        path = sorted(candidates)[-1]
        return ensure_2d_or_time0(extract_from_global_data(path, eco_bound, variable='FCOVER'))
    before_doy = max([d for d in doy_list if d < target_doy], default=None)
    after_doy = min([d for d in doy_list if d > target_doy], default=None)
    if before_doy is None or after_doy is None:
        raise ValueError("Target date is outside the range of observation DOYs.")

    before_date = next(dt for doy, dt in obs_doys if doy == before_doy)
    after_date = next(dt for doy, dt in obs_doys if doy == after_doy)
    year_before, month_before, day_before = before_date.year, before_date.month, before_date.day
    year_after, month_after, day_after = after_date.year, after_date.month, after_date.day

    path_LAI_before = next(p for p in files_fvc if f"_{year_before:04d}{month_before:02d}{day_before:02d}" in p)
    path_LAI_after = next(p for p in files_fvc if f"_{year_after:04d}{month_after:02d}{day_after:02d}" in p)

    ds_before = ensure_2d_or_time0(extract_from_global_data(path_LAI_before, eco_bound, variable='FCOVER'))
    ds_after = ensure_2d_or_time0(extract_from_global_data(path_LAI_after, eco_bound, variable='FCOVER'))
    ds_before.load()
    ds_after.load()

    weight = (target_doy - before_doy) / (after_doy - before_doy)
    ds_interp = ds_before + weight * (ds_after - ds_before)
    ds_filled = ds_interp.fillna(ds_after).fillna(ds_before)
    return ds_filled


def find_lai_path(files_lai, year, month, day):
    date_str = f"_{year:04d}{month:02d}{day:02d}"
    matches = [p for p in files_lai if date_str in p]
    if matches:
        return matches[0]
    else:
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
    """
    target_date = datetime(year, month, day)
    target_date_str = target_date.strftime("%Y%m%d")
    exact_matches = [f for f in files_oco2 if target_date_str in f]
    if exact_matches:
        oco2_path = exact_matches[0]
    else:
        month_day_str = f"{month:02d}{day:02d}"
        date_matches = [f for f in files_oco2 if f"{month_day_str}_" in f]
        if not date_matches:
            print(f"No OCO2 file found for month {month:02d} day {day:02d}")
            return None

        def extract_year(f):
            try:
                date_part = f.split("_day_")[1][:8]
                return int(date_part[:4])
            except Exception:
                return None

        date_matches = [(f, extract_year(f)) for f in date_matches if extract_year(f)]
        date_matches.sort(key=lambda x: abs(x[1] - year))
        oco2_path = date_matches[0][0]
    try:
        ds = extract_from_global_data(oco2_path, eco_bound, variable='XCO2')
        return ds
    except Exception as e:
        print(f"Error reading OCO2 file: {e}")
        return None


def read_ci(year, month, eco_bound, files_CI):
    matches = [f for f in files_CI if f'Month{month}' in f]
    ci_path = matches[0] if matches else None
    ds = extract_from_global_data(ci_path, eco_bound, variable='band_data')
    ds = ds.where(ds != 65535, np.nan)
    ds.data = ds.data / 1000.0
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
    'SHR': 12, 'GRA': 13, 'SAV': 14, 'WSA': 15, 'WET': 18, 'nonV': 0
}


def plot_lucc(data, title="LUCC Classification", label="Class", ax=None):
    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 5))
    value_to_label = {v: k for k, v in GLC30_mapper.items()}
    arr = data.values if isinstance(data, xr.DataArray) else data
    unique_vals = np.unique(arr[~np.isnan(arr)]).astype(int)
    cmap = plt.get_cmap('viridis', len(unique_vals))
    norm = mcolors.BoundaryNorm(boundaries=np.append(unique_vals, unique_vals[-1] + 1), ncolors=len(unique_vals))
    im = ax.imshow(arr, cmap=cmap, norm=norm, aspect='auto')
    ax.set_title(title)
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
        gpp in gC m-2 s-1
    """
    umolPhoton_to_mj = 1 / (4.57 * 1e6)  # Convert from µmolPhoton-1 to MJ-1
    LAI = ds['LAI']
    PAR = ds['PARH'] * umolPhoton_to_mj
    temperature = ds['Tair']
    vpd = ds['VPD']
    LUEmax, Tmax, Tmin, Topt, k = vector

    # Temperature constraint
    numerator = (temperature - Tmin) * (temperature - Tmax)
    denominator = numerator - (temperature - Topt) ** 2
    denominator = np.where(denominator == 0, np.nan, denominator)
    stressor_temp = numerator / denominator
    stressor_temp[(temperature <= Tmin) | (temperature >= Tmax)] = 0
    stressor_temp = np.clip(stressor_temp, 0, 1)

    # VPD constraint
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


def gpp_gs_calculator(gt, CO2, T_C, VPD_hPa):
    """
    GPP = (1-ratio_ci_ca)*gt*CO2_ppm
    Calculate the intercellular CO2 concentration ratio (c_ratio) based on
    Medlyn et al. (2011) using temperature, VPD, and enzymatic kinetics.
    T. F. Keenan (2023) CO2 based GPP https://doi.org/10.1038/s41558-023-01867-2

    Parameters:
        T_C: temperature in celsius
        CO2: concentration in mol dry area, need to convert to ppm (µmol)
        VPD_hPa: Vapor Pressure Deficit in hPa

    Returns:
        GPP in gC

    Note:
        global long-term average value of ci is 0.7, typical for C3 plants
        (Prentice et al., 2014; Wang et al., 2017).
    """
    R = 8.314  # J mol-1 K-1
    mol_to_umol = 1e6
    mp_to_molm2 = 40.088 / 1.6
    umolCO2_to_gC = 12 * 1e-6
    CO2_ppm = CO2 * mol_to_umol

    T = T_C + 273.15
    Kc = 39.97 * np.exp(79.43 * (T - 298.15) / (298.15 * R * T))  # CO2 coefficient
    Ko = 27.48 * np.exp(36.38 * (T - 298.15) / (298.15 * R * T))  # O2 coefficient
    Po = 21  # kPa
    K = Kc * (1 + Po / Ko)  # kPa

    n_star = 0.8903
    b = 0.35651
    Ksi = (b * K / (1.6 * n_star)) ** 0.5

    VPD_kPa = VPD_hPa / 10
    c_ratio = Ksi / (Ksi + np.sqrt(VPD_kPa))

    T_K = T_C + 273.15
    T_star = 4.22 * np.exp(37830 * (T_K - 298.15) / (298.15 * R * T_K))  # Pa
    ci = c_ratio * CO2_ppm
    ci_ratio_refined = ((ci - T_star) / (ci + 2 * T_star) * c_ratio).clip(0.6, 1)

    GPP = (1 - ci_ratio_refined) * (gt * mp_to_molm2) * CO2_ppm * umolCO2_to_gC
    return GPP


def generate_hdf5_file_GPP(directory_output, dict_outputs, key_eco):
    base_filename = f'EEH2STIC_L3_GPP_{key_eco}_0000_00.h5'
    filename = os.path.join(directory_output, base_filename)
    print("Writing output to : " + filename)
    with h5py.File(filename, 'w') as f_h5:
        dset = f_h5.create_dataset('GPPd', data=np.float32(dict_outputs['GPPd']))
        dset.attrs['long_name'] = 'Daily gross primary productivity'
        dset.attrs['units'] = 'gC.m-2.day-1'
        dset.attrs['format'] = 'scaled'
        dset.attrs['coordsys'] = 'cartesian'
        dset.attrs['valid_range'] = np.array([0, 50])
        dset.attrs['fill_value'] = -9999
        dset.attrs['scale_factor'] = 1
        dset.attrs['add_offset'] = 0

        dset = f_h5.create_dataset('WUEd', data=np.float32(dict_outputs['WUEd']))
        dset.attrs['long_name'] = 'Daily water use efficiency'
        dset.attrs['units'] = 'gC.m-2/mm.H2O'
        dset.attrs['format'] = 'scaled'
        dset.attrs['coordsys'] = 'cartesian'
        dset.attrs['valid_range'] = np.array([0, 100])
        dset.attrs['fill_value'] = -9999
        dset.attrs['scale_factor'] = 1
        dset.attrs['add_offset'] = 0

        dset = f_h5.create_dataset('GPP_gs', data=np.float32(dict_outputs['GPP_gs']))
        dset.attrs['long_name'] = 'instantaneous STIC-gs derived GPP'
        dset.attrs['units'] = 'gC.m-2.s-1'
        dset.attrs['format'] = 'scaled'
        dset.attrs['coordsys'] = 'cartesian'
        dset.attrs['valid_range'] = np.array([0, 0.001])
        dset.attrs['fill_value'] = -9999
        dset.attrs['scale_factor'] = 1
        dset.attrs['add_offset'] = 0

        dset = f_h5.create_dataset('GPP_LUE', data=np.float32(dict_outputs['GPP_LUE']))
        dset.attrs['long_name'] = 'instantaneous big-leaf model derived GPP (not considering VPD stress)'
        dset.attrs['units'] = 'gC.m-2.s-1'
        dset.attrs['format'] = 'scaled'
        dset.attrs['coordsys'] = 'cartesian'
        dset.attrs['valid_range'] = np.array([0, 0.001])
        dset.attrs['fill_value'] = -9999
        dset.attrs['scale_factor'] = 1
        dset.attrs['add_offset'] = 0

        dset = f_h5.create_dataset('GPPi', data=np.float32(dict_outputs['GPPmin']))
        dset.attrs['long_name'] = 'STIC-gs calibrated BL-LUE GPP'
        dset.attrs['notes'] = 'min(GPP_gs, GPP_BL) if VPD_hPa>9 else GPP_BL'
        dset.attrs['units'] = 'gC.m-2.s-1'
        dset.attrs['format'] = 'scaled'
        dset.attrs['coordsys'] = 'cartesian'
        dset.attrs['valid_range'] = np.array([0, 0.001])
        dset.attrs['fill_value'] = -9999
        dset.attrs['scale_factor'] = 1
        dset.attrs['add_offset'] = 0

        dset = f_h5.create_dataset('LUCC', data=np.int16(dict_outputs['LUCC']))
        dset.attrs['long_name'] = 'GLC_FCS30D'
        dset.attrs['class_dict'] = "{'CRO': 2, 'EBF': 5, 'DBF': 6, 'ENF': 7, 'DNF': 8, 'MF': 9,'SHR': 12, 'GRA': 13, 'SAV': 14, 'WSA':15, 'WET': 18, 'nonV': 0}"
        dset.attrs['units'] = 'unitless'
        dset.attrs['format'] = 'scaled'
        dset.attrs['coordsys'] = 'cartesian'
        dset.attrs['valid_range'] = np.array([0, 18])
        dset.attrs['fill_value'] = 0
        dset.attrs['scale_factor'] = 1
        dset.attrs['add_offset'] = 0
    print("GPP Output Completed")

    return base_filename
