import sys
import os
import numpy as np
import h5py
import glob
import calendar
import pandas as pd
from datetime import datetime,timedelta
import xarray as xr
import re
import time
import seaborn as sns
import matplotlib.pyplot as plt
from shapely.geometry import Point, Polygon, box



def regex_extract_ecostress(filename):
    # Create a regex pattern to extract relevant details
    pattern = r".*_(\d+)_(\d+)_(\d{4})(\d{2})(\d{2})T(\d{2})(\d{2})(\d{2}).*"

    # Apply regex search
    match = re.search(pattern, filename)

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


def Create_footprint_mask(path_geo, row_site):
    """Create a footprint mask based on geolocation data."""
    
    with h5py.File(path_geo, "r") as f_geo:
        # Extract latitude and longitude information
        eco_lat =  np.array(f_geo['Geolocation']['latitude'])
        eco_lon =  np.array(f_geo['Geolocation']['longitude'])
   
    
    # Extract site coordinates
    lon, lat = row_site['LOCATION_LONG'], row_site['LOCATION_LAT']
    
    # Calculate footprint buffer size (forest = 3km, others = 1km)
    # buffer_degree = 0.0135 if 'F' in row_site['IGBP'] else 0.0045
    buffer_degree = 0.0045

    # Define bounding box for filtering
    min_x, min_y, max_x, max_y = Point(lon, lat).buffer(buffer_degree).bounds

    # Apply cloud and footprint mask
    mask_fp = (min_x <= eco_lon) & (eco_lon <= max_x) & (min_y <= eco_lat) & (eco_lat <= max_y) 
    
    # read time  
    # fp_time = ds_geo['Geolocation']['line_start_time_j2000'].sel(mask_ffp, method="nearest").values
    # linetime_mean = np.mean(fp_time)

    return mask_fp

def Create_cloud_mask(path_cld):
    """Create a cloud mask based on provided data."""
    
    with h5py.File(path_cld, "r") as f_cld:
        # Extract cloud
        eco_cld =  np.array(f_cld['SDS']['Cloud_final']) 
        # Create cloud mask (binary) 0=clear 1= cloud
        cmask = (eco_cld == 0).astype(np.int8)
        
    return cmask

def Extract_ECOSTRESS_from_mask(file_path, mask_ffp, variables=["LST"]):
    """Extract values from ECOSTRESS footprint using provided footprint and cloud masks."""
    
    # Open the HDF5 file using xarray with the appropriate engine
    ds = xr.open_dataset(file_path, engine="h5netcdf", phony_dims="sort")
    # print(ds)

    # Validate variables
    missing_vars = [var for var in variables if var not in ds.variables]
    if missing_vars:
        raise ValueError(f"Variables not found in dataset: {missing_vars}")
    mask_nan = ds[variables[0]]!=-9999 # overall mask 
    mask_cld_ffp = mask_ffp&mask_nan # remove invalid data
    valid_pixel_num = int(np.sum(mask_cld_ffp))
    # Extract values without nan
    return [ds[var].where(mask_cld_ffp&(ds[var]!=-9999)).values.astype(float) for var in variables],valid_pixel_num

def Extract_LSTE(filename_lste, mask_cld_ffp):
    """Extract and process LST and BBE values from ECOSTRESS dataset."""
    variables = ["LST", "BBE"]
    if np.sum(mask_cld_ffp) > 0:  # Ensuring valid data within footprint&cloud mask
        value_list,valid_pixel_num = Extract_ECOSTRESS_from_mask(filename_lste, mask_cld_ffp,variables)
        mean_list = []
        for key,img in zip(variables,value_list):
            value_mean =  np.nanmean(img) if ~np.isnan(img).all() else np.nan
            mean_list.append(value_mean)
        return mean_list  + [valid_pixel_num]
    else:
        return [np.nan for i in range(len(variables))] + [valid_pixel_num]


def Extract_GPP(filename_GPP,mask_fp):
    """Extract and process GPP values from EEH2 dataset."""
    # nan=-9999
    variables = ["GPPd", "WUEd",'GPP_gs','GPP_LUE','GPPi']
    if np.sum(mask_fp) > 0:  # Ensuring valid data within footprint&cloud mask
        value_list,valid_pixel_num = Extract_ECOSTRESS_from_mask(filename_GPP,mask_fp,variables)
        mean_list = []
        for key,img in zip(variables,value_list):
            value_mean =  np.nanmean(img) if ~np.isnan(img).all() else np.nan
            mean_list.append(value_mean)
        valid_pixel_num
        return mean_list + [valid_pixel_num]
    else:
        return [np.nan for i in range(len(variables))] + [valid_pixel_num]
    
def timetransform(time):
    t = datetime(2000, 1, 1, 12, 0) + timedelta(seconds=time)
    
    year = t.strftime('%Y')
    month = t.strftime('%m')
    day = t.strftime('%d')
    hour = t.strftime('%H')
    minute = t.strftime('%M')
    second = t.strftime('%S')   

    return (year,month,day,hour,minute,second)

def run_extraction(df_site_with_orbits,files_GPP,files_lste,files_geo,directory_output):
    # Initialize list to store extracted rows
    data_list = []
    for _, site_info in df_site_with_orbits.iterrows():
        orbit_list = [str(id).strip() for id in str(site_info.orbitID_list).split(',')]

        # Filter relevant LSTE files
        files_GPP_with_sites = [f for f in files_GPP if any(orbit_id in f for orbit_id in orbit_list)]
        files_LSTE_with_sites = [f for f in files_lste if any(orbit_id in f for orbit_id in orbit_list)]
        if not files_GPP_with_sites:
            continue

        for name_gpp in files_GPP_with_sites:  # Process only the first match
            # Site information
            name_site, lat_site, lon_site,maxYear,minYear = site_info.SITE, site_info.LOCATION_LAT, site_info.LOCATION_LONG,site_info.maxYear,site_info.minYear
            print(name_site,maxYear,minYear,name_gpp)

            # Extract orbit & time details
            orbit_str, scene_str, year_str, month_str, day_str, hour_str, min_str, sec_str = regex_extract_ecostress(name_gpp)
            orbit_id = f"{orbit_str}_{scene_str}"

            if (int(year_str)<minYear)|(int(year_str)>maxYear):
                continue
            # Compose filenames for ECOSTRESS data
            key_eco = f"{orbit_str}_{scene_str}_{year_str}{month_str}{day_str}T{hour_str}{min_str}{sec_str}"
            path_geo = next((f for f in files_geo if key_eco in f), None)
            # path_cld = next((f for f in files_cld if key_eco in f), None)
            path_GPP = next((f for f in files_GPP if key_eco in f), None)
            path_LSTE = next((f for f in files_lste if key_eco in f), None)

            # Early exit if required files are missing
            if not all([path_geo, path_GPP,path_LSTE]):
                print("Corresponding data not available")
                continue

            try:
                start_time = time.time()
                mask_fp = Create_footprint_mask(path_geo, site_info) #3km for forest, 1km for non-forest
                print(f"Reading GEO file done in {time.time() - start_time:.2f} seconds")
                if np.sum(mask_fp) <30:
                    print(f'site is not located within the tile or at the boundary')
                    continue
            except Exception as e:
                print(f"Error processing GEO file: {e}")
                continue

            try:
                start_time = time.time()
                # mask_fp_cld = mask_fp & mask_cld # cloud and footprint mask
                GPP_list = Extract_GPP(path_GPP, mask_fp)
                LSTE_list = Extract_LSTE(path_LSTE, mask_fp)
                print(f"Extracting LSTE done in {time.time() - start_time:.2f} seconds")
            except Exception as e:
                print(f"Error extracting LSTE: {e}")
                continue

            # Validate extracted LST data before storing
            if not np.isnan(GPP_list[0]):
                row_data = [name_site, orbit_id, int(year_str), int(month_str), int(day_str),
                            int(hour_str), int(min_str), int(sec_str)]+GPP_list+LSTE_list
                data_list.append(row_data)
                print(row_data)
            clear_output(wait=True)

            # Convert to DataFrame and save results
            columns = ['SITE', 'orbitID', 'Year', 'Month', 'Day', 'Hour', 'Minute', 'Second']+ \
             ["GPPd", "WUEd",'GPP_gs','GPP_LUE','GPPi','GPPd_pixel_num']+['LST','BBE','LST_pixel_num']
        df_eco = pd.DataFrame(data_list, columns=columns)
        df_eco.to_csv(directory_output + "ECOSTRESS_GPP_extraction.csv", index=False)



if __name__ == '__main__':

    #STEP1 Read paths on S3
    # Read paths of GPP product with sites  
    directory_GPP = r'/ECOSTRESS_RW/EEH2/EEHGPP-final/*.h5'
    files_GPP = glob.glob(directory_GPP)

    # Read the list of file paths
    files_geo = [] ??
    files_lste = [] ??
    
    # Load site data
    df_site_with_orbits = pd.read_csv("Sitelist_ICOS_with_orbits-for-EEH.csv")
    directory_output = ??
    #could change df_site_with_orbits to single site dataframe?
    # -->df_single_site_info = df_site_with_orbits[df_site_with_orbits.SITE==site_name]
    run_extraction(df_site_with_orbits,files_GPP,files_lste,files_geo,directory_output)