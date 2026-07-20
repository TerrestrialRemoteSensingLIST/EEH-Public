import sys
sys.path.insert(0,'./utils')
from ReadData_EEH_GPP_final import *
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
from shapely.geometry import Point, Polygon, box
import geopandas as gpd
import re
import sys
from pyhdf.SD import SD, SDC
from scipy.interpolate import griddata

import matplotlib
matplotlib.use('Agg')

import wasdi
from S3_cache import *

def refine_GPP(gpp_file,writepath,map_cld,map_geo,
stic_path,files_PARH,files_fvc,files_era5,files_oco2,directory_cld,directory_geo):

    wasdi.wasdiLog(gpp_file)

    # Extract orbit & time details
    orbit_str, scene_str, year_str, month_str, day_str, hour_str, min_str, sec_str = regex_extract_ecostress(gpp_file)
    year, month, day, hour, minute = int(year_str), int(month_str), int(day_str), int(hour_str), int(min_str)   
    doy = datetime(year, month, day).timetuple().tm_yday
    doy_str = str(int(doy)) 

    # Compose filenames for ECOSTRESS data
    key_eco = f"{orbit_str}_{scene_str}_{year_str}{month_str}{day_str}T{hour_str}{min_str}{sec_str}"
    orbit_id = f"{orbit_str}_{scene_str}"

    path_geo = os.path.join(directory_geo,map_geo['ECOv002_L1B_GEO_'+key_eco])
    path_cld = os.path.join(directory_cld,map_cld['ECOv002_L2_CLOUD_'+key_eco])
     
    stic_files = os.listdir(stic_path)  # liste des fichiers dans le dossier
    stic_file = next((f for f in stic_files if key_eco in f), None)

    stic_file = stic_path+'/'+stic_file

    if stic_file:
        wasdi.wasdiLog(f"Found STIC file: {stic_file}")
    else:
        wasdi.wasdiLog(f"Aucun fichier trouvé contenant {key_eco}")
    #stic_file = [f for f in stic_path if key_eco in f][0]

    # 1.2 Early exit 
    # exit if required files are missing
    if not all([path_geo, path_cld, stic_file]):
        wasdi.wasdiLog("Corresponding data not available")
        return
    # exit if cloudy image(>75%) to save time             
    try:
        start_time = time.time()
        mask_cld = Create_cloud_mask(path_cld)
        wasdi.wasdiLog(f"Reading Cloud file done in {time.time() - start_time:.2f} seconds")
        cld_percent = 1-np.mean(mask_cld)
        if cld_percent >= 0.75:
            wasdi.wasdiLog(f'too many cloud cover {int(cld_percent*100)}%, pass')
            return
    except Exception as e:
        wasdi.wasdiLog(f"Error processing Cloud file: {e}")
        return
    #test_file_GPP(writepath,mask_cld,key_eco)



    # 1.3 get lat, lon, SZA from and ECOSTRESS metadata
    eco_lat,eco_lon,eco_sza = read_geo(path_geo)
    # exit if nighttime: SZA>90 done in first round

    lat_max, lat_min, lon_max, lon_min = np.nanmax(eco_lat),np.nanmin(eco_lat),np.nanmax(eco_lon),np.nanmin(eco_lon)
    eco_bound = (lat_max+0.5, lat_min-0.5, lon_max+0.5, lon_min-0.5) #buffering with 0.5 degree for interpolation
    # wasdi.wasdiLog(eco_bound)

    # 1.3 STIC-gs and ET from EEH2 product    
    try:
        start_time = time.time()
        (ETD, gt) = Extract_STIC(stic_file, mask_cld)
        valid_percent = np.mean(ETD>0)
        wasdi.wasdiLog(valid_percent)
        if valid_percent==0:
            wasdi.wasdiLog("There is no valid ETD>0, GPP output will be full of NaN")
            return

        wasdi.wasdiLog(f"{stic_file} \nExtracting STIC results done in {time.time() - start_time:.2f} seconds")
    except Exception as e:
        wasdi.wasdiLog(f"Error extracting STIC result: {e}")
        return

    #####################################################################################
    # STEP2 read and rescale other inputs
    #####################################################################################
    # 2.1 read other data into xarray and projections by date

    # Read global scale daily/annual dataset:  LAI, FVC, LUCC, OCO2,PAR,PARH
    # (1) PAR daily mean & PARH: 
    path_PARH = [p for p in files_PARH if  f"PARin{year_str}{month_str}{day_str}" in p][0] 

    # find the cloest half-hourly PAR observartions, calculate PARsum based on half-hourly PAR    
    PARH, PARmean = extract_PAR_from_global_data(path_PARH, eco_bound, hour,minute)
    # match ECOSTRESS tiles
    PARH_70m = wrap_to_ECOSTRESS(PARH, eco_lat, eco_lon) 
    PARDmean_70m = wrap_to_ECOSTRESS(PARmean, eco_lat, eco_lon) # in  µmolPhoton m-2 s-1
    PAR_scaler_raw = PARDmean_70m/PARH_70m
    PAR_scaler = xr.where((PAR_scaler_raw < 100) & (PAR_scaler_raw > 1/100), PAR_scaler_raw, np.nan)

    wasdi.wasdiLog("PARH done")

    # (2) Read FCOVER and interpolation by time 
    # RT0: Near Real-Time (NRT) data, often the most recent but less consolidated.
    # RT1 to RT6: Progressively more refined or consolidated versions, 
    # with RT6 being the most accurate and validated, incorporating corrections and additional data.
    FVC_300m = interpolate_FVC_value(year, month, day, files_fvc, eco_bound)
    FVC_70m = wrap_to_ECOSTRESS(FVC_300m, eco_lat, eco_lon)

    wasdi.wasdiLog("FVC done")

    #  (3) Read ERA5 and interpolation by time 
    # read temperature and calculate vpd from ERA5 single level data
    files_era5_single = [f for f in files_era5 if 'single' in f]
    t2m_ERA5_10km,vpd_ERA5_10km = interpolate_ERA5_value(year, month, day, hour, minute,files_era5_single, eco_bound)
    temp_70m = wrap_to_ECOSTRESS(t2m_ERA5_10km, eco_lat, eco_lon)
    vpd_70m = wrap_to_ECOSTRESS(vpd_ERA5_10km, eco_lat, eco_lon)

    wasdi.wasdiLog("ERA5 done")


    # (4) OCO2: match with date 
    oco2_50km = read_oco2(year,month,day,eco_bound, files_oco2)
    oco2_70m = wrap_to_ECOSTRESS(oco2_50km, eco_lat, eco_lon)

    wasdi.wasdiLog("OCO2 done")

    # (5) global GLC30: Read regional tiles and mosaic (done in first round)

    # 2.2 read outputs from last edition of ECOSTRESS L3 product
    (GPP_BL_raw,LUCC) = Extract_ECOSTRESS_from_mask(gpp_file, mask_cld,['GPP_BL','LUCC'])
    GPP_LUE = xr.where((GPP_BL_raw!=-9999), GPP_BL_raw, np.nan)


    # 2.3 stacking inputs for modelling
    inputs = {
        "ETD": ETD,#in mm/day
        'gt': gt,#in m/s
        "FVC": FVC_70m.data,
        'Tair': temp_70m.data, #in °C 
        'VPD': vpd_70m.data, #in hPa
        "PAR_scaler": PAR_scaler.data,#µmolPhoton m−2 s−1
        "CO2": oco2_70m.data,#in mol CO2/mol dry air
    }
   

    #####################################################################################
    # STEP3 calculate hourly GPP and daily GPP&WUE
    #####################################################################################
    # 3.1 hourly GPP: LUE-gs 

    mol_to_umol = 1e6
    umolCO2_to_gC =  12 * 1e-6 # Convert from µmolCO2 to gC
    mp_to_molm2 = 40.088/1.6
    GPPumol_max = 60
    s_to_day = 60 *60* 24 #second to day
    

    GPP_gs =gpp_gs_calculator(inputs['gt'],inputs['CO2'],  inputs['Tair'], inputs['VPD']).clip(0,0.001) # gC m−2  s−1  
    

    # hybrid model
    vpd_mask = (vpd_70m>9) #hPa
    # Element-wise minimum of two arrays, return NaN if  either element is nan
    GPPmin = np.where(vpd_mask,np.minimum(GPP_LUE, GPP_gs),GPP_LUE).clip(0,0.001)  # handling nan, gC m−2 s−1

    # final nan mask
    veg_p25_mask = (inputs['FVC']>0.25)
    eco_valid_mask = (ETD > 0.01) & (~np.isnan(ETD))&veg_p25_mask #make sure ETD is not close to 0
    GPP_LUE,GPP_gs,GPPmin,ETD_safe = np.where(eco_valid_mask,GPP_LUE,np.nan),\
    np.where(eco_valid_mask,GPP_gs,np.nan),\
    np.where(eco_valid_mask,GPPmin,np.nan),\
    np.where(eco_valid_mask,ETD,np.nan)


    # 3.2 Daily GPPd based on PARh/PARd scalar 
    GPPd_gC = ((GPPmin*PAR_scaler_raw)*s_to_day)#.clip(0,30)  # gC m−2  day−1 

    # 3.3 WUE=GPPd/ETd
    # where GPPd  is in gC m-2day-1, ET is in mm day-1 and the WUEd is in gC m-2/mm H2O
    # Compute WUE safely
    WUE_raw = (GPPd_gC / ETD_safe) #gC m-2/mm H2O
    # WUE = xr.where((WUE_raw > 0.05) & (WUE_raw < 10), WUE_raw, np.nan)

    # 3.4 Output: GPPd,  WUEd, landcover type 
    outputs = {
        'GPPd': np.nan_to_num(GPPd_gC, nan=-9999),
        'WUEd': np.nan_to_num(WUE_raw, nan=-9999),
        'GPP_gs': np.nan_to_num(GPP_gs, nan=-9999),
        'GPP_LUE': np.nan_to_num(GPP_LUE, nan=-9999),
        'GPPmin': np.nan_to_num(GPPmin, nan=-9999),
        'LUCC': np.nan_to_num(LUCC, nan=-9999)
    }

   
    # 3.5 set attributes in output files
    generate_hdf5_file_GPP(writepath,outputs,key_eco) 



def run():
    # matching other data with ECOSTRESS 
    # df_bound = pd.read_csv('ECOSTRESS_orbit_bounds_updated.csv')
    df_LUT = pd.read_csv(r'LookUpTable_LUE_HH_gsFULL-globe.csv')

    # Check if 'nonV' is missing and add it
    if 'nonV' not in df_LUT.index:
        # Create a new row with same columns, filled with 0
        new_row = pd.DataFrame([0] * len(df_LUT.columns), 
                               index=df_LUT.columns).T
        new_row['type'] = 'nonV'
        # Append to the DataFrame
        df_LUT = pd.concat([df_LUT, new_row])

    #####################################################################################
    # STEP1 read the meta and geo info of ECOSTRESS observation
    #####################################################################################
    # 1.1 read meta and geo from Cloud

    rp         = wasdi.getSavePath()
    wp         = rp+'ecostress_fast/'

    # GPP file names from first round GPP product (all files under folder ./EEHGPP)
    gpp_file  = wasdi.getPath(wasdi.getParameter('GPP_FILE'))
    # GPP file names from unique_GPP_path_for_ICOS_evaluation.txt
    #gpp_file_for_evaluation  = wasdi.getPath(wasdi.getParameter('GPP_FILE_for_evaluation'))
    output_path= wasdi.getPath(wasdi.getParameter('GPP_OUTPUT_PATH'))

    stic_path   = rp+wasdi.getParameter('STIC_PATH')
    fvc_path   = rp+wasdi.getParameter('FVC_PATH')
    geo_path   = rp+wasdi.getParameter('GEO_PATH')
    cld_path   = rp+wasdi.getParameter('CLD_PATH')
    era5_path  = rp+wasdi.getParameter('ERA5_PATH')
    parh_path  = rp+wasdi.getParameter('PARH_PATH')
    oco2_path  = rp+wasdi.getParameter('OCO2_PATH')

    #map of STIC/CLOUD and GEO
    map_cld   = cache_S3_with_pattern('CLOUD',wp,r'(ECOv002_L2_CLOUD.*)_\d{4}_\d{2}.h5')
    map_geo   = cache_S3_with_pattern('L1B_GEO',wp,r'(ECOv002_L1B_GEO.*)_\d{4}_\d{2}.h5')

    #just list of files
    files_parh = cache_S3_simple('PARH',wp,parh_path)
    files_fvc  = cache_S3_simple('FCOVER',wp,fvc_path)
    files_era5 = cache_S3_simple('ERA5',wp,era5_path)
    files_oco2 = cache_S3_simple('OCO2',wp,oco2_path)

    wasdi.updateStatus("RUNNING", 0)
    refine_GPP(gpp_file,output_path,map_cld,map_geo,stic_path,files_parh,files_fvc,files_era5,files_oco2,cld_path,geo_path)

    wasdi.updateStatus("DONE", 100)

    
if __name__ == '__main__':
    wasdi.init("./config.json")
    run()

