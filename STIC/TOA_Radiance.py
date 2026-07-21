#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Script to run the STIC model for ET estimation

Created on September 1 2024
@author: Tian Hu at LIST

© 2022 – Luxembourg Institute of Science and Technology
Authors : Tian Hu, Kaniska Mallick
Code licensed under MIT
SPDX-License-Identifier: MIT
"""

#Function for calculating TOA solar radiation
import numpy as np

#Calculate solar zenith angle
def f_SolarZenith(doy, time_decimal, latitude_deg, longitude_deg):
    #Convert degrees to radians
    latitude = np.radians(latitude_deg)
    
    #Calculate day angle
    gamma = 2 * np.pi * (doy - 1) / 365
    
    #Calculate solar declination
    delta = (0.006918 
         - 0.399912 * np.cos(gamma) 
         + 0.070257 * np.sin(gamma)
         - 0.006758 * np.cos(2 * gamma)
         + 0.000907 * np.sin(2 * gamma)
         - 0.002697 * np.cos(3 * gamma)
         + 0.00148  * np.sin(3 * gamma))
    
    #Calculate equation of time
    eot = 229.18 * (0.000075 + 0.001868 * np.cos(gamma)
                    - 0.032077 * np.sin(gamma)
                    - 0.014615 * np.cos(2 * gamma)
                    - 0.040849 * np.sin(2 * gamma))
    
    #Estimate time zone from longitude
    #timezone_offset = round(longitude_deg / 15)
    timezone_offset = np.round(longitude_deg / 15)
    
    standard_meridian = 15 * timezone_offset
    
    #Time correction
    tc = 4 * (longitude_deg - standard_meridian) + eot
    
    #Local solar time
    solar_time = time_decimal + tc / 60.0
    
    #Hour angle
    hra = np.radians(15 * (solar_time - 12))
    
    #Solar zenith angle
    cos_theta = (np.sin(latitude) * np.sin(delta)
             + np.cos(latitude) * np.cos(delta) * np.cos(hra))
    
    #cos_theta = min(max(cos_theta, -1), 1)  # Clamp to [-1, 1]
    cos_theta = np.clip(cos_theta, -1, 1)

    zenith_angle_rad = np.arccos(cos_theta)
    zenith_angle_deg = np.degrees(zenith_angle_rad)
    
    return zenith_angle_deg

#Calculate instantaneous TOA radiance
def f_RgTOAi(doy, sza_deg):
    RgTOAi = 1367 * (1 + 0.033 * np.cos(2 * np.pi * doy / 365)) * np.sin(np.radians(90 - sza_deg)) 
    RgTOAiMJ = RgTOAi * 60 * 30 * 1e-6; #Exoatmospheric shortwave radiation instantaneous (MJ)
    
    return RgTOAiMJ

#Calculate daily TOA radiance
def f_enot(doy, latitude_deg):
    day_ang = (2 * np.pi * (doy - 1)) / 365
    
    E0 = 1.00011 + 0.034221 * np.cos(day_ang) + 0.00128 * np.sin(day_ang) + \
          0.000719 * np.cos(2 * day_ang) + 0.000077 * np.sin(2 * day_ang)
    
    soldecl = (0.006918 - 0.399912 * np.cos(day_ang) + 0.070257 * np.sin(day_ang) - \
          0.006758 * np.cos(2 * day_ang) + 0.000907 * np.sin(2 * day_ang) - \
          0.002697 * np.cos(3 * day_ang) + 0.00148 * np.sin(3 * day_ang)) * 180 / np.pi
    
    sha_deg = np.arccos(-np.tan(np.pi * latitude_deg / 180) * np.tan(np.pi * soldecl / 180)) * (180 / np.pi)
    
    #Exoatmospheric shortwave radiation integrated (MJ)
    RexoInt = (24 / np.pi) * 1367 * 0.0036 * E0 * np.cos(latitude_deg * np.pi / 180) * np.cos(soldecl * np.pi / 180) * \
              (np.sin(sha_deg * np.pi / 180) - (sha_deg * np.pi / 180) * np.cos(sha_deg * np.pi / 180)) 
    
    return RexoInt

def f_TOARadiance(doy, time_decimal, latitude_deg, longitude_deg):
    #Calculate SZA
    sza_deg = f_SolarZenith(doy, time_decimal, latitude_deg, longitude_deg)
    
    #Calculate instantaneous TOA radiance
    RgTOAiMJ = f_RgTOAi(doy, sza_deg)
    
    #Calculate daily integrated TOA radiance
    RexoInt = f_enot(doy, latitude_deg)
    
    return (RgTOAiMJ, RexoInt)
    
    
    
    
    
    

