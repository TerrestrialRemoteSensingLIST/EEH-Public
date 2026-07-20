#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Script to run the Temperature Emissivity Separation (TES) algorithm for LST estimation

Created on April 15 2021
@author: Tian Hu at LIST

© 2022 – Luxembourg Institute of Science and Technology
Authors : Tian Hu, tian.hu@list.lu
Code licensed under MIT
SPDX-License-Identifier: MIT
"""

# Atmospheric correction
import pyrttov
import sys
import numpy as np

rttov_installdir = '/root/rttov/' #preinstalled RTTOV lib

# Transform 1D array (same for all profiles) to a [nprof, nlevels] array
def expand2nprofiles(n, nprof):
    outp = np.ones((nprof, len(n)))*n
    return outp

def runRTTOV(p_era5, new_t, new_q, new_sp, new_q2m, new_t2m, new_skt, vza, sza, watermask,
             lat_eco, lon_eco, alt_eco, year_str, month_str, day_str):
 
    nlevels = new_t.shape[0] #level number of atmospheric profiles
    nprofiles = new_t.shape[1] #profile number for the entire ECOSTRESS subimage
    nchans = 3 #channel numer for ECOSTRESS TIR bands
    chan_list = (2,4,5)
    myProfiles = pyrttov.Profiles(nprofiles, nlevels)
    
    #Associate the atmospheric profile data with myProfiles
    myProfiles.GasUnits = 1 #2: ppmv over moist air (default)    1: kg/kg over moist air     0 (or less): ppmv over dry air

    new_p = expand2nprofiles(p_era5, nprofiles)
    myProfiles.P = np.flip(new_p,axis=1) #pressure levels (hPa)

    myProfiles.T = np.flip(new_t.transpose(),axis=1) #atmospheric temperature (K)

    myProfiles.Q = np.flip(new_q.transpose(),axis=1) #water vapour (unit as per GasUnits)

    #[p,t,q,u,v,wfetc]
    s2m = np.zeros((nprofiles,6))
    s2m[:,0] = new_sp
    s2m[:,1] = new_t2m
    s2m[:,2] = new_q2m
    myProfiles.S2m = s2m

    #[t,salinity,snowfraction,foamfraction,fastem(1:5),specularity]
    skin = np.zeros((nprofiles,10))
    #skin = np.zeros((nprofiles,9))
    skin[:,0] = new_skt
    print(skin.shape)
    myProfiles.Skin = skin

    #[satzenith,satazimuth,sunzenith,sunazimuth]
    angle = np.zeros((nprofiles,4))
    angle[:,0] = vza.reshape(vza.shape[0]*vza.shape[1])
    angle[:,2] = sza.reshape(sza.shape[0]*sza.shape[1])
    myProfiles.Angles = angle

    #[surfacetype,watertype]
    surftype = np.zeros((nprofiles,2))
    surftype[:,0] = watermask.reshape(watermask.shape[0]*watermask.shape[1])
    myProfiles.SurfType = surftype

    #[latitute(deg),longitude(deg),elevation(km)]
    surfgeom = np.zeros((nprofiles,3))
    surfgeom[:,0] = lat_eco.reshape(lat_eco.shape[0]*lat_eco.shape[1]) #[-90,90]
    surfgeom[:,1] = lon_eco.reshape(lon_eco.shape[0]*lon_eco.shape[1]) + 180 #[0,360]
    surfgeom[:,2] = alt_eco.reshape(alt_eco.shape[0]*alt_eco.shape[1])/1000 #km
    myProfiles.SurfGeom = surfgeom

    #[year,month,day,hour,minute,second], not currently used by RTTOV, only month information are used for emissivity selection
    time = np.array([int(year_str),int(month_str),int(day_str),0,0,0],dtype=np.int32)
    myProfiles.DateTimes = expand2nprofiles(time, nprofiles) 

    print('Step 1 finished!')

    #######################################################
    #2. Set up Rttov instance for the ECOSTRESS instrument
    #######################################################
    Rttov = pyrttov.Rttov()
    Rttov.FileCoef = '{}/{}'.format(rttov_installdir,"rtcoef_rttov12/rttov8pred54L/rtcoef_iss_1_ecostres.dat")
    #Rttov.FileCoef = '{}/{}'.format(rttov_installdir,"rtcoef_rttov13/rttov8pred54L/rtcoef_iss_1_ecostres.dat")
    Rttov.Options.AddInterp = True
    Rttov.Options.StoreTrans = True
    # Rttov.Options.StoreRad = True
    Rttov.Options.StoreRad2 = True
    Rttov.Options.VerboseWrapper = True
    Rttov.Options.ApplyRegLimits = True
    Rttov.Options.Nthreads = 4
    Rttov.Options.NprofsPerCall = 15

    try:
        Rttov.loadInst(chan_list)
    except pyrttov.RttovError as e:
        sys.stderr.write("Error loading instrument(s): {!s}".format(e))
        sys.exit(1)

    Rttov.Profiles = myProfiles

    print('Step 2 finished!')

    ##########################################
    #3. Load the emissivity and BRDF atlases
    ##########################################

    # Set up the surface emissivity/reflectance arrays and associate with the Rttov objects
    surfemisrefl = np.zeros((2,nprofiles,nchans), dtype=np.float64)
    Rttov.SurfEmisRefl = surfemisrefl

    print('Step 3 finished!')

    ############################################################################################
    #4. Run RTTOV at specific VZAs to retrieve transmittance and atmospheric upwelling radiance
    ############################################################################################
    try:
        Rttov.runDirect()
    except pyrttov.RttovError as e:
        sys.stderr.write("Error running RTTOV direct model: {!s}".format(e))

    print('Step 4 finished!')

    trans = Rttov.TauTotal  #transmittance, output atmospheric parameter 1

    upclear = Rttov.Rad2UpClear  #output atmospheric parameter 2

    #converse all the radiances from mW/cm-1/sr/m2 to W/um/sr/m2
    #wn1 = 0.1207322786E+04
    wn2 = 0.1138166913E+04
    #wn3 = 0.1087276567E+04
    wn4 = 0.9534730659E+03
    wn5 = 0.8272245003E+03
    wn = np.array([[wn2,wn4,wn5]])

    upclear1 = upclear*wn**2*1E-7

    ########################################################################################################
    #5. Run RTTOV at VZA of 53 deg to retrieve hemispherically effective atmospheric downwelling irradiance
    ########################################################################################################
    angle[:,0] = expand2nprofiles(np.array([53.0]), nprofiles)[:,0]

    try:
        Rttov.runDirect()
    except pyrttov.RttovError as e:
        sys.stderr.write("Error running RTTOV direct model: {!s}".format(e))

    print('Step 5 finished!')

    #atmospheric downwelling irradiance
    dnclear = Rttov.Rad2DnClear  #output atmospheric parameter 3

    #converse all the radiances from mW/cm-1/sr/m2 to W/um/sr/m2
    dnclear1 = dnclear*wn**2*1E-7
    
    return (trans, upclear1, dnclear1)
