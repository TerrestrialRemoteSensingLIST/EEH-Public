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

# Component function of STIC

import numpy as np
from Radiation import *
from SMWetness import *

# Constants
sigma = 5.67e-8     # Stefan-Boltzmann constant
gamma = 0.67        # Psychrometric constant (hpa/K)
Rd = 286.9          # Gas constant for dry air
Rw = 461.5          # Gas constant for wet air
P = 101325.         # atmosphere pressure (Pascal)
Cpw = 1846.         # Specific heat of moist air
Cpd = 1005.         # Specific heat of air at constant pressure (J/kg/K)
threshold = 0.01    # Threshold for termination
alfa0 = 1.26         # Priestley-Taylor coefficient
kPAR = 0.5
kRN = 0.6   

def f_pSYCHROMETRICS(TS,TA,RH,P,Rd,Rw,Cpw,Cpd,sigma):  
    esstar = 6.13753*np.exp((17.27*TS)/(TS+237.3)) # saturation vapor pressure at surface temperature, TS (unit hPa)  
    eastar = 6.13753*np.exp((17.27*TA)/(TA+237.3)) # saturation vapor pressure at air temperature, TA (unit hPa)  
    ea = (RH/100)*(eastar) # actual vapor pressure of air (unit hPa)

    DA = eastar - ea # vapor pressure deficit of air (hPa) 
    slopeTA = 4098*eastar/(TA + 237.3)**2 # slope of saturation vapor pressure versus air temperature (hPa)       
    TD = 237.3*np.log(ea/6.13753)/(17.27 - np.log(ea/6.13753)) # dewpoint temperature (deg C)
    TD[TD < -30] = np.nan
    DPD = TA - TD # dewpoint depression
    
    slopeTD = 4098*ea/(TD + 237.3)**2
    slopeTaTd = (eastar - ea)/(TA - TD) # slope of saturation vapor pressure at dewpoint temperature
    
    qref = (0.622*ea*100)/(101325-(0.387*ea*100))
    rhoD = P/(Rd*(TA+273.15)) # density of dry air
    rho = rhoD*((1+qref)/(1+qref*(Rw/Rd))) # density of air
    CP = qref*Cpw + (1-qref)*Cpd # specific heat of air
    
    gR = ((4*sigma*(TA+273.15)**3)/CP)*((Rd*(TA+273.15))/P) # radiative conductance (m/s)(momentarily not needed)

    dTS = TS - TA # difference between TS and TA
    
    return (esstar,eastar,ea,DA,slopeTA,TD,DPD,slopeTD,slopeTaTd,rho,CP,gR,dTS)

def f_StateEQ(rho,CP,gamma,alfa,slope,PHI,e0,ea,e0star,M):
    # Aerodynamic conductance
    gB = (2* PHI* alfa* slope* gamma)/ (2* CP* slope* e0* rho -
                             2* CP* slope* ea* rho - 2* CP* ea* gamma* rho +
                             CP* e0* gamma* rho + CP* e0star* gamma* rho -
                             CP* M* e0* gamma* rho + CP* M* e0star* gamma* rho)
    
    # Surface conductance
    denominator = (CP* e0star** 2* gamma* rho - CP* e0** 2* gamma* rho -
        2* CP* slope* e0 ** 2* rho + 2* CP* slope* ea* e0* rho -
        2* CP* slope* ea* e0star* rho + 2* CP* slope* e0* e0star* rho +
        2* CP* ea* e0* gamma* rho - 2* CP* ea* e0star* gamma* rho +
        CP* M* e0 ** 2* gamma* rho + CP* M* e0star ** 2* gamma* rho -
        2* CP* M* e0* e0star* gamma* rho)
    denominator[denominator == 0] = 0.00001  #when e0star == e0
    gS = -(2* (PHI* alfa* slope* ea* gamma - PHI* alfa* slope* e0* gamma))/ denominator
    
    # T0 - TA
    dT = (2* slope* e0 - 2* slope* ea - 2* ea* gamma + e0* gamma +
          e0star* gamma - M* e0* gamma + M* e0star* gamma +
          2* alfa* slope* ea - 2* alfa* slope* e0)/ (2* alfa* slope* gamma)
    
    # Evaporative fraction
    EF = -(2* alfa* slope* ea - 2* alfa* slope* e0)/ (2* slope* e0 -
                                      2* slope* ea - 2* ea* gamma + e0* gamma + e0star* gamma -
                                      M* e0* gamma + M* e0star* gamma)
    
    # Adjust the abnormal conductances
    gB[gB < 0] = .0001  
    gB[gB > 0.2] = 0.2
    gS[gS < 0] = .0001  
    gS[gS > 0.06] = 0.06

    # Maximum surface-air temperature difference rarely overpasses 20 degC
    dT[dT < -10] = -10    
    dT[dT > 20] = 20
    EF[EF < 0] = .0001       
    EF[EF > 1] = 1

    return (gB,gS,dT,EF)

def STIC(TS,etas,TA,RH,sr_dir,sr_dif,alb_dir,alb_hem,fc,ttSEC):       
    # Mask for valid pixels
    mask = np.ones(TS.shape)
    mask[np.logical_or.reduce((np.isnan(TS),np.isnan(etas),np.isnan(alb_dir),np.isnan(alb_hem),np.isnan(fc)))] = 0
    
    # Calculate the total solar insolation (W.m-2)
    RG = sr_dir + sr_dif  
    
    ## Psychrometrics
    # esstar, saturated vapour pressure at surface temperature (hPa)
    # eastar, saturated vapour pressure at air temperature (hPa)
    # ea, atmosphere vapour pressure (hPa)
    # DA, atmosphere vapour pressure deficit (hPa)
    # slopeTA, slope of saturation vapor pressure versus air temperature at TA (hPa/degC)
    # TD, dewpoint temperature (degC)
    # DPD, dewpoint depression, difference between TA and TD (used for cross verification of evaporation, degC)
    # slopeTD, slope of saturation vapor pressure versus temperature at TD (hPa/degC)
    # slopeTaTd, not used in the model (used for checking different slopes, hPa/degC)
    # rho, air density (kg.m-3)
    # CP, specific heat of air at constant pressure (MJ.kg-1.K-1) 
    # dTS, TS - TA (degC)
    (esstar,eastar,ea,DA,slopeTA,TD,_,slopeTD,slopeTaTd,rho,CP,gR,dTS) = f_pSYCHROMETRICS(TS,TA,RH,P,Rd,Rw,Cpw,Cpd,sigma)
    
    ## Hypothetic wet conditions used for calculating potential evapotranspiration
    # Twet, wet surface temperature (degC)
    # eWETstar, saturated vapour pressure at Twet (hPa)
    # sWET, slope of saturation vapor pressure versus temperature at Twet (hPa/degC) 
    # taWB, wet bulb temperature (degC)
    # tdWB, dewpoint wet bulb temperature (degC)
    # eastarWB, saturation vapor pressure at taWB (hPa)
    # slopeWB, slope of saturation vapor pressure versus air temperature at taWB (hPa/degC) 
    # eaWB, vapour pressure at taWB (degC)
    # dTwet, Twet - taWB (degC)
    # DAwet, vapour pressure deficit at taWB (hPa)       
    #(Twet,eWETstar,_,taWB,tdWB,eastarWB,slopeWB,eaWB,dTwet,DAwet) = f_Twet(eastar,slopeTA,TA,slopeTD,TD,gamma)

    #(_,Lout_wet,Lnet_wet,RN_wet) = f_Radiation_wetsurface(RG,taWB,eaWB,Twet,sigma)

    #(Mwet,_,MsoilWET,MsWET,MrzWET,s11WET,s22WET,s33WET,s44WET,
         #esWET,_,DsWET,TdewIndexWET) = f_SoilMoisture_INITIALIZE(gamma,slopeWB,Twet,taWB,tdWB,dTwet,                                                                                     RG,RN_wet,Lnet_wet,fc,DA,eastarWB,eaWB,eWETstar)

    #(_, _, _, _, PHI_wet) = f_G_PHI_wetsurface(RN_wet,ttSEC,MsoilWET) # PHI = RN - G
    
    #(gBwet,gSwet,dTwet,_) = f_StateEQ(rho,CP,gamma,alfa0,slopeWB,PHI_wet,esWET,eaWB,eWETstar,Mwet)

    #EpWET = (slopeWB*PHI_wet + rho*CP*gBwet*DAwet)/(slopeWB + gamma) # Potential evaporation (Penman)

    #TpWET = (slopeWB*RN_wet + rho*CP*gBwet*DAwet)/(slopeWB + gamma*(1 + MrzWET*gBwet/gSwet)) # Potential transpiration

    ## Actual conditions
    # etaa, air emissivity (0-1, unitless)
    # Lin, downwelling atmosphere longwave radiation (W.m-2)
    # Lout, surface upwelling longwave radiation (W.m-2)
    # Lnet, net longwave radiation (W.m-2)
    # RN, net radiation (W.m-2)
    (_,Lin,Lout,Lnet,Sin,Sout,Snet,RN) = f_NetRadiation(sigma,TA,ea,TS,etas,sr_dir,sr_dif,alb_dir,alb_hem)   
    
    # Initialize
    (M,Mcan,Msoil,Ms,Mrz,s11,s22,s33,s44,
             es,TSD,Ds,TdewIndex) = f_SoilMoisture_INITIALIZE(gamma,slopeTA,TS,TA,TD,dTS,RG,RN,Lnet,fc,DA,eastar,ea,esstar)
    
    T0Dold = TSD.copy()  
    e0star = esstar.copy()
    e0 = es.copy()
    slope = slopeTA.copy()

    LAI = -np.log( 1- fc)/kPAR #Correction on 12/11/2021
    RNsoil = RN*np.exp(-kRN*LAI) #Correction on 12/11/2021
    
    (_, _, _, G, PHI) = f_G_PHI_actualsurface(RN,RNsoil,ttSEC,Msoil) #Correction on 12/11/2021
    
    alfa = 1.26
    EpWET = (alfa*slope*PHI)/(slope + gamma)

    (gB,gS,dT,EF) = f_StateEQ(rho,CP,gamma,alfa,slope,PHI,e0,ea,e0star,M)
    
    gB[gB < 0] = .0001   
    gB[gB > 0.2] = 0.2
    gS[gS < 0] = .0001   
    gS[gS > 0.06] = 0.06

    dT[dT < -10] = -10     
    dT[dT > 20] = 20
    EF[EF < 0] = .0001      
    EF[EF > 1] = 1
    
    gBB = gB.copy()
    gSS = gS.copy()
    dTT = dT.copy()
    EFF = EF.copy()

    T0 = dTT + TA
    slope00 = (e0star - ea)/(T0 - TD) 
    #slope00 = slope.copy()
    
    # Calculate ET and H based on initial results from state eqs.
    # McNaughton and Jarvis (1986)
    omega1 = ((slope/gamma) + 1)/((slope/gamma) + 1 + gBB/gSS) 
    LHFeq1 = (PHI*(slope/gamma))/((slope/gamma) + 1)
    LHFimp1 = (CP*0.0289644/gamma)*gSS*40*DA
    LHF1 = omega1*LHFeq1 + (1 - omega1)*LHFimp1

    # Alternative form to initilize LHF, from Martin (1988, AFM)
    gR = (4*sigma*etas*(TA+273.15)**3)/(rho*CP) #Radiative conductance
    omega2 = ((slope/gamma) + 1 + gR/gBB)/((slope/gamma) + 1 + gBB/gSS + gR/gBB + gR/gSS)
    LHFeq2 = (PHI*(slope/gamma))/((slope/gamma) + 1 + gR/gBB)
    LHFimp2 = (CP*0.0289644/gamma)*gSS*40*DA
    LHF2 = omega2*LHFeq2 + (1 - omega2)*LHFimp2

    LHF = LHF1.copy()
    SHF = (gamma*PHI*(1 + gBB/gSS) - rho*CP*gBB*DA)/(slope + gamma*(1+gBB/gSS)) # Deduced from the PM equation

    Ep = (slope*PHI + rho*CP*gBB*DA)/(slope + gamma) # Potential evaporation (Penman)
    Tp = (slope*RN + rho*CP*gBB*DA)/(slope + gamma*(1 + M*gBB/gSS)) # Potential transpiration            

    Es = (1-fc)*Ms*EpWET # Soil evaporation
    Ei = Mcan*EpWET # Canopy/interception evaporation
    Ei[(Ep > PHI) & (Ds-DA > 0) & (dTS > 0)] = 0

    EVAP = Es + Ei # Total evaporation                                     
    TRANSP = LHF - EVAP # Transpiration

    #rDIFF = PHI/(SHF + LHF)
    
    # Surface dryness constraint, correct the initialization
    #index = np.logical_and.reduce((PHI > 0,rDIFF > 1,Ds-DA > 0,TD < 0,dTS > 0))
    #LHF[index] = (1-omega2[index])*LHFimp2[index]
    #SHF[index] = PHI[index] - LHF[index]
    
    # Surface dryness constraint, correct the initialization
    #index = np.logical_and.reduce((PHI > 0,rDIFF < 1,Ds-DA > 0,TdewIndex > 1))
    #LHF[index] = (1-omega2[index])*LHFimp2[index]
    #SHF[index] = PHI[index] - LHF[index]
    
    # Wet or water unlimited condition
    #index = np.logical_and.reduce((PHI > 0,rDIFF < 1,Ds-DA < 0))
    #LHF[index] = omega2[index]*LHFeq2[index]
    #SHF[index] = PHI[index] - LHF[index]
    
    # Wet or water unlimited condition
    #index = np.logical_and.reduce((PHI > 0,rDIFF > 1,Ds-DA < 0,dTS < 0))
    #LHF[index] = omega2[index]*LHFeq2[index]
    #SHF[index] = PHI[index] - LHF[index]
    
    # Wet or water unlimited condition
    #index = np.logical_and(PHI > 0,Ds-DA < 0)
    #LHF[index] = omega2[index]*LHFeq2[index]
    #SHF[index] = PHI[index] - LHF[index]
    
    # Re-constraining conductance based on the corrected LHF and SHF initilization
    #EFF = LHF/(SHF + LHF)
    #EFF[EFF < 0] = 0.01
    #EFF[EFF > 1] = 0.99
    
    #dTT = ((1-EFF)/EFF)*((e0 - ea)/gamma)
    #dTT[dTT < -10] = -10    
    #dTT[dTT > 20] = 20
    #T0 = dTT + TA

    #gBB = SHF/(rho*CP*dTT)
    #gBB[gBB < 0] = .0001
    #gBB[gBB > 0.1] = 0.1

    #gSS = (gamma*gBB*LHF)/(slope*PHI + rho*CP*gBB*DA - LHF*(slope + gamma)) # Deducted from the PM equation
    #gSS[gSS < 0] = .0001   
    #gSS[gSS > 0.1] = 0.1
    
    # Iteration
    LHFold = LHF.copy()
    LHFnew = LHF.copy()
    error_LHF = np.ones(TS.shape)*0.05
    steps = 0
    converged = np.asarray(np.zeros(TS.shape),dtype=bool) # Create the status array
    
    # Create the variable arrays
    [D0,T0D,alfaN] = [np.zeros(TS.shape,np.float32)+np.nan for i in range(3)]
    
    while(np.nanmax(error_LHF) > threshold and steps < 15):    
        index = np.where(converged)
        print('converged pixels: ',index[0].shape[0])
        
        # Selecting valid and non-converged pixels
        i = np.logical_and(mask,~converged)
        
        # Re-estimate saturated vapor pressure at source/sink height
        e0star[i] = ea[i] + (gamma*LHFnew[i]*(gBB[i] + gSS[i]))/(rho[i]*CP[i]*gBB[i]*gSS[i])       
        e0star[np.logical_and(i,e0star < 0)] = esstar[np.logical_and(i,e0star < 0)]
        e0star[np.logical_and(i,e0star > 250)] = 250
        
        # Re-estimate vapor pressure at source/sink height
        D0[i] = (gBB[i]/gSS[i])*(gamma/(slope[i] + gamma*(1 + gBB[i]/gSS[i])))*(DA[i] + ((slope[i]*PHI[i])/(rho[i]*CP[i]*gBB[i])))
        D0[np.logical_and(i,D0 < 0)] = Ds[np.logical_and(i,D0 < 0)]
        
        e0[i] = e0star[i] - D0[i]        
        e0[np.logical_and(i,e0 < 0)] = es[np.logical_and(i,e0 < 0)]
        e0[np.logical_and(i,e0 < ea)] = es[np.logical_and(i,e0 < ea)]
        e0[np.logical_and(i,e0 > e0star)] = es[np.logical_and(i,e0 > e0star)]
        
        # Re-estimate M (direct LST feedback into M computation)     
        T0D[i] = TD[i] + (gamma*LHFnew[i])/(rho[i]*CP[i]*gBB[i]*s11[i])    
        T0D[np.logical_and(i,T0D < TD)] = TD[np.logical_and(i,T0D < TD)]
        T0D[np.logical_and(i,T0D > TS)] = T0Dold[np.logical_and(i,T0D > TS)]
        
        (M[i],Ms[i],Mcan[i],Msoil[i],Mrz[i]) = f_SoilMoisture_ITERATE(gamma,slope[i],s11[i],s22[i],s33[i],s44[i],TS[i],TA[i],dTS[i],
                                         TD[i],T0D[i],RG[i],RN[i],Lnet[i],fc[i],DA[i],D0[i],eastar[i],ea[i],e0star[i],esstar[i])
        
        # Re-estimate PT coefficient
        alfaN[i] = (gSS[i]*slope00[i]*(T0[i] - TD[i])*(2*slope[i] + 2*gamma + gamma*(gBB[i]/gSS[i])*(1 + M[i])))/ \
                (2*slope[i]*(gamma*(T0[i] - TA[i])*(gBB[i] + gSS[i]) + gSS[i]*slope00[i]*(T0[i] - TD[i])))
        alfaN[np.logical_and(i,alfaN < 0)] = 1
        alfaN[np.logical_and(i,alfaN > 2)] = 2
        
        # Re-estimate net available energy
        (_, _, _, G[i], PHI[i]) = f_G_PHI_actualsurface(RN[i],RNsoil[i],ttSEC[i],Ms[i]) #Correction on 15/11/2021
        
        # Re-estimate conductances and states
        (gBB[i],gSS[i],dTT[i],EFF[i]) = f_StateEQ(rho[i],CP[i],gamma,alfaN[i],slope[i],PHI[i],e0[i],ea[i],e0star[i],M[i])
        
        gBB[np.logical_and(i,gBB < 0)] = 0.0001
        gBB[np.logical_and(i,gBB > 0.2)] = 0.2
        gSS[np.logical_and(i,gSS < 0)] = 0.0001
        gSS[np.logical_and(i,gSS > 0.06)] = 0.06
        dTT[np.logical_and(i,dTT < -10)] = -10
        dTT[np.logical_and(i,dTT > 20)] = 20
        EFF[np.logical_and(i,EFF < 0)] = 0.0001
        EFF[np.logical_and(i,EFF > 1)] = 1
        
        T0[i] = dTT[i] + TA[i]
    
        # Potential evaporation and potential transpiration 
        Ep[i] = (slope[i]*PHI[i] + rho[i]*CP[i]*gBB[i]*DA[i])/(slope[i] + gamma) # Potential evaporation (Penman)

        Ep[np.logical_and(i,Ep > RG)] = RG[np.logical_and(i,Ep > RG)]   

        Tp[i] = (slope[i]*RN[i] + rho[i]*CP[i]*gBB[i]*DA[i])/(slope[i] + gamma*(1 + M[i]*gBB[i]/gSS[i])) # Potential transpiration

        Tp[np.logical_and(i,Tp > RG)] = RG[np.logical_and(i,Tp > RG)]

        # Re-estimate LHF and SHF 
        LHFnew[i] = (rho[i]*CP[i]/gamma)*((gBB[i]*gSS[i])/(gBB[i] + gSS[i]))*(slope[i]*(T0[i] - TA[i]) + DA[i]) #An alternative is to use the PM equation to calculate LHF

        index = np.logical_and.reduce((LHFnew < 0,LHFnew < PHI,i))
        LHFnew[index] = rho[index]*CP[index]*gBB[index]*(e0[index] - ea[index])/gamma

        SHF[i] = (gamma*PHI[i]*(1 + gBB[i]/gSS[i]) - rho[i]*CP[i]*gBB[i]*DA[i])/(slope[i] + 
                                         gamma*(1+(gBB[i]/gSS[i]))) # surface heat flux, SHF = PHI - LHF, PM equation is only used for calculating SHF in the current version of STIC
        
        slope00[i] = ((gamma*LHFnew[i])/(rho[i]*CP[i]*gSS[i]) + (e0[i] - ea[i]))/(T0[i] - TD[i])
        
        # Re-conduct ET partitioning
        Es[i] = (1 - fc[i])*Ms[i]*EpWET[i]
        Ei[i] = Mcan[i]*EpWET[i]
        Ei[np.logical_and.reduce((i,Ep > RN,D0-DA > 0,dTS > 0))] = 0
        EVAP[i] = Es[i] + Ei[i]
        TRANSP[i] = LHFnew[i] - EVAP[i]
        
        # Error
        error_LHF[i] = np.abs(LHFold[i] - LHFnew[i])
        LHFold[i] = LHFnew[i]
        steps = steps + 1
        
        converged[i] = error_LHF[i] < threshold
        
    LHF = LHFnew.copy()

    print('Number of iterations = ', steps)

    # Final output from the STIC model
    RN_STIC = RN
    G_STIC = G
    H_STIC = SHF
    LE_STIC = LHF
    gah_STIC = gBB
    gsc_STIC = gSS
    T0_STIC = T0
    Ms_STIC = Ms
    Mrz_STIC = Mrz
        
    return (RN_STIC,G_STIC,H_STIC,LE_STIC,gah_STIC,gsc_STIC,T0_STIC,Ms_STIC,Mrz_STIC,converged,Lin,Lout,Lnet,Sin,Sout)
