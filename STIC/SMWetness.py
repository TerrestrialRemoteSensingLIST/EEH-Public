# -*- coding: utf-8 -*-

import numpy as np

def f_SoilMoisture_INITIALIZE(gamma,slope,TS,TA,TD,dTS,RG,RN,Lnet,fc,DA,eastar,ea,esstar):

    '''
        This function estimates the soil moisture availability 
        (M or wetness, 0-1) based on thermal and meteorological
        information. 
        However, this M will be treated as initial M, 
        which will be later on estimated through iteration in the actual ET estimation loop to
        establish feedback between M and biophysical states
    '''
    
    # TU computation (surface dewpoint temperature)
    s11 = (45.03 + 3.014*TD + 0.05345*TD**2 + 0.00224*TD**3)*1e-2 # slope of SVP at TD (hpa/K)
    s22 = (esstar - ea)/(TS - TD)
    s33 = (45.03 + 3.014*TS + 0.05345*TS**2 + 0.00224*TS**3)*1e-2 # slope of SVP at TS (hpa/K)
    s44 = (eastar - ea)/(TA - TD)                                                            

    # Surface dewpoint temperature (degC)
    T0D = (esstar - ea - s33*TS + s11*TD)/(s11 - s33) 

    # Surface moisture (Msurf)
    Msurf = (s11/s22)*((T0D - TD)/(TS - TD)) # Surface wetness  #######To be checked

    Msurf[Msurf > 1] = 0.9999
    Msurf[Msurf < 0] = 0.0001

    # Surface vapor pressure and deficit
    esurf = ea + Msurf*(esstar - ea)
    Dsurf = esurf - ea

    # Separating soil and canopy wetness to form a composite surface moisture
    Ms = Msurf.copy()
    Mcan = fc*Msurf
    Msoil = (1-fc)*Msurf

    TdewIndex = (TS - T0D)/(TA - TD) # TdewIndex > 1 signifies super dry condition
    Ep_PT = (1.26*slope*RN)/(slope + gamma) # Potential evaporation (Priestley-Taylor eqn.)                              

    # Surface wetness comes from the soil, vegetation contribution is neglegible
    Ms[(fc <= 0.25) & (TdewIndex < 1)] = Msoil[(fc <= 0.25) & (TdewIndex < 1)] 
    Mcan[(fc <= 0.25) & (TdewIndex < 1)] = 0

    Ms[(fc <= 0.25) & (TA > 10) & (TD < 0) & (Lnet < -125)] = Msoil[(fc <= 0.25) & (TA > 10) & (TD < 0) & (Lnet < -125)]
    Mcan[(fc <= 0.25) & (TA > 10) & (TD < 0) & (Lnet < -125)] = 0
    
    # Root zone moisture (Mrz)
    Mrz = (gamma*s11*(T0D - TD))/(slope*s33*(TS - TD) + gamma*s44*(TA - TD) - slope*s11*(T0D - TD))
    
    Mrz[Mrz > 1] = 0.9999
    Mrz[Mrz < 0] = 0.0001

    # Combine M to account for Hysteresis and initial estimation of surface vapor pressure
    M = Ms.copy()
    M[(Ep_PT > RN) & (dTS > 0)] = Mrz[(Ep_PT > RN) & (dTS > 0)]
    M[(Ep_PT > RN) & (fc <= 0.25)] = Mrz[(Ep_PT > RN) & (fc <= 0.25)]
    M[(Ep_PT > RN) & (Dsurf > DA)] = Mrz[(Ep_PT > RN) & (Dsurf > DA)] 

    M[(fc <= 0.25) & (dTS > 0) & (TA > 10) & (TD < 0) & (Lnet < -125)] = Mrz[(fc <= 0.25) & (dTS > 0) & (TA > 10) & (TD < 0) & (Lnet < -125)]
    M[(fc <= 0.25) & (dTS > 0) & (TA > 10) & (TD < 0) & (Dsurf > DA)] = Mrz[(fc <= 0.25) & (dTS > 0) & (TA > 10) & (TD < 0) & (Dsurf > DA)]
    M[(Ep_PT < RN) & (fc <= 0.25) & (Dsurf > DA)] = Mrz[(Ep_PT < RN) & (fc <= 0.25) & (Dsurf > DA)]

    es = ea + M*(esstar - ea)
    
    # vapor pressure deficit at surface 
    Ds = esstar - es 

    return (M,Mcan,Msoil,Ms,Mrz,s11,s22,s33,s44,es,T0D,Ds,TdewIndex)

def f_SoilMoisture_ITERATE(gamma,slope,s1,s2,s3,s4,TS,TA,dTS,TD,T0D,RG,RN,Lnet,fc,DA,D0,eastar,ea,e0star,esstar):
    '''
      This functions estimates the soil moisture availability (M) (or
      wetness) (value 0 to 1) based on thermal IR and meteorological
      information. However, this M will be treated as initial M, which will be
      later on estimated through iteration in the actual ET estimation loop to
      establish feedback between M and biophysical states
    '''
    # Surface moisture (Msurf)
    k = (e0star - ea)/(esstar - ea)
    Msurf = (s1/(k*s2))*((T0D - TD)/(TS - TD)) # Surface wetness

    Msurf[Msurf > 1] = 0.9999
    Msurf[Msurf < 0] = 0.0001
    
    # Separating soil and canopy wetness to form a composite surface moisture
    Ms = Msurf.copy()
    Mcan = fc*Msurf
    Msoil = (1-fc)*Msurf

    TdewIndex = (TS - T0D)/(TA - TD)
    Ep_PT = (1.26*slope*RN)/(slope + gamma) # Potential evaporation (Priestley-Taylor eqn.)                              

    Ms[(fc <= 0.25) & (TdewIndex < 1)] = Msoil[ (fc <= 0.25) & (TdewIndex < 1)]
    Ms[(fc <= 0.25) & (TA > 10) & (TD < 0) & (Lnet < -125)] = Msoil[(fc <= 0.25) & (TA > 10) & (TD < 0) & (Lnet < -125)]
    
    Mcan[(fc <= 0.25) & (TdewIndex < 1)] = 0
    Mcan[(fc <= 0.25) & (TA > 10) & (TD < 0) & (Lnet < -125)] = 0
    
    # Rootzone moisture (Mrz)                                
    Mrz = (gamma*s1*(T0D - TD))/(slope*s3*(TS - TD) + gamma*s4*(TA - TD) - slope*s1*(T0D - TD))
    
    Mrz[Mrz > 1] = 0.9999
    Mrz[Mrz < 0] = 0.0001

    TdewIndex = (TS - T0D)/(TA - TD)
    Ep_PT = (1.26*slope*RN)/(slope + gamma) # Potential evaporation (Priestley-Taylor eq.)  
    
    # Combine M to account for Hysteresis and initial estimation of surface vapor pressure
    M = Ms.copy()
    M[(Ep_PT > RN) & (dTS > 0) & (fc <= 0.25) & (D0 > DA) & (TdewIndex < 1)] = Mrz[(Ep_PT > RN) & (dTS > 0) & 
                                                              (fc <= 0.25) & (D0 > DA) & (TdewIndex < 1)]
    M[(Lnet < -125) & (dTS > 0) & (fc <= 0.25) & (D0 > DA) & (TA > 10) & (TD < 0)] = Mrz[(Lnet < -125) & (dTS > 0) & (fc <= 0.25) & 
                                                              (D0 > DA) & (TA > 10) & (TD < 0)]
    
    return (M,Ms,Mcan,Msoil,Mrz)