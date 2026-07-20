# -*- coding: utf-8 -*-

import numpy as np

def f_Twet(eastar,slopeTA,TA,slopeTD,TD,gamma):                           
    # Initialize
    alfaWET = 1.26 # wet-surface Priestley-Taylor coefficient
    eWETstar = eastar.copy() # saturation vapor pressure of wet surface
    sWET = slopeTA.copy() # slope of saturation vapor presure for wet surface temperature

    RHaWET = 0.98 # relative humidity air for wet surface
    taWB = (gamma*TA + TD*slopeTD)/(gamma + slopeTD) # wet bulb temperature
    eastarWB = 6.13753*np.exp((17.27*taWB)/(taWB + 237.3)) # saturation vapor pressure at wet bulb temperature
    slopeWB = 4098.*eastarWB/(taWB + 237.3)**2 # slope of saturation vapor presure for wet surface temperature
    eaWB = eastarWB*RHaWET # vapour pressure at wet bulb  temperature
    tdWB = 237.3*np.log(eaWB/6.13753)/(17.27 - np.log(eaWB/6.13753)) # dewpoint temperature at wet bulb temperature
    DAwet = eastarWB - eaWB # vapour pressure deficit at wet bulb temperature
    
    TwetNEW = TA.copy()
    TwetOLD = TA.copy()
    error_Twet = 0.05
    steps = 0
    
    while(np.amax(error_Twet) > 0.01 and steps < 10000):  
        TwetNEW = taWB + ((eWETstar - eaWB)*(slopeWB + gamma - alfaWET*slopeWB)/(alfaWET*slopeWB*gamma)) 

#         # New formulation of STIC 1.3, combining PT and Bowen ratio equation for wet surface
#         TwetNEW = TA + ((eWETstar - ea)*(sWET + gamma - alfaWET*sWET)/(alfaWET*sWET*gamma))

        eWETstar = 6.13753*np.exp((17.27*TwetNEW)/(TwetNEW + 237.3))
        sWET = 4098*eWETstar/(TwetNEW + 237.3)**2

        # Error
        error_Twet = np.abs(TwetOLD - TwetNEW)
        TwetOLD = TwetNEW
        steps = steps + 1  

    print('Number of iterations for Twet calculation: ',steps)

    Twet = TwetNEW.copy()
    dTwet = Twet - taWB
    
    return (Twet,eWETstar,sWET,taWB,tdWB,eastarWB,slopeWB,eaWB,dTwet,DAwet)

def f_Radiation_wetsurface(RG,taWB,eaWB,Twet,sigma):   
    etaaWET = 1.24*(eaWB/(taWB + 273.15))**(1./7) # air emissivity above wet surface
    etasWET = 0.99 # surface emissivity of a wet surface
    RSn = RG - (0.06*RG) # Net shortwave radiation, with 6% albedo of water
    Lin_wet = sigma*etaaWET*(taWB + 273.15)**4 # Downwelling Longwave over a wet surface
    Lout_wet = sigma*etasWET*(Twet + 273.15)**4 # Upwelling Longwave over a wet surface
    Lnet_wet = etasWET*Lin_wet - Lout_wet # Net Longwave over a wet surface
    RN_wet = RSn + Lnet_wet # Net radiation over a wet surface
    
    return (Lin_wet,Lout_wet,Lnet_wet,RN_wet)

def f_G_PHI_wetsurface(RN_wet,ttSEC,Mwet):
    cgMIN = 0.05 # for wet surface
    cgMAX = 0.10 # for dry surface 
    tgMIN = 74000 # for wet surface
    tgMAX = 80000 # for dry surface 
    solNooN = 12.*60.*60
    tg0 = solNooN - ttSEC

    # Estimating GHF according to Santanello and Friedl (2003)
    cg = (1-Mwet)*cgMAX + Mwet*cgMIN # considers G and net radiation weighted by soil wetness
    tg = (1-Mwet)*tgMAX + Mwet*tgMIN # the lag between max. net radiation and G weighted by soil wetness

    G_wet = RN_wet*cg*np.cos(2.*np.pi*(tg0 + 10800)/tg) # Wet surface soil heat flux
    G_wet[RN_wet < 0] = -G_wet[RN_wet < 0]

    PHI_wet = RN_wet - G_wet
    
    return (cg,tg,tg0,G_wet,PHI_wet)

def f_NetRadiation(sigma,TA,ea,TS,etas,sr_dir,sr_dif,alb_dir,alb_hem):
    etaa = 1.24*(ea/(TA + 273.15))**(1./7) # air emissivity
    Lin = sigma*etaa*(TA + 273.15)**4
    Lout = sigma*etas*(TS + 273.15)**4
    Lnet = etas*Lin - Lout
    

    Sin = sr_dir + sr_dif
    Sout = sr_dir*alb_dir + sr_dif*alb_hem
    Snet = Sin - Sout

    RN = sr_dir*(1 - alb_dir) + sr_dif*(1 - alb_hem) + Lnet
    
    return (etaa,Lin,Lout,Lnet,Sin,Sout,Snet,RN)

def f_G_PHI_actualsurface(RN,RNsoil,ttSEC,M):
    cgMIN = 0.05 # for wet surface
    cgMAX = 0.35 # for dry surface, in water-controlled ecosystems
#     cgMAX = 0.20 # for dry surface, in radiation-controlled ecosystems

    tgMIN = 74000 # for wet surface
    tgMAX = 100000 # for dry surface 
    
    solNooN = 12*60*60
    tg0 = solNooN - ttSEC

    # Estimating GHF according to Santanello and Friedl (2003)
    cg = (1-M)*cgMAX + M*cgMIN
    tg = (1-M)*tgMAX + M*tgMIN

    #G = RN*cg*np.cos(2*np.pi*(tg0 + 10800)/tg)
    #G[RN < 0] = -G[RN < 0]

    G = RNsoil*cg*np.cos(2*np.pi*(tg0 + 10800)/tg) #Correction on 12/11/2021
    G[RNsoil < 0] = -G[RNsoil < 0] #Correction on 12/11/2021

    PHI = RN - G
    
    return (cg,tg,tg0,G,PHI)
