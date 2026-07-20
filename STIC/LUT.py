# -*- coding: utf-8 -*-

import numpy as np

#Calculate daily ET using the calibrated LUT method
def f_ETDaily(ET, hour, ta, ta_max, lwd, lwu, swd, swu, RgTOAiMJ, RgTOAInt, biome):
    
    #Adjust the early and late hours
    hour = max(8, min(hour, 16))
    
    #ET upscaling using the LUT method
    #1 for forest, 2 for cropland, 3 for grassland, 4 for shrubland, 5 for wetland, 6 for savanna
    rn = swd - swu + lwd - lwu
    
    #Initialize output array with NaNs or zeros
    ET_DLY = np.full_like(ET, -9999)  # or use np.zeros_like if appropriate

    #Forest (biome == 1)
    forest = (biome == 1)
    forest_early = forest & np.isin(hour, [8, 9, 10, 11, 12])
    forest_late  = forest & np.isin(hour, [13, 14, 15, 16])

    ET_DLY[forest_early] = (
        ET[forest_early] * RgTOAInt[forest_early] / RgTOAiMJ[forest_early] *
        (rn[forest_early] / (swd[forest_early] - swu[forest_early] + lwd[forest_early]))
    )

    ET_DLY[forest_late] = (
        ET[forest_late] * RgTOAInt[forest_late] / RgTOAiMJ[forest_late]
    )

    # Wetland (biome == 5)
    wetland = (biome == 5)
    wetland_lwd = wetland & np.isin(hour, [8, 9, 10, 15, 16])
    wetland_ta  = wetland & np.isin(hour, [11, 12, 13, 14])

    ET_DLY[wetland_lwd] = (
        ET[wetland_lwd] * RgTOAInt[wetland_lwd] / RgTOAiMJ[wetland_lwd] *
        (lwd[wetland_lwd] / lwu[wetland_lwd])
    )

    ET_DLY[wetland_ta] = (
        ET[wetland_ta] * RgTOAInt[wetland_ta] / RgTOAiMJ[wetland_ta] *
        (ta[wetland_ta] / ta_max[wetland_ta])
    )

    # Savanna (biome == 6)
    savanna = (biome == 6)
    savanna_rn  = savanna & (hour == 8)
    savanna_lwd = savanna & np.isin(hour, [9, 10])
    savanna_ta  = savanna & np.isin(hour, [11, 12, 13, 14, 15, 16])

    ET_DLY[savanna_rn] = (
        ET[savanna_rn] * RgTOAInt[savanna_rn] / RgTOAiMJ[savanna_rn] *
        (rn[savanna_rn] / (swd[savanna_rn] - swu[savanna_rn] + lwd[savanna_rn]))
    )

    ET_DLY[savanna_lwd] = (
        ET[savanna_lwd] * RgTOAInt[savanna_lwd] / RgTOAiMJ[savanna_lwd] *
        (lwd[savanna_lwd] / lwu[savanna_lwd])
    )

    ET_DLY[savanna_ta] = (
        ET[savanna_ta] * RgTOAInt[savanna_ta] / RgTOAiMJ[savanna_ta] *
        (ta[savanna_ta] / ta_max[savanna_ta])
    )

    # Grassland, cropland, shrubland
    others = ~(forest | wetland | savanna)

    ET_DLY[others] = (
        ET[others] * RgTOAInt[others] / RgTOAiMJ[others]
    )
    
    return ET_DLY
