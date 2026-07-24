# Scope
This code is to estimate instantaneous and daily evapotranspiration from the ECOSTRESS LST and ancillary data using the STIC model and the temporal upscaling method.

# Model
The original STIC (STIC 1.0 - 1.3) model was developed by Mallick et al. (2014, 2016, 2018). Hu and Mallick et al. (2023) adapted the STIC model (Matlab version) to implement it on the ECOSTRESS data (Python version, STICPy).

More information about the model can be found in the papers below:
1) Mallick, Kaniska, Jarvis A.J., Boegh E., Fisher J.B., Drewry D.T., et al. (2014). A Surface Temperature Initiated Closure (STIC) for surface energy balance fluxes. Remote Sensing of Environment, 141, 243-61. https://doi.org/10.1016/j.rse.2013.10.022
2) Mallick, Kaniska, Trebs, I., Boegh, E., Giustarini, L., Schlerf, M., Drewry, D. T., et al. (2016). Canopy-scale biophysical controls of transpiration and evaporation in the Amazon Basin. Hydrology and Earth System Sciences, 20(10), 4237–4264. https://doi.org/10.5194/hess-20-4237-2016
3) Mallick, Kaniska, Toivonen, E., Trebs, I., Boegh, E., Cleverly, J., Eamus, D., et al. (2018). Bridging Thermal Infrared Sensing and Physically‐Based Evapotranspiration Modeling: From Theoretical Implementation to Validation Across an Aridity Gradient in Australian Ecosystems. Water Resources Research, 54(5), 3409–3435. https://doi.org/10.1029/2017WR021357
4) Hu, Tian, Mallick, Kaniska, Hitzelberger, P., Didry, Y., Boulet, G., Szantoi, Z., et al. (2023). Evaluating European ECOSTRESS Hub Evapotranspiration Products Across a Range of Soil-Atmospheric Aridity and Biomes Over Europe. Water Resources Research, 59(8), e2022WR034132. https://doi.org/https://doi.org/10.1029/2022WR034132

# Inputs
This step used to require a fixed 9-column input.txt config file, processed
one line per image (see legacy format below). It has been replaced by
command-line arguments (see Usage below):

1) `--input-files` / `--input-dir` : ECOSTRESS L2_LSTE file(s) to process (TES output)
2) `--geo-dir`          : directory for ECOSTRESS L1B_GEO
3) `--cloud-dir`        : directory for ECOSTRESS cloud mask
4) `--fcover-dir`       : directory for FVC data
5) `--albedo-mota-dir`  : directory for blacksky and whitesky albedo
6) `--lulc-dir`         : directory for LULC data
7) `--era5-dir`         : directory for ERA5 data
8) `--output-dir`       : directory for storing output HDF5 data

# Usage
To run the code, type the following in command line:

python myProcessor.py \
    --input-files /data/output_tes/EEH2TES_L2_LSTE_..._0000_00.h5 \
    --geo-dir /data/geo \
    --cloud-dir /data/cloud \
    --fcover-dir /data/fcover \
    --albedo-mota-dir /data/mota \
    --lulc-dir /data/lulc \
    --era5-dir /data/era5 \
    --output-dir /data/output_stic
    
# Citation
The following paper should be cited when using this code:

Hu, Tian, Mallick, Kaniska, Hitzelberger, P., Didry, Y., Boulet, G., Szantoi, Z., et al. (2023). Evaluating European ECOSTRESS Hub Evapotranspiration Products Across a Range of Soil-Atmospheric Aridity and Biomes Over Europe. Water Resources Research, 59(8), e2022WR034132. https://doi.org/https://doi.org/10.1029/2022WR034132