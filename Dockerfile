FROM python:3.10-bookworm

# ---- system dependencies (RTTOV) ----
RUN apt-get update && apt-get install -y --no-install-recommends \
        gfortran make perl \
        libhdf5-dev \
    && rm -rf /var/lib/apt/lists/*

# numpy must be installed before RTTOV f2py wrapper build
RUN pip install --no-cache-dir numpy

# ---- RTTOV 12.3 (from user-supplied tarball) ----
COPY rttov123.tar.gz /opt/
RUN mkdir -p /opt/rttov123 \
    && cd /opt/rttov123 \
    && tar xzf /opt/rttov123.tar.gz \
    && rm /opt/rttov123.tar.gz

# Patch RTTOV wrapper for Python 3.10+ (collections.Iterable → collections.abc.Iterable)
RUN sed -i 's/collections\.Iterable/collections.abc.Iterable/g' /opt/rttov123/wrapper/pyrttov/__init__.py

RUN cd /opt/rttov123/src \
    && ../build/Makefile.PL RTTOV_HDF=0 RTTOV_F2PY=1 RTTOV_USER_LAPACK=0 \
    && make ARCH=gfortran INSTALLDIR=./

ENV RTTOV_INSTALLDIR=/opt/rttov123
ENV PYTHONPATH="/opt/rttov123/wrapper:${PYTHONPATH}"

# ---- GDAL 3.11.1 from source (needed for COG INTERLEAVE=BAND) ----
RUN apt-get update && apt-get install -y --no-install-recommends \
        cmake swig libproj-dev libgeos-dev libcurl4-openssl-dev \
        libsqlite3-dev libtiff-dev libgeotiff-dev libopenjp2-7-dev \
        libpng-dev libjpeg-dev zlib1g-dev pkg-config \
    && rm -rf /var/lib/apt/lists/*

RUN curl -fsSL https://github.com/OSGeo/gdal/releases/download/v3.11.1/gdal-3.11.1.tar.gz \
        | tar xz -C /tmp \
    && cd /tmp/gdal-3.11.1 \
    && mkdir build && cd build \
    && cmake .. -DCMAKE_INSTALL_PREFIX=/usr \
        -DCMAKE_BUILD_TYPE=Release \
        -DBUILD_APPS=ON \
        -DBUILD_PYTHON_BINDINGS=OFF \
        -DGDAL_BUILD_OPTIONAL_DRIVERS_HDF5=ON \
    && make -j$(nproc) \
    && make install \
    && ldconfig \
    && rm -rf /tmp/gdal-3.11.1

# ---- Python dependencies ----
RUN pip install --no-cache-dir \
        h5py scipy cfgrib \
        opencv-python-headless \
        netCDF4 pyhdf \
        GDAL==$(gdal-config --version) \
        pandas xarray rioxarray rasterio \
        matplotlib shapely geopandas \
        h5netcdf pyresample \
        python-dotenv requests cdsapi eumdac boto3

# ---- smoke test: verify RTTOV + all imports work ----
RUN python -c "\
import sys, os; sys.path.insert(0,'/opt/rttov123/wrapper'); \
import pyrttov; \
import h5py, numpy, scipy, cv2, cfgrib; \
from osgeo import gdal; \
import pandas, xarray, rasterio, shapely, geopandas; \
print('OK — all dependencies verified, RTTOV', pyrttov.pyrttov_version if hasattr(pyrttov,'pyrttov_version') else 'loaded'); \
print('GDAL', gdal.__version__); \
os._exit(0)"

# ---- application code ----
WORKDIR /app
COPY TES/ TES/
COPY STIC/ STIC/
COPY GPP/ GPP/
COPY tools/ tools/
COPY run_pipeline.py .
COPY .env.example .env

ENTRYPOINT ["python", "run_pipeline.py"]
