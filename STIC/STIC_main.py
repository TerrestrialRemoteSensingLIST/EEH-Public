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

# Main function

import argparse
import glob
import json
import os
import re
import sys
import traceback
from datetime import datetime

import h5py
import numpy as np

from ReadECOSTRESSData import Read_ECOSTRESS
from ReadAncillaryData import Read_Ancillary
from ReadERA5Data import Read_ERA5
from PySTIC import STIC
from TOA_Radiance import f_TOARadiance, f_ETDaily
from LUT import *  # noqa: F401,F403  (kept for parity with original imports)


# ---------------------------------------------------------------------------
# Local file-name caching (replaces S3_cache.cache_S3_with_pattern)
# ---------------------------------------------------------------------------
def build_local_pattern_map(directory, pattern, label=None):
    """
    Local, pure-Python replacement for cache_S3_with_pattern().

    Original S3 version read a pre-generated list of filenames from
    '<mount_folder>utils/<name>_S3_cached.txt' and applied re.search(pattern,
    filename) to each line, using capturing group 1 as the lookup key.

    This function reproduces the exact same key -> filename mapping logic
    (re.search, group(1)), but scans a local directory directly instead of
    reading a pre-generated S3 file list.

    Parameters
    ----------
    directory : str
        Local folder to scan.
    pattern : str
        Regular expression with exactly one capturing group, applied with
        re.search() (i.e. matched anywhere in the filename, not anchored).
    label : str, optional
        Purely cosmetic, used only for log messages.

    Returns
    -------
    dict : {captured_group: actual_filename}
    """
    if not os.path.isdir(directory):
        raise FileNotFoundError(f"Directory not found for '{label}': {directory}")

    regex = re.compile(pattern)
    file_map = {}
    for fname in sorted(os.listdir(directory)):
        if not os.path.isfile(os.path.join(directory, fname)):
            continue
        m = regex.search(fname)
        if not m:
            continue
        key = m.group(1)
        if key in file_map:
            print(f'Warning: duplicate key "{key}" for "{label}" -> keeping '
                  f'"{file_map[key]}", ignoring "{fname}"')
            continue
        file_map[key] = fname

    if not file_map:
        print(f'Warning: no files matching pattern for "{label}" found in {directory}')
    return file_map


# ---------------------------------------------------------------------------
# Raster export helper (kept for optional GIS export; requires GDAL)
# ---------------------------------------------------------------------------
def generate_raster_file(driver, filename, data, x_dim, y_dim,
                          geo_transform, proj_wkt,
                          no_data_value, data_type):
    """
    Description:
        Creates a raster file on disk for the specified data, using the
        specified driver.
    Note: It is assumed that the driver supports setting of the no data
          value. It is the caller's responsibility to fix it if it does not.
    Note: It is assumed that the caller specified the correct file
          extension in the filename parameter for the specified driver.
    """
    try:
        raster = driver.Create(filename, x_dim, y_dim, 1, data_type)
        raster.SetGeoTransform(geo_transform)
        raster.SetProjection(proj_wkt)
        raster.GetRasterBand(1).WriteArray(data)
        raster.GetRasterBand(1).SetNoDataValue(no_data_value)
        raster.FlushCache()
        del raster
    except Exception:
        print(f'Failed to generate file {filename}')


def adjust_lulc(lulc):
    lulc[(lulc == 80) | (lulc == 200)] = 0   # water
    lulc[(lulc == 111) | (lulc == 121)] = 1  # evergreen conifer
    lulc[(lulc == 112) | (lulc == 122)] = 1  # evergreen broadleaf
    lulc[(lulc == 113) | (lulc == 123)] = 1  # deciduous conifer
    lulc[(lulc == 114) | (lulc == 124)] = 1  # deciduous broadleaf
    lulc[(lulc == 115) | (lulc == 116) | (lulc == 125) | (lulc == 126)] = 1  # mixed forest
    lulc[lulc == 20] = 6   # woody savanna
    lulc[(lulc == 30) | (lulc == 100)] = 3   # grassland
    lulc[lulc == 90] = 5   # wetland
    lulc[lulc == 40] = 2   # cropland
    lulc[lulc == 50] = 8   # urban
    lulc[lulc == 70] = 9   # snow
    lulc[lulc == 60] = 7   # barren surfaces
    lulc[lulc == 0] = 0    # exclusion of invalid pixels
    return lulc


def generate_hdf5_file(directory_output, ET, H, G, Rn, gah, gsc, Ms, Mrz, ETD,
                        orbit_str, date_str, hour_str, min_str, sec_str):
    base_filename = ('EEH2STIC_L3_ET_' + orbit_str + '_' + date_str + 'T'
                      + hour_str + min_str + sec_str + '_0000_00.h5')
    filename = os.path.join(directory_output, base_filename)
    print("Writing output to: " + filename)
    f_h5 = h5py.File(filename, 'w')

    dset = f_h5.create_dataset('LE', data=np.float32(ET))
    dset.attrs['long_name'] = 'Latent heat flux'
    dset.attrs['units'] = 'W.m-2'
    dset.attrs['format'] = 'scaled'
    dset.attrs['coordsys'] = 'cartesian'
    dset.attrs['valid_range'] = np.array([0, 2000])
    dset.attrs['fill_value'] = -9999
    dset.attrs['scale_factor'] = 1
    dset.attrs['add_offset'] = 0

    dset = f_h5.create_dataset('H', data=np.float32(H))
    dset.attrs['long_name'] = 'Sensible heat flux'
    dset.attrs['units'] = 'W.m-2'
    dset.attrs['format'] = 'scaled'
    dset.attrs['coordsys'] = 'cartesian'
    dset.attrs['valid_range'] = np.array([0, 2000])
    dset.attrs['fill_value'] = -9999
    dset.attrs['scale_factor'] = 1
    dset.attrs['add_offset'] = 0

    dset = f_h5.create_dataset('G', data=np.float32(G))
    dset.attrs['long_name'] = 'Soil heat flux'
    dset.attrs['units'] = 'W.m-2'
    dset.attrs['format'] = 'scaled'
    dset.attrs['coordsys'] = 'cartesian'
    dset.attrs['valid_range'] = np.array([0, 2000])
    dset.attrs['fill_value'] = -9999
    dset.attrs['scale_factor'] = 1
    dset.attrs['add_offset'] = 0

    dset = f_h5.create_dataset('Rn', data=np.float32(Rn))
    dset.attrs['long_name'] = 'Net radiation'
    dset.attrs['units'] = 'W.m-2'
    dset.attrs['format'] = 'scaled'
    dset.attrs['coordsys'] = 'cartesian'
    dset.attrs['valid_range'] = np.array([0, 2000])
    dset.attrs['fill_value'] = -9999
    dset.attrs['scale_factor'] = 1
    dset.attrs['add_offset'] = 0

    dset = f_h5.create_dataset('gah', data=np.float32(gah))
    dset.attrs['long_name'] = 'Aerodynamic conductance'
    dset.attrs['units'] = 'm.s-1'
    dset.attrs['format'] = 'scaled'
    dset.attrs['coordsys'] = 'cartesian'
    dset.attrs['valid_range'] = np.array([0, 1])
    dset.attrs['fill_value'] = -9999
    dset.attrs['scale_factor'] = 1
    dset.attrs['add_offset'] = 0

    dset = f_h5.create_dataset('gsc', data=np.float32(gsc))
    dset.attrs['long_name'] = 'Surface conductance'
    dset.attrs['units'] = 'm.s-1'
    dset.attrs['format'] = 'scaled'
    dset.attrs['coordsys'] = 'cartesian'
    dset.attrs['valid_range'] = np.array([0, 1])
    dset.attrs['fill_value'] = -9999
    dset.attrs['scale_factor'] = 1
    dset.attrs['add_offset'] = 0

    dset = f_h5.create_dataset('Ms', data=np.float32(Ms))
    dset.attrs['long_name'] = 'Water stress'
    dset.attrs['units'] = 'unitless'
    dset.attrs['format'] = 'scaled'
    dset.attrs['coordsys'] = 'cartesian'
    dset.attrs['valid_range'] = np.array([0, 1])
    dset.attrs['fill_value'] = -9999
    dset.attrs['scale_factor'] = 1
    dset.attrs['add_offset'] = 0

    dset = f_h5.create_dataset('Mrz', data=np.float32(Mrz))
    dset.attrs['long_name'] = 'Water stress root zone'
    dset.attrs['units'] = 'unitless'
    dset.attrs['format'] = 'scaled'
    dset.attrs['coordsys'] = 'cartesian'
    dset.attrs['valid_range'] = np.array([0, 1])
    dset.attrs['fill_value'] = -9999
    dset.attrs['scale_factor'] = 1
    dset.attrs['add_offset'] = 0

    dset = f_h5.create_dataset('ETD', data=np.float32(ETD))
    dset.attrs['long_name'] = 'Daily ET'
    dset.attrs['units'] = 'mm.day-1'
    dset.attrs['format'] = 'scaled'
    dset.attrs['coordsys'] = 'cartesian'
    dset.attrs['valid_range'] = np.array([0, 2000])
    dset.attrs['fill_value'] = -9999
    dset.attrs['scale_factor'] = 1
    dset.attrs['add_offset'] = 0

    f_h5.close()
    print("STIC Model Running Completed")

    return base_filename


# ---------------------------------------------------------------------------
# Main processing loop
# ---------------------------------------------------------------------------
def run_STIC(lste_files, directory_geo, directory_cld, directory_fvc,
             directory_alb_mota, directory_lulc, directory_era5,
             directory_output, geo_pattern, cloud_pattern, fcover_pattern,
             mota_pattern, error_log_dir='.'):
    """
    Run the STIC model on a list of ECOSTRESS L2_LSTE files.

    Parameters
    ----------
    lste_files : list of str
        Paths to the ECOSTRESS L2_LSTE .h5 files to process
        (e.g. the EEHTES output of run_TES.py).
    directory_geo, directory_cld : str
        Directories containing the matching L1B_GEO / L2_CLOUD .h5 files.
    directory_fvc : str
        Directory containing the CGLS FCOVER .nc files.
    directory_alb_mota : str
        Directory containing the MCD43C3 (.hdf) albedo files, used for both
        directional and hemispherical albedo (mirrors the original script,
        where map_albdir2 == map_albhem2).
    directory_lulc : str
        Directory containing the PROBAV LC100 land-cover raster (fixed
        filename, no lookup pattern needed -- see Read_LULC()).
    directory_era5 : str
        Directory containing the ERA5 files.
    directory_output : str
        Directory where the output ET HDF5 files are written.
    geo_pattern, cloud_pattern, fcover_pattern, mota_pattern : str
        Regex patterns (with one capturing group) used to build the local
        file-name lookup maps. Defaults match the original ECOSTRESS/CGLS/
        MODIS naming conventions -- override only if your local filenames
        differ.
    error_log_dir : str
        Directory where the execution error JSON log is written, if any.
    """
    map_error = dict()

    print('Building local file-name caches...')
    map_geo = build_local_pattern_map(directory_geo, geo_pattern, label='L1B_GEO')
    map_cloud = build_local_pattern_map(directory_cld, cloud_pattern, label='CLOUD')
    map_fvc = build_local_pattern_map(directory_fvc, fcover_pattern, label='FCOVER')
    # NOTE: kept empty, mirroring the original script -- ReadAncillaryData's
    # Read_ALB_DIR/Read_ALB_HEM only use map_albdir1/map_albhem1 in a
    # permanently disabled (`if False:`) legacy branch (pre-2020.7 CGLS
    # Albedo_Directional/Albedo_Hemispherical .nc products).
    map_albdir1 = {}
    map_albhem1 = {}
    map_albdir2 = build_local_pattern_map(directory_alb_mota, mota_pattern, label='MOTA')
    map_albhem2 = map_albdir2
    print('Done caching.')

    for filename_lste in lste_files:
        print('Run STIC model on ECOSTRESS data ' + filename_lste)

        try:
            (lst_eco, lse_eco, lat_eco, lon_eco, _, watermask, cloudmask,
             year, month, day, hour, minute, second,
             date_str, hour_str, min_str, sec_str, orbit_str, _) = Read_ECOSTRESS(
                filename_lste, directory_geo, directory_cld, map_geo, map_cloud)
        except Exception:
            traceback.print_exc()
            print('Error while reading ECOSTRESS data -> stopping there')
            map_error[filename_lste] = traceback.format_exc()
            continue

        # Converting from K to Celsius degree
        lst_eco = lst_eco - 273.15

        # Calculating day of year and decimal time
        date = datetime(int(year), int(month), int(day))
        doy = date.timetuple().tm_yday
        time_decimal = float(hour) + float(minute) / 60. + float(second) / 3600.

        # Converting time from UTC to local solar time (in seconds)
        time_delta = lon_eco / 15.
        time_ls = float(hour) + float(minute) / 60. + float(second) / 3600. + time_delta
        time_ls[time_ls < 0] += 24
        time_ls[time_ls >= 24] = time_ls[time_ls >= 24] % 24
        time_ls = time_ls * 3600

        # Mask out pixels covered by cloud
        lst_eco[cloudmask == 1] = np.nan
        lse_eco[cloudmask == 1] = np.nan

        # Mask out pixels covered by water
        lst_eco[watermask == 1] = np.nan
        lse_eco[watermask == 1] = np.nan

        print('Reading ECOSTRESS data finished!')

        try:
            (fvc, alb_dir, alb_hem, lulc) = Read_Ancillary(
                directory_fvc, directory_alb_mota, directory_alb_mota, directory_lulc,
                lat_eco, lon_eco, year, month, day,
                map_fvc, map_albdir1, map_albdir2, map_albhem1, map_albhem2)
            lulc_adjusted = adjust_lulc(lulc)
            print('Reading ancillary CGLS data finished!')
        except Exception:
            traceback.print_exc()
            print('Error while reading ancillary data -> stopping there')
            map_error[filename_lste] = traceback.format_exc()
            continue

        try:
            (t_s, rh_s, _, sr_dir, sr_dif, ta_max) = Read_ERA5(
                directory_era5, lat_eco, lon_eco, year, month, day, hour, minute, second)
            t_s = t_s - 273.15   # Converting from K to Celsius degree
            rh_s = rh_s * 100    # Converting from 0-1 to percentage
            print('Reading ERA5 data finished!')
        except Exception:
            traceback.print_exc()
            print('Error while reading ERA5 data -> stopping there')
            map_error[filename_lste] = traceback.format_exc()
            continue

        try:
            (RN_STIC, G_STIC, H_STIC, LE_STIC, gah_STIC, gsc_STIC,
             T0_STIC, Ms_STIC, Mrz_STIC, converged, Lin, Lout, Lnet, Sin, Sout) = STIC(
                lst_eco, lse_eco, t_s, rh_s, sr_dir, sr_dif, alb_dir, alb_hem, fvc, time_ls)
            print('Running the STIC model finished!')

            # Mask out the invalid pixels
            RN_STIC[np.isnan(RN_STIC)] = -9999
            G_STIC[np.isnan(G_STIC)] = -9999
            H_STIC[np.isnan(H_STIC)] = -9999
            LE_STIC[np.isnan(LE_STIC)] = -9999
            gah_STIC[np.isnan(gah_STIC)] = -9999
            gsc_STIC[np.isnan(gsc_STIC)] = -9999
            Ms_STIC[np.isnan(Ms_STIC)] = -9999
            Mrz_STIC[np.isnan(Mrz_STIC)] = -9999

            # Mask out the non-converged pixels
            RN_STIC[~converged] = -9999
            G_STIC[~converged] = -9999
            H_STIC[~converged] = -9999
            LE_STIC[~converged] = -9999
            gah_STIC[~converged] = -9999
            gsc_STIC[~converged] = -9999
            Ms_STIC[~converged] = -9999
            Mrz_STIC[~converged] = -9999

            # Mask out the night-time pixels
            G_STIC[RN_STIC < 0] = -9999
            H_STIC[RN_STIC < 0] = -9999
            LE_STIC[RN_STIC < 0] = -9999
            gah_STIC[RN_STIC < 0] = -9999
            gsc_STIC[RN_STIC < 0] = -9999
            Ms_STIC[RN_STIC < 0] = -9999
            Mrz_STIC[RN_STIC < 0] = -9999
            RN_STIC[RN_STIC < 0] = -9999

            # Mask out ET for non-vegetated land surface types
            LE_STIC[(lulc_adjusted < 1) | (lulc_adjusted > 6)] = -9999
            H_STIC[(lulc_adjusted < 1) | (lulc_adjusted > 6)] = -9999
            G_STIC[(lulc_adjusted < 1) | (lulc_adjusted > 6)] = -9999
            RN_STIC[(lulc_adjusted < 1) | (lulc_adjusted > 6)] = -9999
            gah_STIC[(lulc_adjusted < 1) | (lulc_adjusted > 6)] = -9999
            gsc_STIC[(lulc_adjusted < 1) | (lulc_adjusted > 6)] = -9999
            Ms_STIC[(lulc_adjusted < 1) | (lulc_adjusted > 6)] = -9999
            Mrz_STIC[(lulc_adjusted < 1) | (lulc_adjusted > 6)] = -9999

            # Calculate daily ET from the instantaneous estimates
            (RgTOAiMJ, RgTOAInt) = f_TOARadiance(doy, time_decimal, lat_eco, lon_eco)
            # Unit conversion (W.m-2 -> mm over a 30-min / 1800s window,
            # using the latent heat of vaporization, 2264.76 kJ/kg)
            LE_STIC1 = LE_STIC * 30 * 60 / (2264.76 * 1e3)
            ET_DLY = f_ETDaily(LE_STIC1, float(hour), t_s, ta_max, Lin, Lout, Sin, Sout,
                                RgTOAiMJ, RgTOAInt, lulc_adjusted)
            ET_DLY[LE_STIC == -9999] = -9999

            # Output the variables in HDF5 file
            generate_hdf5_file(directory_output, LE_STIC, H_STIC, G_STIC, RN_STIC,
                                gah_STIC, gsc_STIC, Ms_STIC, Mrz_STIC, ET_DLY,
                                orbit_str, date_str, hour_str, min_str, sec_str)
            print('Outputting to HDF5 file finished')
        except Exception:
            traceback.print_exc()
            print('Error -> stopping there')
            map_error[filename_lste] = traceback.format_exc()

    if map_error:
        pid = str(os.getpid())
        os.makedirs(error_log_dir, exist_ok=True)
        error_log_path = os.path.join(error_log_dir, f'execution_{pid}_errors.txt')
        with open(error_log_path, 'w') as file:
            file.write(json.dumps(map_error))
        print(f'Some files failed to process, see {error_log_path}')


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def build_arg_parser():
    parser = argparse.ArgumentParser(
        prog='STIC_main.py',
        description='Run the STIC model for ET estimation on ECOSTRESS L2_LSTE data.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        '--input-files', '-i', nargs='+', metavar='LSTE_FILE',
        help='One or more ECOSTRESS L2_LSTE .h5 files to process (e.g. EEHTES output).'
    )
    input_group.add_argument(
        '--input-dir', metavar='DIR',
        help='Directory to scan for L2_LSTE .h5 files (used together with --pattern).'
    )
    parser.add_argument(
        '--pattern', default='*L2_LSTE*.h5',
        help='Glob pattern used to select LSTE files when --input-dir is given.'
    )

    parser.add_argument('--geo-dir', required=True,
                         help='Directory containing the ECOSTRESS L1B_GEO .h5 files.')
    parser.add_argument('--cloud-dir', required=True,
                         help='Directory containing the ECOSTRESS L2_CLOUD .h5 files.')
    parser.add_argument('--fcover-dir', required=True,
                         help='Directory containing the CGLS FCOVER .nc files.')
    parser.add_argument('--albedo-mota-dir', required=True,
                         help='Directory containing the MCD43C3 (.hdf) albedo files '
                              '(used for both directional and hemispherical albedo).')
    parser.add_argument('--lulc-dir', required=True,
                         help='Directory containing the PROBAV LC100 land-cover raster.')
    parser.add_argument('--era5-dir', required=True,
                         help='Directory containing the ERA5 reanalysis files.')
    parser.add_argument('--output-dir', required=True,
                         help='Directory where the output ET HDF5 files are written.')

    pattern_group = parser.add_argument_group(
        'File-name matching patterns (override only if your local filenames '
        'differ from the standard ECOSTRESS/CGLS/MODIS conventions)')
    pattern_group.add_argument(
        '--geo-pattern', default=r'(ECOv002_L1B_GEO.*)_\d{4}_\d{2}\.h5')
    pattern_group.add_argument(
        '--cloud-pattern', default=r'(ECOv002_L2_CLOUD.*)_\d{4}_\d{2}\.h5')
    pattern_group.add_argument(
        '--fcover-pattern', default=r'c_gls_FCOVER300(?:-RT\d+)?_(\d{8})0000_GLOBE_.*\.nc')
    pattern_group.add_argument(
        '--mota-pattern', default=r'MCD43C3\.A(\d{7}).*\.hdf')

    parser.add_argument('--error-log-dir', default='.',
                         help='Directory to write the execution error JSON log to, if any.')

    return parser


def main(argv=None):
    parser = build_arg_parser()
    args = parser.parse_args(argv)

    if args.input_files:
        lste_files = args.input_files
    else:
        lste_files = sorted(glob.glob(os.path.join(args.input_dir, args.pattern)))
        if not lste_files:
            parser.error(f"No files matching '{args.pattern}' found in {args.input_dir}")

    os.makedirs(args.output_dir, exist_ok=True)

    run_STIC(
        lste_files=lste_files,
        directory_geo=args.geo_dir,
        directory_cld=args.cloud_dir,
        directory_fvc=args.fcover_dir,
        directory_alb_mota=args.albedo_mota_dir,
        directory_lulc=args.lulc_dir,
        directory_era5=args.era5_dir,
        directory_output=args.output_dir,
        geo_pattern=args.geo_pattern,
        cloud_pattern=args.cloud_pattern,
        fcover_pattern=args.fcover_pattern,
        mota_pattern=args.mota_pattern,
        error_log_dir=args.error_log_dir,
    )


if __name__ == '__main__':
    main()