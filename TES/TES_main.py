#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Script to run the Temperature Emissivity Separation (TES) algorithm for LST estimation

Created on April 15 2021
@author: Tian Hu, Yoanne Didry at LIST

© 2026 – Luxembourg Institute of Science and Technology
Authors : Tian Hu (tian.hu@list.lu), Yoanne Didry (yoanne.didry@list.lu)
Code licensed under MIT
SPDX-License-Identifier: MIT
"""

# Main function

import argparse
import glob
import hashlib
import json
import os
import re
import traceback

import cv2
import h5py
import numpy as np
from scipy.interpolate import griddata

from AtmCorrection import runRTTOV
from ReadECOL1BData import Read_L1B_Data
from ReadERA5Data import Read_ERA5
from TES_vec import LST_Estimate


# Default TES coefficients 
DEFAULT_ALPHA1 = 0.9895
DEFAULT_ALPHA2 = 0.7994
DEFAULT_ALPHA3 = 0.8572


def _cache_file_path(cache_dir, label, directory, key_index1, key_index2):
    """Build a stable, unique cache filename for a given (directory, key slice) pair."""
    key = hashlib.sha1(f"{directory}|{key_index1}|{key_index2}".encode('utf-8')).hexdigest()[:16]
    safe_label = re.sub(r'[^A-Za-z0-9_-]', '_', label or 'cache')
    return os.path.join(cache_dir, f'geo_cache_{safe_label}_{key}.json')


def build_local_geo_map(directory_geo, key_index1=0, key_index2=41,
                         exclude_substrings=('.xml', '.dmrpp'),
                         cache_dir=None, force_rebuild=False):
    """
    This function scans a local directory directly and reproduces that
    exact keying. Results are optionally cached to disk (JSON) to avoid
    rescanning slow/networked storage on every run.

    When multiple files share the same key, the file whose name is
    lexicographically greatest is kept deterministically.

    Parameters
    ----------
    directory_geo : str
        Local folder containing the GEO .h5 files.
    key_index1, key_index2 : int
        Character slice [key_index1:key_index2] of the filename used as
        the lookup key. Defaults (0, 41) match the ECOSTRESS GEO naming
        convention and Read_L1B_Data()'s key construction.
    exclude_substrings : tuple of str
        Files containing any of these substrings are skipped before any
        further processing (filters out sidecar metadata files such as
        .xml/.dmrpp).
    cache_dir : str, optional
        Directory used to persist the resulting {key: filename} map to
        disk (JSON). If None, no caching is performed (directory is
        rescanned every call).
    force_rebuild : bool
        If True, ignore any existing cache file and rescan the directory.

    Returns
    -------
    dict : {filename[key_index1:key_index2]: actual_filename}
    """
    cache_path = None
    if cache_dir:
        os.makedirs(cache_dir, exist_ok=True)
        cache_path = _cache_file_path(cache_dir, 'GEO', directory_geo, key_index1, key_index2)
        if not force_rebuild and os.path.isfile(cache_path):
            try:
                with open(cache_path, 'r') as f:
                    geo_map = json.load(f)
                print(f'[GEO] Using cached file-name map ({len(geo_map)} entries) from {cache_path}')
                return geo_map
            except Exception as e:
                print(f'[GEO] Warning: failed to read cache {cache_path} ({e}) -> rescanning directory')

    if not os.path.isdir(directory_geo):
        raise FileNotFoundError(f"GEO directory not found: {directory_geo}")

    print(f'[GEO] Scanning {directory_geo} (this may take a while on slow/networked storage)...')

    # Sort entries by name first for a deterministic scan order, independent
    # of the underlying filesystem's (unordered) os.scandir() iteration order.
    with os.scandir(directory_geo) as it:
        entries = sorted(it, key=lambda e: e.name)

    geo_map = {}
    for entry in entries:
        fname = entry.name
        # Cheap in-memory string checks first: avoids a stat() syscall
        # (potentially slow on networked/mounted storage) for every entry
        # that is trivially excluded anyway.
        if any(sub in fname for sub in exclude_substrings):
            continue
        if not entry.is_file(follow_symlinks=False):
            continue
        key = fname[key_index1:key_index2]
        if key in geo_map:
            if fname > geo_map[key]:
                print(f'[GEO] Note: preferring "{fname}" over "{geo_map[key]}" for key "{key}"')
                geo_map[key] = fname
            continue
        geo_map[key] = fname

    if not geo_map:
        print(f'[GEO] Warning: no GEO files found in {directory_geo}')
    else:
        print(f'[GEO] Found {len(geo_map)} matching files in {directory_geo}')

    if cache_path:
        try:
            with open(cache_path, 'w') as f:
                json.dump(geo_map, f)
            print(f'[GEO] Cached {len(geo_map)} entries to {cache_path}')
        except Exception as e:
            print(f'[GEO] Warning: failed to write cache {cache_path}: {e}')

    return geo_map


def generate_hdf5_file(directory_output, lst, emib2, emib4, emib5, bbe, mask, qa,
                        orbit_str, date_str, hour_str, min_str, sec_str):
    lst_scaled = np.uint16(lst / 0.02)
    emi_b2_scaled = np.uint8((emib2 - 0.489999) / 0.002)
    emi_b4_scaled = np.uint8((emib4 - 0.489999) / 0.002)
    emi_b5_scaled = np.uint8((emib5 - 0.489999) / 0.002)
    bbe_scaled = np.uint8((bbe - 0.489999) / 0.002)
    lst_scaled[mask] = 0
    emi_b2_scaled[mask] = 0
    emi_b4_scaled[mask] = 0
    emi_b5_scaled[mask] = 0
    bbe_scaled[mask] = 0

    base_filename = ('EEH2TES_L2_LSTE_' + orbit_str + '_' + date_str + 'T'
                      + hour_str + min_str + sec_str + '_0000_00.h5')
    filename = os.path.join(directory_output, base_filename)

    with h5py.File(filename, 'w') as f_lst:
        dset_names = ['Emis2', 'Emis4', 'Emis5', 'BBE', 'LST', 'qa']
        datasets = {
            'Emis2': emi_b2_scaled,
            'Emis4': emi_b4_scaled,
            'Emis5': emi_b5_scaled,
            'BBE': bbe_scaled,
            'LST': lst_scaled,
            'qa': np.int16(qa)
        }
        attributes = {
            'Emis2': {'long_name': 'Band 2 Emissivity', 'units': 'n/a', 'format': 'scaled', 'coordsys': 'cartesian', 'valid_range': [1, 255], 'fill_value': 0, 'scale_factor': 0.002, 'add_offset': 0.489999},
            'Emis4': {'long_name': 'Band 4 Emissivity', 'units': 'n/a', 'format': 'scaled', 'coordsys': 'cartesian', 'valid_range': [1, 255], 'fill_value': 0, 'scale_factor': 0.002, 'add_offset': 0.489999},
            'Emis5': {'long_name': 'Band 5 Emissivity', 'units': 'n/a', 'format': 'scaled', 'coordsys': 'cartesian', 'valid_range': [1, 255], 'fill_value': 0, 'scale_factor': 0.002, 'add_offset': 0.489999},
            'BBE': {'long_name': 'Broad Band Emissivity', 'units': 'n/a', 'format': 'scaled', 'coordsys': 'cartesian', 'valid_range': [1, 255], 'fill_value': 0, 'scale_factor': 0.002, 'add_offset': 0.489999},
            'LST': {'long_name': 'Land Surface Temperature', 'units': 'K', 'format': 'scaled', 'coordsys': 'cartesian', 'valid_range': [7500, 65535], 'fill_value': 0, 'scale_factor': 0.02, 'add_offset': 0},
            'qa': {'long_name': 'Quality Assurance', 'units': 'N/A', 'format': 'scaled', 'coordsys': 'cartesian', 'valid_range': [-5, 5], 'fill_value': -9999, 'scale_factor': 1, 'add_offset': 0}
        }
        for name in dset_names:
            dset = f_lst.create_dataset(name, data=datasets[name], compression='gzip', compression_opts=4)
            for key, value in attributes[name].items():
                dset.attrs[key] = value

    return base_filename


def run_TES(rad_files, directory_geo, directory_era5, directory_output,
            rttov_installdir, rttov_wrapper_dir, rttov_coef_file=None,
            alpha1=DEFAULT_ALPHA1, alpha2=DEFAULT_ALPHA2, alpha3=DEFAULT_ALPHA3,
            rttov_lib_preload=None, cache_dir=None, force_rebuild_cache=False,
            error_log_dir='.'):
    """
    Run the TES algorithm on a list of ECOSTRESS L1B RAD files, using a
    single (alpha1, alpha2, alpha3) coefficient set.

    Parameters
    ----------
    rad_files : list of str
        Paths to the ECOSTRESS L1B RAD .h5 files to process.
    directory_geo : str
        Directory containing the matching L1B GEO .h5 files.
    directory_era5 : str
        Directory containing the ERA5 reanalysis files.
    directory_output : str
        Directory where the output LST HDF5 files are written.
    rttov_installdir : str
        RTTOV installation directory (forwarded to AtmCorrection.runRTTOV).
    rttov_wrapper_dir : str
        Directory of the RTTOV python wrapper (pyrttov).
    rttov_coef_file : str, optional
        Explicit path to the RTTOV coefficient file.
    alpha1, alpha2, alpha3 : float
        TES coefficients used in LST_Estimate(). Defaults match the
        original single-coefficient-set variant of this script.
    rttov_lib_preload : str, optional
        Optional explicit path to a compiled RTTOV f2py wrapper shared
        library (.so) to preload via ctypes before importing pyrttov.
    cache_dir : str, optional
        Directory used to persist the GEO file-name lookup map to disk
        (JSON), avoiding a costly re-scan of slow/networked storage on
        every run. If None, caching is disabled.
    force_rebuild_cache : bool
        If True, ignore any existing cache and rescan the GEO directory.
    error_log_dir : str
        Directory where the execution error JSON log is written, if any.
    """
    if rttov_lib_preload:
        import ctypes
        print(f'Preloading RTTOV wrapper library from {rttov_lib_preload}')
        ctypes.CDLL(rttov_lib_preload)

    map_error = dict()

    map_geo = build_local_geo_map(directory_geo, cache_dir=cache_dir,
                                   force_rebuild=force_rebuild_cache)

    for filename_rad in rad_files:
        print('Run Temperature Emissivity Separation algorithm on ECOSTRESS data ' + filename_rad)
        try:
            (lat_eco, lon_eco, alt_eco, lf, _, vza, sza, year, month, day, hour, minute, second,
             r2, rqa2, r4, rqa4, r5, rqa5,
             date_str, hour_str, min_str, sec_str, orbit_str) = Read_L1B_Data(filename_rad, directory_geo, map_geo)
        except Exception:
            traceback.print_exc()
            print('Error while opening L1B Data -> stopping there')
            map_error[filename_rad] = traceback.format_exc()
            continue
        print('Reading ECOSTRESS L1B data completed')

        # Resize the ECOSTRESS data with a 10-by-10 window
        nsubset = 10
        nrow = np.ceil(lat_eco.shape[0] / nsubset).astype(int)
        ncol = np.ceil(lat_eco.shape[1] / nsubset).astype(int)

        lat_eco_agg = cv2.resize(lat_eco, dsize=(ncol, nrow), interpolation=cv2.INTER_LINEAR)
        lon_eco_agg = cv2.resize(lon_eco, dsize=(ncol, nrow), interpolation=cv2.INTER_LINEAR)
        alt_eco_agg = cv2.resize(alt_eco, dsize=(ncol, nrow), interpolation=cv2.INTER_LINEAR)
        lf_agg = cv2.resize(lf, dsize=(ncol, nrow), interpolation=cv2.INTER_LINEAR)
        vza_agg = cv2.resize(vza, dsize=(ncol, nrow), interpolation=cv2.INTER_LINEAR)
        sza_agg = cv2.resize(sza, dsize=(ncol, nrow), interpolation=cv2.INTER_LINEAR)
        watermask_agg = np.ones((nrow, ncol))
        watermask_agg[lf_agg > 0.5] = 0  # Land

        print('Aggregating ECOSTRESS L1B data completed')

        try:
            (p_era5, new_t, new_q, new_q2m, new_t2m, new_skt, _) = Read_ERA5(
                directory_era5, lat_eco_agg, lon_eco_agg, year, month, day, hour, minute, second)
            new_q[new_q < 0.1e-10] = 0.1e-10
        except Exception:
            traceback.print_exc()
            print('Error while opening ERA5 Data -> stopping there')
            map_error[filename_rad] = traceback.format_exc()
            continue

        print('Reading and interpolating ERA5 data completed')

        # Calculate surface pressure based on the ECOSTRESS altitude data
        sp = 1013.25 * (1.0 - 2.225577e-5 * alt_eco_agg) ** 5.25588  # Altitude unit: m
        new_sp = sp.reshape(sp.shape[0] * sp.shape[1])

        try:
            (trans0, upclear0, dnclear0) = runRTTOV(
                p_era5, new_t, new_q, new_sp, new_q2m, new_t2m,
                new_skt, vza_agg, sza_agg, watermask_agg,
                lat_eco_agg, lon_eco_agg, alt_eco_agg, year, month, day,
                rttov_installdir=rttov_installdir,
                rttov_wrapper_dir=rttov_wrapper_dir,
                rttov_coef_file=rttov_coef_file,
            )
        except Exception:
            traceback.print_exc()
            print('Error while running RTTOV -> stopping there')
            map_error[filename_rad] = traceback.format_exc()
            continue

        print('Running RTTOV completed')

        # Convert from 1-dimension to 3-dimension
        upclear = upclear0.reshape(lat_eco_agg.shape[0], lat_eco_agg.shape[1], 3)
        dnclear = dnclear0.reshape(lat_eco_agg.shape[0], lat_eco_agg.shape[1], 3)
        trans = trans0.reshape(lat_eco_agg.shape[0], lat_eco_agg.shape[1], 3)

        # Interpolate the atmospheric parameters
        upclear_f = np.zeros((lat_eco.shape[0], lat_eco.shape[1], 3))
        dnclear_f = np.zeros((lat_eco.shape[0], lat_eco.shape[1], 3))
        trans_f = np.zeros((lat_eco.shape[0], lat_eco.shape[1], 3))

        xy = np.vstack((np.ravel(lon_eco_agg), np.ravel(lat_eco_agg))).T

        for band in range(3):
            z = np.ravel(upclear[:, :, band])
            upclear_f[:, :, band] = griddata(xy, z, (lon_eco, lat_eco), method='nearest')
        print('Interpolating RTTOV upL output completed')

        for band in range(3):
            z = np.ravel(dnclear[:, :, band])
            dnclear_f[:, :, band] = griddata(xy, z, (lon_eco, lat_eco), method='nearest')
        print('Interpolating RTTOV dnL output completed')

        for band in range(3):
            z = np.ravel(trans[:, :, band])
            trans_f[:, :, band] = griddata(xy, z, (lon_eco, lat_eco), method='nearest')
        print('Interpolating RTTOV trans output completed')

        # Single-coefficient-set TES estimation (no A1/A2/A3 loop)
        try:
            (lst, emib2, emib4, emib5, _, qa) = LST_Estimate(
                r2, r4, r5, upclear_f, dnclear_f, trans_f, alpha1, alpha2, alpha3)

            bbe = np.zeros(emib2.shape)
            mask = np.logical_and.reduce((emib2 > 0, emib4 > 0, emib5 > 0))
            bbe[mask] = 0.3287 * emib2[mask] + 0.3783 * emib4[mask] + 0.3158 * emib5[mask] - 0.0255

            # Creating mask for the invalid values
            mask = np.logical_or.reduce((rqa2 > 1, rqa4 > 1, rqa5 > 1, lst == 0, bbe == 0))

            print('Estimating LST completed')
            generate_hdf5_file(directory_output, lst, emib2, emib4, emib5, bbe, mask, qa,
                                orbit_str, date_str, hour_str, min_str, sec_str)
            print('Outputting to HDF5 completed')
        except Exception:
            traceback.print_exc()
            print('Error -> stopping there')
            map_error[filename_rad] = traceback.format_exc()

    if map_error:
        pid = str(os.getpid())
        os.makedirs(error_log_dir, exist_ok=True)
        error_log_path = os.path.join(error_log_dir, f'execution_{pid}_errors.txt')
        with open(error_log_path, 'w') as file:
            file.write(json.dumps(map_error))
        print(f'Some files failed to process, see {error_log_path}')


def build_arg_parser():
    parser = argparse.ArgumentParser(
        prog='TES_main.py',
        description=('Run the Temperature Emissivity Separation (TES) algorithm '
                      'for LST estimation on ECOSTRESS L1B data (single coefficient set).'),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        '--input-files', '-i', nargs='+', metavar='RAD_FILE',
        help='One or more ECOSTRESS L1B RAD .h5 files to process.'
    )
    input_group.add_argument(
        '--input-dir', metavar='DIR',
        help='Directory to scan for ECOSTRESS L1B RAD .h5 files (used together with --pattern).'
    )
    parser.add_argument(
        '--pattern', default='*L1B_RAD*.h5',
        help='Glob pattern used to select RAD files when --input-dir is given.'
    )

    parser.add_argument('--geo-dir', required=True,
                         help='Directory containing the ECOSTRESS L1B GEO .h5 files.')
    parser.add_argument('--era5-dir', required=True,
                         help='Directory containing the ERA5 reanalysis files.')
    parser.add_argument('--output-dir', required=True,
                         help='Directory where the output LST HDF5 files are written.')

    rttov_group = parser.add_argument_group('RTTOV configuration')
    rttov_group.add_argument(
        '--rttov-installdir', required=True,
        help='RTTOV installation directory (contains e.g. rtcoef_rttov12/...).'
    )
    rttov_group.add_argument(
        '--rttov-wrapper-dir', default=None,
        help=("Directory of the RTTOV python wrapper (pyrttov), added to sys.path. "
              "Defaults to '<rttov-installdir>/wrapper' if not given.")
    )
    rttov_group.add_argument(
        '--rttov-coef-file', default=None,
        help=('Explicit path to the RTTOV coefficient file. If not given, defaults to '
              '<rttov-installdir>/rtcoef_rttov12/rttov8pred54L/rtcoef_iss_1_ecostres.dat')
    )
    rttov_group.add_argument(
        '--rttov-lib-preload', default=None,
        help=('Optional explicit path to a compiled RTTOV f2py wrapper shared library '
              '(.so) to preload via ctypes before importing pyrttov. Only needed as a '
              'workaround if pyrttov fails to load its wrapper on its own; leave unset '
              'otherwise.')
    )

    tes_group = parser.add_argument_group('TES coefficients')
    tes_group.add_argument('--alpha1', type=float, default=DEFAULT_ALPHA1,
                            help='alpha1 coefficient used in LST_Estimate().')
    tes_group.add_argument('--alpha2', type=float, default=DEFAULT_ALPHA2,
                            help='alpha2 coefficient used in LST_Estimate().')
    tes_group.add_argument('--alpha3', type=float, default=DEFAULT_ALPHA3,
                            help='alpha3 coefficient used in LST_Estimate().')

    cache_group = parser.add_argument_group('GEO file-name cache (avoids rescanning slow/networked storage)')
    cache_group.add_argument(
        '--cache-dir', default=None,
        help=('Directory used to persist the GEO file-name lookup map to disk (JSON). '
              'If not given, defaults to <output-dir>/.tes_filename_cache. '
              'Pass an empty string ("") to disable caching entirely and always rescan.')
    )
    cache_group.add_argument(
        '--force-rebuild-cache', action='store_true',
        help='Ignore any existing GEO file-name cache and rescan the directory.'
    )

    parser.add_argument('--error-log-dir', default='.',
                         help='Directory to write the execution error JSON log to, if any.')

    return parser


def main(argv=None):
    parser = build_arg_parser()
    args = parser.parse_args(argv)

    if args.input_files:
        rad_files = args.input_files
    else:
        rad_files = sorted(glob.glob(os.path.join(args.input_dir, args.pattern)))
        if not rad_files:
            parser.error(f"No files matching '{args.pattern}' found in {args.input_dir}")

    rttov_wrapper_dir = args.rttov_wrapper_dir or os.path.join(args.rttov_installdir, 'wrapper')

    os.makedirs(args.output_dir, exist_ok=True)

    if args.cache_dir is None:
        cache_dir = os.path.join(args.output_dir, '.tes_filename_cache')
    elif args.cache_dir == '':
        cache_dir = None  # caching explicitly disabled
    else:
        cache_dir = args.cache_dir

    run_TES(
        rad_files=rad_files,
        directory_geo=args.geo_dir,
        directory_era5=args.era5_dir,
        directory_output=args.output_dir,
        rttov_installdir=args.rttov_installdir,
        rttov_wrapper_dir=rttov_wrapper_dir,
        rttov_coef_file=args.rttov_coef_file,
        alpha1=args.alpha1,
        alpha2=args.alpha2,
        alpha3=args.alpha3,
        rttov_lib_preload=args.rttov_lib_preload,
        cache_dir=cache_dir,
        force_rebuild_cache=args.force_rebuild_cache,
        error_log_dir=args.error_log_dir,
    )


if __name__ == '__main__':
    main()
