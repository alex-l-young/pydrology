# =======================================================
# NEXRAD Level-II Scan processing.
# =======================================================

# Library imports.
import numpy as np
import os
from pathlib import Path
import re
import netCDF4 as nc4
from metpy.io import Level2File
import cftime
from multiprocessing import Pool
from datetime import datetime
import argparse

def nexrad_sweep_to_array(nexrad_filepath):
    """
    Extract sweep from nexrad and output at REF and RHO with x and y coordinates.
    :param nexrad_filepath: Path to NEXRAD file.
    :return: Dictionary containing data. Keys: REF_DATA, REF_X, REF_Y, RHO_DATA, RHO_X, RHO_Y.
    """

    # Load the file.
    f = Level2File(nexrad_filepath)

    # File datetime.
    file_dt = nexrad_datetime(nexrad_filepath)

    # Pull data out of the file
    sweep = 0

    # First item in ray is header, which has azimuth angle
    az = np.array([ray[0].az_angle for ray in f.sweeps[sweep]])

    # 5th item is a dict mapping a var name (byte string) to a tuple
    # of (header, data array)
    ref_hdr = f.sweeps[sweep][0][4][b'REF'][0]
    ref_range = np.arange(ref_hdr.num_gates) * ref_hdr.gate_width + ref_hdr.first_gate
    ref = np.array([ray[4][b'REF'][1] for ray in f.sweeps[sweep]])

    # Process data into array.
    ref_array = np.ma.array(ref)
    ref_array[np.isnan(ref_array)] = np.ma.masked
    ref_xlocs = ref_range * np.sin(np.deg2rad(az[:, np.newaxis]))
    ref_ylocs = ref_range * np.cos(np.deg2rad(az[:, np.newaxis]))

    # rho_hdr = f.sweeps[sweep][0][4][b'RHO'][0]
    # rho_range = (np.arange(rho_hdr.num_gates + 1) - 0.5) * rho_hdr.gate_width + rho_hdr.first_gate
    # rho = np.array([ray[4][b'RHO'][1] for ray in f.sweeps[sweep]])
    #
    # # Process data into array.
    # rho_array = np.ma.array(rho)
    # rho_array[np.isnan(rho_array)] = np.ma.masked
    # rho_xlocs = rho_range * np.sin(np.deg2rad(az[:, np.newaxis]))
    # rho_ylocs = rho_range * np.cos(np.deg2rad(az[:, np.newaxis]))

    # Data dictionary.
    nexrad_output = {
        'REF_DATA': ref_array,
        'REF_X': ref_xlocs,
        'REF_Y': ref_ylocs,
        # 'RHO_DATA': rho_array,
        # 'RHO_X': rho_xlocs,
        # 'RHO_Y': rho_ylocs,
        'DATETIME': file_dt,
    }
    print('Extracted array from {}'.format(nexrad_filepath))

    return nexrad_output


def nexrad_to_netcdf(nexrad_scans, output_directory, date, site) -> str:
    """
    Converts a list of output nexrad files to a single netcdf file.
    :param nexrad_scans: List of outputs from the nexrad_sweep_to_array function.
    :param output_directory: Directory in which to save the netcdf file.
    :param date: Day of scan files.
    :param site: NEXRAD location code.
    :return: Path to newly created netcdf file.
    """

    # Sort scans by datetime.
    nexrad_scans = sorted(nexrad_scans, key=lambda d: d['DATETIME'])

    # Load data from first nexrad file for netcdf population.
    nexdata = nexrad_scans[0]
    N_time = len(nexrad_scans)
    N_ref_lat = nexdata['REF_Y'].shape[0]
    N_ref_lon = nexdata['REF_X'].shape[1]
    # N_rho_lat = nexdata['RHO_Y'].shape[0]
    # N_rho_lon = nexdata['RHO_X'].shape[1]

    # Create netcdf shell.
    nc_fname = f'NEXRAD_{site}_{date}.nc'
    rootgrp = nc4.Dataset(output_directory / nc_fname, "w", format="NETCDF4")

    # Create dimensions in root group.
    time = rootgrp.createDimension("time", N_time)
    ref_lat = rootgrp.createDimension("refLat", N_ref_lat)
    ref_lon = rootgrp.createDimension("refLon", N_ref_lon)
    # rho_lat = rootgrp.createDimension("rhoLat", N_rho_lat)
    # rho_lon = rootgrp.createDimension("rhoLon", N_rho_lon)

    # Create variables.
    times = rootgrp.createVariable("time", "f8", ("time",))
    times.units = "hours since 0001-01-01 00:00:00.0"
    times.calendar = "gregorian"
    ref_latitudes = rootgrp.createVariable("refLat", "f4", ("refLat","refLon",))
    ref_longitudes = rootgrp.createVariable("refLon", "f4", ("refLat","refLon",))
    # rho_latitudes = rootgrp.createVariable("rhoLat", "f4", ("rhoLat", "rhoLon",))
    # rho_longitudes = rootgrp.createVariable("rhoLon", "f4", ("rhoLat", "rhoLon",))
    ref_data = rootgrp.createVariable("ref", "f4", ("time", "refLat", "refLon",))
    # rho_data = rootgrp.createVariable("rho", "f4", ("time", "rhoLat", "rhoLon",))

    nexrad_times = []
    ref_array = np.zeros((N_time, N_ref_lat, N_ref_lon))
    # rho_array = np.zeros((N_time, N_rho_lat, N_rho_lon))
    for i, scan in enumerate(nexrad_scans):
        print(f'Processing date: {datetime.strftime(scan["DATETIME"], "%Y-%m-%d %H:%M")}')

        # Data dictionary from sweep.
        try:
            ref_array[i, :, :] = scan['REF_DATA']
        except Exception as e:
            print(e)
            print(f'Could not process array for time {scan["DATETIME"]}')

        # File datetime.
        file_dt = scan["DATETIME"]
        nexrad_times.append(file_dt)


        # rho_array[i,:,:] = data['RHO_DATA']

    # Populate netcdf file.
    ref_latitudes[:] = nexdata['REF_Y']
    ref_longitudes[:] = nexdata['REF_X']
    # rho_latitudes[:] = nexdata['RHO_Y']
    # rho_longitudes[:] = nexdata['RHO_X']
    times[:] = cftime.date2num(nexrad_times,units=times.units,calendar=times.calendar)
    ref_data[:] = ref_array
    # rho_data[:] = rho_array

    # Close file.
    rootgrp.close()

    print(f'Created file: {nc_fname}')
    nc_path = output_directory / nc_fname

    return nc_path


def nexrad_datetime(filename:str) -> datetime:
    """
    Create a datetime object from NEXRAD filename.
    :param filename: NEXRAD filename.
    :return: Datetime object.
    """
    # File basename.
    basename = os.path.basename(filename)

    # Extract the datetime info from the basename.
    date_search = re.findall(r'(\d.*)_', basename)
    date_match = date_search[0]

    # Create datetime object.
    dt = datetime.strptime(date_match, "%Y%m%d_%H%M%S")

    return dt

def parse_arguments():
    # Parse the arguments.
    parser = argparse.ArgumentParser(description='Process NEXRAD Data from a single day.')
    parser.add_argument('--dir', help='Directory containing NEXRAD file day sub-directories.')
    parser.add_argument('--pool', help='Number of processors to pool.')
    parser.add_argument('--date', help='Date of the scans yyyy-mm-dd')
    parser.add_argument('--site', help='NEXRAD site code.')
    args = parser.parse_args()

    Npool = int(args.pool)
    date = args.date
    site = args.site
    dir = Path(args.dir) / date

    # Get the files to be processed.
    ld = dir.glob('**/*')
    output_files = [dir / f.name for f in ld if str(f)[-3:] == 'V06']

    parsed_args = {
        'Npool': Npool,
        'date': date,
        'site': site,
        'output_files': output_files,
        'dir': dir,
    }

    return parsed_args


if __name__ == "__main__":
    # Parse arguments.
    parsed_args = parse_arguments()
    Npool = parsed_args['Npool']
    date = parsed_args['date']
    site = parsed_args['site']
    output_files = parsed_args['output_files']
    dir = parsed_args['dir']

    # Multiprocess scan array extraction.
    with Pool(Npool) as p:
        nexrad_scans = p.map(nexrad_sweep_to_array, output_files)

    # Create netCDF file from scans.
    nexrad_to_netcdf(nexrad_scans, dir, date, site)