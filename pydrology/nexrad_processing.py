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
from metpy.units import units
from metpy.calc.tools import azimuth_range_to_lat_lon
from scipy.interpolate import griddata
import cftime
from multiprocessing import Pool
from datetime import datetime
import argparse

def nexrad_sweep_to_array(nexrad_filepath, sweep=0):
    """
    Extract sweep from nexrad and output at REF and RHO with x and y coordinates.
    :param nexrad_filepath: Path to NEXRAD file.
    :param sweep: Sweep level to extract.
    :return: Dictionary containing data. Keys: REF_DATA, REF_X, REF_Y, RHO_DATA, RHO_X, RHO_Y.
    """

    # Load the file.
    f = Level2File(nexrad_filepath)

    # File datetime.
    file_dt = nexrad_datetime(nexrad_filepath)

    # First item in ray is header, which has azimuth angle
    az = np.array([ray[0].az_angle for ray in f.sweeps[sweep]])
    diff = np.diff(az)
    diff[diff > 180] -= 360.
    diff[diff < -180] += 360.
    avg_spacing = diff.mean()
    az = (az[:-1] + az[1:]) / 2
    az = np.concatenate(([az[0] - avg_spacing], az, [az[-1] + avg_spacing]))
    az = units.Quantity(az, 'degrees')

    # Center lat and lon.
    center_lat = f.sweeps[sweep][0][1].lat
    center_lon = f.sweeps[sweep][0][1].lon

    # 5th item is a dict mapping a var name (byte string) to a tuple
    # of (header, data array)
    ref_hdr = f.sweeps[sweep][0][4][b'REF'][0]
    ref_range = (np.arange(ref_hdr.num_gates + 1) - 0.5) * ref_hdr.gate_width + ref_hdr.first_gate
    ref_range = units.Quantity(ref_range, 'kilometers')
    ref = np.array([ray[4][b'REF'][1] for ray in f.sweeps[sweep]])

    # Process data into array.
    ref_array = np.ma.array(ref)
    ref_array[np.isnan(ref_array)] = np.ma.masked
    ref_xlocs, ref_ylocs = azimuth_range_to_lat_lon(az, ref_range, center_lon, center_lat)

    # Center of ref lat and lon cells.
    ref_xlocs_grid = (ref_xlocs[1:, 1:] + ref_xlocs[:-1, :-1]) / 2
    ref_ylocs_grid = (ref_ylocs[1:, 1:] + ref_ylocs[:-1, :-1]) / 2

    # Extract RHO.
    rho_hdr = f.sweeps[sweep][0][4][b'RHO'][0]
    rho_range = (np.arange(rho_hdr.num_gates + 1) - 0.5) * rho_hdr.gate_width + rho_hdr.first_gate
    rho_range = units.Quantity(rho_range, 'kilometers')
    rho = np.array([ray[4][b'RHO'][1] for ray in f.sweeps[sweep]])

    # Process data into array.
    rho_array = np.ma.array(rho)
    rho_array[np.isnan(rho_array)] = np.ma.masked
    rho_xlocs, rho_ylocs = azimuth_range_to_lat_lon(az, rho_range, center_lon, center_lat)

    # Center of rho lat and lon cells.
    rho_xlocs_grid = (rho_xlocs[1:, 1:] + rho_xlocs[:-1, :-1]) / 2
    rho_ylocs_grid = (rho_ylocs[1:, 1:] + rho_ylocs[:-1, :-1]) / 2

    # Data dictionary.
    nexrad_output = {
        'REF_DATA': ref_array,
        'REF_X': ref_xlocs_grid,
        'REF_Y': ref_ylocs_grid,
        'RHO_DATA': rho_array,
        'RHO_X': rho_xlocs_grid,
        'RHO_Y': rho_ylocs_grid,
        'DATETIME': file_dt,
    }
    print('Extracted array from {}'.format(nexrad_filepath))

    return nexrad_output

class NexNCDF():
    def __init__(self, nexrad_scans:list, output_directory:Path, date:str, site:str):
        """
        Nexrad NETCDF Class.
        Creates an instance of a netcdf file and populates it with nexrad scan data.
        :param nexrad_scans: List of nexrad scan dictionaries generated by nexrad_sweep_to_array.
        :param output_directory: Directory to save the netCDF file.
        :param date: Scan date. 'YYYY-MM-DD'
        :param site: Site code.
        """

        self.output_directory = output_directory
        self.date = date
        self.site = site

        # Sort scans by datetime.
        self.nexrad_scans = sorted(nexrad_scans, key=lambda d: d['DATETIME'])

        # Instantiate a netcdf file and open root group.
        self.nc_fname = f'NEXRAD_{site}_{date}.nc'
        self.nc_path = output_directory / self.nc_fname
        self.rootgrp = nc4.Dataset(self.nc_path, "w", format="NETCDF4")

        # NetCDF dimensions.
        self._scan_dimensions()

        # Set up and populate the time variable.
        self._setup_time()
        self._add_time_data()


    def __del__(self):
        if self.rootgrp._isopen:
            print('...Closing Net CDF File.')
            self.rootgrp.close()


    def close_file(self):
        if self.rootgrp._isopen:
            print('...Closing Net CDF File.')
            self.rootgrp.close()


    def _scan_dimensions(self):
        # Number of time steps.
        self.N_time = len(self.nexrad_scans)

        # Use first scan for dimensions.
        nexdata = self.nexrad_scans[0]

        # Reflectivity dimensions.
        self.N_ref_lat = nexdata['REF_Y'].shape[0]
        self.N_ref_lon = nexdata['REF_X'].shape[1]

        # Rho dimensions.
        self.N_rho_lat = nexdata['RHO_Y'].shape[0]
        self.N_rho_lon = nexdata['RHO_X'].shape[1]


    def _setup_time(self):

        # Time dimension.
        time = self.rootgrp.createDimension("time", self.N_time)

        # Time variable.
        self.times = self.rootgrp.createVariable("time", "f8", ("time",))
        self.times.units = "hours since 0001-01-01 00:00:00.0"
        self.times.calendar = "gregorian"


    def _setup_ref(self):

        # Coordinate dimensions.
        ref_lat = self.rootgrp.createDimension("refLat", self.N_ref_lat)
        ref_lon = self.rootgrp.createDimension("refLon", self.N_ref_lon)

        # Coordinate variables.
        self.ref_latitudes = self.rootgrp.createVariable("refLat", "f4", ("time", "refLat", "refLon",))
        self.ref_longitudes = self.rootgrp.createVariable("refLon", "f4", ("time", "refLat", "refLon",))

        # Data variable.
        self.ref_data = self.rootgrp.createVariable("ref", "f4", ("time", "refLat", "refLon",))


    def _setup_rho(self):

        # Coordinate dimensions.
        rho_lat = self.rootgrp.createDimension("rhoLat", self.N_rho_lat)
        rho_lon = self.rootgrp.createDimension("rhoLon", self.N_rho_lon)

        # Coordinate variables.
        self.rho_latitudes = self.rootgrp.createVariable("rhoLat", "f4", ("time", "rhoLat", "rhoLon",))
        self.rho_longitudes = self.rootgrp.createVariable("rhoLon", "f4", ("time", "rhoLat", "rhoLon",))

        # Data variable.
        self.rho_data = self.rootgrp.createVariable("rho", "f4", ("time", "rhoLat", "rhoLon",))


    def _add_time_data(self):
        print('Adding time data.')
        # Time array instantiation.
        nexrad_times = []

        # Loop through the scans and populate the data arrays.
        for i, scan in enumerate(self.nexrad_scans):
            print(f'Processing date: {datetime.strftime(scan["DATETIME"], "%Y-%m-%d %H:%M")}')
            # Extract data from the scan dictionary.
            try:
                # File datetime.
                file_dt = scan["DATETIME"]
                nexrad_times.append(file_dt)
            except Exception as e:
                print(e)
                print(f'Could not process array for time {scan["DATETIME"]}')

        # Populate netcdf file.
        self.times[:] = cftime.date2num(nexrad_times, units=self.times.units, calendar=self.times.calendar)


    def add_ref_data(self):
        # Set up the ref dimensions and variables.
        self._setup_ref()

        # REF data array instantiation.
        ref_dims = (self.N_time, self.N_ref_lat, self.N_ref_lon)
        ref_array = np.zeros(ref_dims)
        ref_xlocs_array = np.zeros(ref_dims)
        ref_ylocs_array = np.zeros(ref_dims)

        # Loop through the scans and populate the data arrays.
        print('Adding REF data...')
        for i, scan in enumerate(self.nexrad_scans):
            print(f'Processing date: {datetime.strftime(scan["DATETIME"], "%Y-%m-%d %H:%M")}')
            # Extract data from the scan dictionary.
            try:
                # REF.
                ref_array[i, :, :] = scan['REF_DATA']
                ref_xlocs_array[i, :, :] = scan['REF_X']
                ref_ylocs_array[i, :, :] = scan['REF_Y']
            except Exception as e:
                print(e)
                print(f'Could not process array for time {scan["DATETIME"]}')

        # Populate netcdf file.
        self.ref_data[:] = ref_array
        self.ref_latitudes[:] = ref_ylocs_array
        self.ref_longitudes[:] = ref_xlocs_array

        print(f'Added REF data to {self.nc_fname}')


    def add_rho_data(self):
        # Set up the rho dimensions and variables.
        self._setup_rho()

        # RHO data array instantiation.
        rho_dims = (self.N_time, self.N_rho_lat, self.N_rho_lon)
        rho_array = np.zeros(rho_dims)
        rho_xlocs_array = np.zeros(rho_dims)
        rho_ylocs_array = np.zeros(rho_dims)

        # Loop through the scans and populate the data arrays.
        print('Adding RHO data...')
        for i, scan in enumerate(self.nexrad_scans):
            print(f'Processing date: {datetime.strftime(scan["DATETIME"], "%Y-%m-%d %H:%M")}')
            # Extract data from the scan dictionary.
            try:
                # RHO.
                rho_array[i, :, :] = scan['RHO_DATA']
                rho_xlocs_array[i, :, :] = scan['RHO_X']
                rho_ylocs_array[i, :, :] = scan['RHO_Y']
            except Exception as e:
                print(e)
                print(f'Could not process array for time {scan["DATETIME"]}')

        # Populate netcdf file.
        self.rho_data[:] = rho_array
        self.rho_latitudes[:] = rho_ylocs_array
        self.rho_longitudes[:] = rho_xlocs_array

        print(f'Added RHO data to {self.nc_fname}')


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


def grid_interpolate_to_lat_lon(scan, lat, lon, glat, glon):
    """
    Interpolates the scan to a specified lat/lon grid.
    Uses linear interpolation.
    :param scan: The scan to interpolate.
    :param lat: Latitude grid corresponding to the scan.
    :param lon: Longitude grid corresponding to the scan.
    :param glat: Latitude grid to interpolate to.
    :param glon: Longitude grid to interpolate to.
    :return: Interpolated scan.
    """
    # Reshape to column vectors.
    scan_reshape = np.reshape(scan, (-1, 1))
    lat_reshape = np.reshape(lat, (-1, 1))
    lon_reshape = np.reshape(lon, (-1, 1))
    latlon = np.concatenate((lon_reshape, lat_reshape), axis=1)

    # Grid interpolation.
    gscan = griddata(scan_reshape, latlon, (glon, glat), method='linear')

    return gscan


def parse_arguments():
    # Parse the arguments.
    parser = argparse.ArgumentParser(description='Process NEXRAD Data from a single day.')
    parser.add_argument('--dir', help='Directory containing NEXRAD file day sub-directories.')
    parser.add_argument('--pool', help='Number of processors to pool.')
    parser.add_argument('--date', help='Date of the scans yyyy-mm-dd')
    parser.add_argument('--site', help='NEXRAD site code.')
    parser.add_argument('--var', help='Radar variables to include. One of: ("REF", "RHO", "ALL")')
    args = parser.parse_args()

    Npool = int(args.pool)
    date = args.date
    site = args.site
    dir = Path(args.dir) / date
    var = args.var

    # Get the files to be processed.
    ld = dir.glob('**/*')
    output_files = [dir / f.name for f in ld if str(f)[-3:] == 'V06']

    parsed_args = {
        'Npool': Npool,
        'date': date,
        'site': site,
        'output_files': output_files,
        'dir': dir,
        'var': var
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
    var = parsed_args['var']

    # Multiprocess scan array extraction.
    with Pool(Npool) as p:
        nexrad_scans = p.map(nexrad_sweep_to_array, output_files)

    # Create netCDF file from scans.
    nexrad_netcdf = NexNCDF(nexrad_scans, dir, date, site)
    if var in ['REF', 'ALL']:
        nexrad_netcdf.add_ref_data()
    if var in ['RHO', 'ALL']:
        nexrad_netcdf.add_rho_data()
    nexrad_netcdf.close_file()