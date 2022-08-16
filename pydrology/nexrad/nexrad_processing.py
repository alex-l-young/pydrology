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
from itertools import repeat
from datetime import datetime
import argparse
import yaml

# Local imports.
from nexrad_config import NexradConfig

def nexrad_sweep_to_array(nexrad_filepath, cfg):
    """
    Extract sweep from nexrad and output at REF and RHO with x and y coordinates.
    :param nexrad_filepath: Path to NEXRAD file.
    :param sweep: Sweep level to extract.
    :return: Dictionary containing data. Keys: REF_DATA, REF_X, REF_Y, RHO_DATA, RHO_X, RHO_Y.
    """

    # Load the file.
    f = Level2File(nexrad_filepath)

    # Sweep level.
    sweep = cfg.sweep

    # File datetime.
    file_dt = nexrad_datetime(nexrad_filepath)

    # Data dictionary.
    nexrad_output = {
        'DATETIME': file_dt,
    }

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
    if cfg.variable in ['REF', 'ALL']:
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

        if cfg.grid_interp is True:
            # Create interpolation grid.
            glat, glon = NexNCDF.create_lat_lon_grid(cfg.bounds, cfg.spacing, cfg.num_points)

            # Interpolate scan to lat/lon.
            ref_array = NexNCDF.grid_interpolate_to_lat_lon(ref_array, ref_ylocs_grid, ref_xlocs_grid, glat, glon)

            # Update lat/lon grid.
            ref_xlocs_grid = glon
            ref_ylocs_grid = glat

        # Populate dictionary.
        nexrad_output['REF_DATA'] = ref_array
        nexrad_output['REF_X'] = ref_xlocs_grid
        nexrad_output['REF_Y'] = ref_ylocs_grid


    if cfg.variable in ['RHO', 'ALL']:
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

        if cfg.grid_interp is True:
            # Create interpolation grid.
            glat, glon = NexNCDF.create_lat_lon_grid(cfg.bounds, cfg.spacing, cfg.num_points)

            # Interpolate scan to lat/lon.
            rho_array = NexNCDF.grid_interpolate_to_lat_lon(rho_array, rho_ylocs_grid, rho_xlocs_grid, glat, glon)

            # Update lat/lon grid.
            rho_xlocs_grid = glon
            rho_ylocs_grid = glat

        # Populate dictionary.
        nexrad_output['RHO_DATA'] = rho_array
        nexrad_output['RHO_X'] = rho_xlocs_grid
        nexrad_output['RHO_Y'] = rho_ylocs_grid

    print('Extracted array from {}'.format(nexrad_filepath))

    return nexrad_output


class NexNCDF():
    def __init__(self, cfg:NexradConfig, nexrad_scans:list):
        """
        Nexrad NETCDF Class.
        Creates an instance of a netcdf file and populates it with nexrad scan data.
        :param cfg: NEXRAD configuration object.
        :param nexrad_scans: Nexrad scan data in a dictionary.
        """
        # Add configuration to class instance.
        self.cfg = cfg

        # Extract configuration variables.
        self.output_directory = self.cfg.nexrad_directory # Output directory for netCDF file.
        self.date = self.cfg.date # Date of the scans.
        self.site = self.cfg.site # Site code of scans.

        # Sort scans by datetime.
        self.nexrad_scans = sorted(nexrad_scans, key=lambda d: d['DATETIME'])

        # Instantiate a netcdf file and open root group.
        self.nc_fname = f'NEXRAD_{self.cfg.site}_{self.cfg.date}.nc'
        self.nc_path = self.output_directory / self.nc_fname
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
        if cfg.variable in ['REF', 'ALL']:
            self.N_ref_lat = nexdata['REF_Y'].shape[0]
            self.N_ref_lon = nexdata['REF_X'].shape[1]

        # Rho dimensions.
        if cfg.variable in ['RHO', 'ALL']:
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
        if self.cfg.grid_interp is True:
            # Single set of lat and lon if interpolating to grid.
            self.ref_latitudes = self.rootgrp.createVariable("refLat", "f4", ("refLat", "refLon",))
            self.ref_longitudes = self.rootgrp.createVariable("refLon", "f4", ("refLat", "refLon",))
            coord_dims = (self.N_ref_lat, self.N_ref_lon) # Dimensions of the lat/lon coordinates.
        else:
            # Separate lat/lon for each time step if not interpolating to grid.
            self.ref_latitudes = self.rootgrp.createVariable("refLat", "f4", ("time", "refLat", "refLon",))
            self.ref_longitudes = self.rootgrp.createVariable("refLon", "f4", ("time", "refLat", "refLon",))
            coord_dims = (self.N_time, self.N_ref_lat, self.N_ref_lon) # Dimensions of the lat/lon coordinates.

        # Data variable.
        self.ref_data = self.rootgrp.createVariable("ref", "f4", ("time", "refLat", "refLon",))

        # Dimensions of the data.
        ref_dims = {"Data": (self.N_time, self.N_ref_lat, self.N_ref_lon), "Coords": coord_dims}

        return ref_dims


    def _setup_rho(self):

        # Coordinate dimensions.
        rho_lat = self.rootgrp.createDimension("rhoLat", self.N_rho_lat)
        rho_lon = self.rootgrp.createDimension("rhoLon", self.N_rho_lon)

        # Coordinate variables.
        if self.cfg.grid_interp is True:
            # Single set of lat and lon if interpolating to grid.
            self.rho_latitudes = self.rootgrp.createVariable("rhoLat", "f4", ("rhoLat", "rhoLon",))
            self.rho_longitudes = self.rootgrp.createVariable("rhoLon", "f4", ("rhoLat", "rhoLon",))
            coord_dims = (self.N_rho_lat, self.N_rho_lon)  # Dimensions of the lat/lon coordinates.
        else:
            # Separate lat/lon for each time step if not interpolating to grid.
            self.rho_latitudes = self.rootgrp.createVariable("rhoLat", "f4", ("time", "rhoLat", "rhoLon",))
            self.rho_longitudes = self.rootgrp.createVariable("rhoLon", "f4", ("time", "rhoLat", "rhoLon",))
            coord_dims = (self.N_time, self.N_rho_lat, self.N_rho_lon)  # Dimensions of the lat/lon coordinates.

        # Data variable.
        self.rho_data = self.rootgrp.createVariable("rho", "f4", ("time", "rhoLat", "rhoLon",))

        # Dimensions of the data.
        rho_dims = {"Data": (self.N_time, self.N_rho_lat, self.N_rho_lon), "Coords": coord_dims}

        return rho_dims


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
        ref_dims = self._setup_ref()

        # REF data array instantiation.
        ref_array = np.zeros(ref_dims["Data"])
        ref_xlocs_array = np.zeros(ref_dims["Coords"])
        ref_ylocs_array = np.zeros(ref_dims["Coords"])

        # Save interpolation grid in NetCDF arrays.
        if self.cfg.grid_interp is True:
            ref_xlocs_array = self.nexrad_scans[0]['REF_X']
            ref_ylocs_array = self.nexrad_scans[0]['REF_Y']

        # Loop through the scans and populate the data arrays.
        print('Adding REF data...')
        for i, scan in enumerate(self.nexrad_scans):
            print(f'Processing date: {datetime.strftime(scan["DATETIME"], "%Y-%m-%d %H:%M")}')
            # Store data in arrays.
            try:
                if self.cfg.grid_interp is False:
                    # Save individual coordinate arrays.
                    ref_xlocs_array[i, :, :] = scan['REF_X']
                    ref_ylocs_array[i, :, :] = scan['REF_Y']

                # Save reflectivity data in array.
                ref_array[i, :, :] = scan['REF_DATA']

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
        rho_dims = self._setup_rho()

        # RHO data array instantiation.
        rho_array = np.zeros(rho_dims["Data"])
        rho_xlocs_array = np.zeros(rho_dims["Coords"])
        rho_ylocs_array = np.zeros(rho_dims["Coords"])

        # Save interpolation grid in NetCDF arrays.
        if self.cfg.grid_interp is True:
            rho_xlocs_array = self.nexrad_scans[0]['RHO_X']
            rho_ylocs_array = self.nexrad_scans[0]['RHO_Y']

        # Loop through the scans and populate the data arrays.
        print('Adding RHO data...')
        for i, scan in enumerate(self.nexrad_scans):
            print(f'Processing date: {datetime.strftime(scan["DATETIME"], "%Y-%m-%d %H:%M")}')
            # Store data in arrays.
            try:
                if self.cfg.grid_interp is False:
                    # Save individual coordinate arrays.
                    rho_xlocs_array[i, :, :] = scan['RHO_X']
                    rho_ylocs_array[i, :, :] = scan['RHO_Y']

                # Save reflectivity data in array.
                rho_array[i, :, :] = scan['RHO_DATA']
            except Exception as e:
                print(e)
                print(f'Could not process array for time {scan["DATETIME"]}')

        # Populate netcdf file.
        self.rho_data[:] = rho_array
        self.rho_latitudes[:] = rho_ylocs_array
        self.rho_longitudes[:] = rho_xlocs_array

        print(f'Added RHO data to {self.nc_fname}')


    @staticmethod
    def create_lat_lon_grid(bounds:tuple, spacing:tuple=None, npoints:tuple=None) -> tuple:
        """
        Creates an evenly spaced lat/lon grid using with specified spacing or number of points.
        Must specify either spacing or npoints.
        :param bounds: Grid bounds (N, E, S, W) in units of the current lat/lon.
        :param spacing: Grid spacing to use.
        :param npoints: Number of points to use in the grid.
        :return: Creates an evenly spaced grid of points to interpolate to. glat, glon
        """
        # Expand bounds.
        N, E, S, W = bounds

        # Create 1D vectors of lat and lon.
        if spacing is not None:
            lat_spacing = spacing[1]
            lon_spacing = spacing[0]
            lat_vec = np.arange(S, N, lat_spacing)
            lon_vec = np.arange(W, E, lon_spacing)
            # glon, glat = np.mgrid[W:E:lon_spacing, S:N:lat_spacing]
        else:
            lat_points = npoints[1]
            lon_points = npoints[0]
            lat_vec = np.linspace(S, N, lat_points)
            lon_vec = np.linspace(W, E, lon_points)
            # glon, glat = np.mgrid[W:E:lon_points*1j, S:N:lat_points*1j]

        # Make meshgrid.
        glon, glat = np.meshgrid(lon_vec, lat_vec)

        return glat, glon


    @staticmethod
    def grid_interpolate_to_lat_lon(scan:np.array, lat:np.array, lon:np.array, glat:np.array, glon:np.array) -> np.array:
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
        scan_reshape = np.reshape(scan, (-1,))
        lat_reshape = np.reshape(lat, (-1,1))
        lon_reshape = np.reshape(lon, (-1,1))
        latlon = np.concatenate((lon_reshape, lat_reshape), axis=1)

        # Grid interpolation.
        gscan = griddata(latlon, scan_reshape, (glon, glat), method='nearest')

        return gscan


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
    parser.add_argument('--config', help='Path to NEXRAD processing configuration file.')
    args = parser.parse_args()

    # Create path from configuration file path.
    cfg_filepath = args.config

    # Create configuration instance.
    cfg = NexradConfig(cfg_filepath)

    return cfg


if __name__ == "__main__":
    # Parse arguments.
    cfg = parse_arguments()

    # Multiprocess scan array extraction.
    with Pool(cfg.npool) as p:
        nexrad_scans = p.starmap(nexrad_sweep_to_array, zip(cfg.output_files, repeat(cfg)))

    # Create netCDF file from scans.
    nexrad_netcdf = NexNCDF(cfg, nexrad_scans)
    if cfg.variable in ['REF', 'ALL']:
        nexrad_netcdf.add_ref_data()
    if cfg.variable in ['RHO', 'ALL']:
        nexrad_netcdf.add_rho_data()
    nexrad_netcdf.close_file()