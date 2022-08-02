# =======================================================
# Data request functions.
# =======================================================

# Library imports.
import pandas as pd
import numpy as np
import requests
import os
from pathlib import Path
from io import StringIO
import netCDF4 as nc4
import nexradaws
import six
from metpy.io import Level2File
from datetime import datetime
import cftime
import re
import threading

# =======================================================
# NOAA SC-ACIS Meteorology.
# =======================================================
def request_acis_data(met_elements, site_id, start_date, end_date):
    """
    Request meteorology data from NOAA's ACIS. Data is daily and the following meteorology variables can be requested.
        - Maximum Temperature (maxt)
        - Minimum Temperature (mint)
        - Average Temperature (avgt)
        - Precipitation (pcpn)
        - Snow Depth (snow)
    :param met_elements [list]: Names of the elements/variables to be requested.
    :param site_id [string]: Site ID name. Add <space>6 to the end of GHCN site IDs. Add <space>2 for coop site ids.
        (E.g., "USS0020A23S" => "USS0020A23S 6")
    :param start_date [string]: Start date of request. Format yyyy-mm-dd. E.g., '2022-01-01'
    :param end_date [string]:  End date of request. Format yyyy-mm-dd. E.g., '2022-01-01'
    :return: Data Frame of the requested ACIS data.
    """

    # Meteorology elements reference.
    met_elements_ref = {'maxt': 'MaxTemp', 'mint': 'MinTemp', 'avgt': 'AvgTemp', 'pcpn': 'Precipitation',
                        'snow': 'SnowDepth'}

    # Populate json request.
    json_req = {"elems":[]}

    # Add in met elements.
    # Also build final dictionary.
    for elem in met_elements:
        json_req["elems"].append({"name": elem})

    # Add date range.
    json_req["sDate"] = start_date
    json_req["eDate"] = end_date

    # Add site ID.
    json_req['sid'] = site_id

    # Make request.
    url = r'https://data.rcc-acis.org/StnData'
    r = requests.post(url, json=json_req)

    # Get dictionary response.
    try:
        met_data = r.json()
        # print(met_data)
    except Exception as e:
        print(e)
        print('Could not make MET Dictionary from response.')
        print(r.text)

    # Check if all values are present. Number of met elements in request plus one date value.
    if len(met_data['data'][0]) != len(met_elements) + 1:
        print('Met data length does not agree with number of elements requested.')
        print(met_data)

    # Extract the met data into a data dictionary.
    extract_data = {}
    extract_data['Date'] = []
    for elem in met_elements:
        # Get parameter name from met element abbr.
        param_name = met_elements_ref[elem]
        extract_data[param_name] = []

    for entry in met_data['data']:
        extract_data['Date'].append(entry[0])
        for i in range(1, len(entry)):
            param_name = met_elements_ref[met_elements[i-1]] # Get the full name of the met element.
            extract_data[param_name].append(entry[i])

    # Create data frame.
    extract_data['site_id'] = np.repeat(site_id.replace(' ', '-'), len(extract_data['Date']))
    extract_data['name'] = np.repeat(met_data['meta']['name'], len(extract_data['Date']))
    extract_data['latitude'] = np.repeat(met_data['meta']['ll'][1], len(extract_data['Date']))
    extract_data['longitude'] = np.repeat(met_data['meta']['ll'][0], len(extract_data['Date']))
    acis_df = pd.DataFrame(extract_data)

    return acis_df


def acis_to_csv(acis_df, csv_dir):
    """
    Save the ACIS DataFrame as a csv.
    :param acis_df [Data Frame]: ACIS DataFrame. Constructed from request_acis_data().
    :param csv_dir [string]: Directory in which to save the ACIS data.
    :return: None
    """

    # Save to csv file "<site_name>_<site_id>_<lat>_<lon>_<start_date>_<end_date>.csv".
    site_name = acis_df.loc[0, 'name']
    site_id = acis_df.loc[0, 'site_id']
    lat = acis_df.loc[0, 'latitude']
    lon = acis_df.loc[0, 'longitude']
    start_date = acis_df.loc[0, 'Date']
    end_date = acis_df.loc[acis_df.index[-1], 'Date']
    df_fname = f'{site_name}_{site_id}_{lat}_{lon}_{start_date}_{end_date}.csv'

    acis_df.to_csv(os.path.join(csv_dir, df_fname), index=False)


# =======================================================
# USGS Stream Gage Data.
# =======================================================
def request_usgs_data(gage_id, parameter, start_date, start_time, end_date, end_time, gmt_offset):
    """

    :param gage_id [string]: ID of gage to request data for.
    :param start_date [string]: Start date in format yyyy-mm-dd. "2022-06-24"
    :param start_time [string]: Local start time in format HH:MM:SS.mmm. "11:17:05.203"
    :param end_date [string]: End date in format yyyy-mm-dd.
    :param end_time [string]: Local end time in format HH:MM:SS.mmm
    :param gmt_offset [string]: Number of hour offset from gmt (+ or -) in format +/-HH:MM. "-04:00"
    :param parameter [string]: 'discharge' or 'height'
    :return:
    """
    # Convert parameter name to the code.
    if parameter == 'discharge':
        code = '00060'
    else:
        code = '00065'

    # Populate the request URL with ID and dates.
    usgs_url_format = ("https://waterservices.usgs.gov/nwis/iv/?sites={}&parameterCd={}&"
                       "startDT={}T{}{}&endDT={}T{}{}&"
                       "siteStatus=all&format=rdb")
    usgs_url = usgs_url_format.format(gage_id, code, start_date, start_time, gmt_offset, end_date, end_time, gmt_offset)

    r = requests.get(usgs_url)
    r_text = r.text

    gage_df = create_discharge_df(r_text, parameter)

    return gage_df


def create_discharge_df(usgs_text: str, parameter: str) -> pd.DataFrame:
    cols = ['agency', 'site_no', 'datetime', 'tz_cd', parameter, 'provis_accept']
    skiprows = 25

    # Keep incrementing the number of rows to skip until the correct number is reached.
    while True:
        try:
            df = pd.read_csv(StringIO(usgs_text), header=None, sep='\t', skiprows=skiprows)
            df.columns = cols

            # There are two rows of headers that also need to be skipped.
            if df.loc[0, 'agency'] != 'USGS':
                skiprows += 1
                continue

            break
        except:
            skiprows += 1
            pass

    return df

# =======================================================
# NEXRAD Doppler Data.
# =======================================================
def request_nexrad_data(site, date, output_directory):
    """

    :param site: NEXRAD call code. E.g., "KBGM"
    :param date: Date of request. "YYYY-MM-DD"
    :param output_directory: Directory in which to output the scans.
    :return:
    """
    # Split the date into year, month, and day.
    date_split = date.split('-')
    year, month, day = (date_split[0], date_split[1], date_split[2])

    conn = nexradaws.NexradAwsInterface()
    scans = conn.get_avail_scans(year, month, day, site)
    print(f'Found {len(scans)} scans for this date...', '\n')

    # Sort scans by date.
    times = []
    scans_to_download = []
    for scan in scans:
        scan_str = str(scan)
        scan_split = scan_str.split('_')

        # Don't process MDM files.
        if 'MDM' not in scan_split[-1]:
            time = scan_split[-2]
            times.append(int(time))
            scans_to_download.append(scan)
        else:
            continue

    scan_and_date = zip(scans_to_download,times)
    scan_sort = [item[0] for item in sorted(scan_and_date, key=lambda x: x[1])]

    # Save all nexrad scans.
    output_files = []
    c = 1
    for scan in scan_sort:
        print(f'SCAN {c}')
        c += 1
        file_name = str(scan).split('/')[-1][:-1]
        output_file = os.path.join(output_directory, file_name)
        output_files.append(os.path.join(output_file, file_name))

        localfiles = conn.download(scan, output_file)
        six.print_(localfiles.success)
        six.print_(localfiles.success[0].filepath)
        # create_download_thread(scan, output_file)

    return output_files


def download_nexrad(scan, output_file):
    conn = nexradaws.NexradAwsInterface()
    localfiles = conn.download(scan, output_file)
    six.print_(localfiles.success)
    six.print_(localfiles.success[0].filepath)

    return localfiles


def create_download_thread(scan, output_file):
    download_thread = threading.Thread(target=download_nexrad, args=(scan, output_file))
    download_thread.start()

def nexrad_sweep_to_array(nexrad_filepath):
    """
    Extract sweep from nexrad and output at REF and RHO with x and y coordinates.
    :param nexrad_filepath: Path to NEXRAD file.
    :return: Dictionary containing data. Keys: REF_DATA, REF_X, REF_Y, RHO_DATA, RHO_X, RHO_Y.
    """

    # Load the file.
    f = Level2File(nexrad_filepath)

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
    }

    return nexrad_output


def nexrad_to_netcdf(nexrad_scan_files, output_directory, date, site) -> str:
    """
    Converts a list of output nexrad files to a single netcdf file.
    :param nexrad_scan_files: List of file paths to add to the netcdf file.
    :param output_directory: Directory in which to save the netcdf file.
    :param date: Day of scan files.
    :param site: NEXRAD location code.
    :return: Path to newly created netcdf file.
    """

    # Load data from first nexrad file for netcdf population.
    nexdata = nexrad_sweep_to_array(nexrad_scan_files[0])
    N_time = len(nexrad_scan_files)
    N_ref_lat = nexdata['REF_Y'].shape[0]
    N_ref_lon = nexdata['REF_X'].shape[1]
    # N_rho_lat = nexdata['RHO_Y'].shape[0]
    # N_rho_lon = nexdata['RHO_X'].shape[1]

    # Create netcdf shell.
    nc_fname = f'NEXRAD_{site}_{date}.nc'
    rootgrp = nc4.Dataset(os.path.join(output_directory, nc_fname), "w", format="NETCDF4")

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
    for i, file in enumerate(nexrad_scan_files):
        print(f'Processing file: {file}')

        # File datetime.
        file_dt = nexrad_datetime(file)
        nexrad_times.append(file_dt)

        # Data dictionary from sweep.
        data = nexrad_sweep_to_array(file)
        ref_array[i,:,:] = data['REF_DATA']
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


if __name__ == '__main__':
    site = 'KBGM'
    date = '2022-07-25'
    output_directory = Path.home() / 'Desktop' / 'NEXRAD_Output'
    output_files = request_nexrad_data(site, date, output_directory)
    print(output_files)
    # output_files = ['/Users/alexyoung/Desktop/NEXRAD_Output/KBGM20220725_000141_V06/KBGM20220725_000141_V06', '/Users/alexyoung/Desktop/NEXRAD_Output/KBGM20220725_000711_V06/KBGM20220725_000711_V06', '/Users/alexyoung/Desktop/NEXRAD_Output/KBGM20220725_001241_V06/KBGM20220725_001241_V06', '/Users/alexyoung/Desktop/NEXRAD_Output/KBGM20220725_001812_V06/KBGM20220725_001812_V06', '/Users/alexyoung/Desktop/NEXRAD_Output/KBGM20220725_002342_V06/KBGM20220725_002342_V06', '/Users/alexyoung/Desktop/NEXRAD_Output/KBGM20220725_002914_V06/KBGM20220725_002914_V06', '/Users/alexyoung/Desktop/NEXRAD_Output/KBGM20220725_003458_V06/KBGM20220725_003458_V06', '/Users/alexyoung/Desktop/NEXRAD_Output/KBGM20220725_004042_V06/KBGM20220725_004042_V06', '/Users/alexyoung/Desktop/NEXRAD_Output/KBGM20220725_004639_V06/KBGM20220725_004639_V06', '/Users/alexyoung/Desktop/NEXRAD_Output/KBGM20220725_005224_V06/KBGM20220725_005224_V06']
    #
    #
    nexrad_to_netcdf(output_files, output_directory, date, site)