# =======================================================
# USGS Stream Gage Data.
# =======================================================

# Library imports.
import pandas as pd
import requests
from io import StringIO
from typing import Union

def request_usgs_data(gage_id, parameter, start_date, start_time, end_date, end_time, gmt_offset, timeout=20):
    """

    :param gage_id [string]: ID of gage to request data for.
    :param start_date [string]: Start date in format yyyy-mm-dd. "2022-06-24"
    :param start_time [string]: Local start time in format HH:MM:SS.mmm. "11:17:05.203"
    :param end_date [string]: End date in format yyyy-mm-dd.
    :param end_time [string]: Local end time in format HH:MM:SS.mmm
    :param gmt_offset [string]: Number of hour offset from gmt (+ or -) in format +/-HH:MM. "-04:00"
    :param parameter [string]: 'discharge' or 'height'
    :param timeout [int or float]: Request timeout in seconds.
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
    print(usgs_url)

    r = requests.get(usgs_url, timeout=timeout)
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

        if skiprows > 1000:
            raise ValueError("USGS text can not be formed into data frame.")

    return df