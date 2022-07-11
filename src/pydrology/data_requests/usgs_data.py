# ======================================================
# Utility functions for processing stream data.
# ======================================================

# Library imports.
import requests
import numpy as np
import pandas as pd
from io import StringIO
from datetime import timedelta

# Local imports.


# =====================================================
# Functions.
def download_usgs_gage_data(gage_id, start_date, start_time, end_date, end_time, gmt_offset):
    """

    :param gage_id [string]: ID of gage to request data for.
    :param start_date [string]: Start date in format yyyy-mm-dd. "2022-06-24"
    :param start_time [string]: Local start time in format HH:MM:SS.mmm. "11:17:05.203"
    :param end_date [string]: End date in format yyyy-mm-dd.
    :param end_time [string]: Local end time in format HH:MM:SS.mmm
    :param gmt_offset [string]: Number of hour offset from gmt (+ or -) in format +/-HH:MM. "-04:00"
    :return:
    """

    # Populate the request URL with ID and dates.
    usgs_url_format = ("https://waterservices.usgs.gov/nwis/iv/?sites={}&parameterCd=00065&"
                       "startDT={}T{}{}&endDT={}T{}{}&"
                       "siteStatus=all&format=rdb")
    usgs_url = usgs_url_format.format(gage_id, start_date, start_time, gmt_offset, end_date, end_time, gmt_offset)

    r = requests.get(usgs_url)
    r_text = r.text

    gage_df = create_discharge_df(r_text)

    return gage_df


def create_discharge_df(usgs_text):
    cols = ['agency', 'site_no', 'datetime', 'tz_cd', 'discharge', 'provis_accept']
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


def resample_gage_data(gage_df, dt, gage_data_col_name):
    # TODO: Make this a general resampling algorithm.

    # Convert datetime column to datetime object.
    gage_df['datetime'] = pd.to_datetime(gage_df['datetime'])

    # Delta T of the data. Converted to minutes.
    dt_data = (gage_df.loc[1, 'datetime'] - gage_df.loc[0, 'datetime']).total_seconds() / 60

    # Resampling factor.
    dt_factor = dt_data / dt

    # Create the new datetime array.
    cur_dt = gage_df.loc[0, 'datetime']
    end_dt = gage_df.loc[gage_df.index[-1], 'datetime']
    datetime_ar = []
    while cur_dt <= end_dt:
        datetime_ar.append(cur_dt)
        cur_dt = cur_dt + timedelta(minutes=dt)

    # Use trapezoidal integration to downsample.
    gage_resamp = []
    if dt_factor < 1:
        for i in range(1, len(datetime_ar)):
            prev_ts = datetime_ar[i-1]
            cur_ts = datetime_ar[i]
            # Collect time steps to integrate over.
            # TODO: Could be made faster with sequential indexing isntead of searching for >= and <.
            int_data = gage_df.loc[(gage_df.datetime >= prev_ts) & (gage_df.datetime < cur_ts), gage_data_col_name]

            # Trapezoidal integration. Use seconds for the units of dx since discharge is in cfs.
            cumu_gage = np.trapz(int_data, dx=dt_data * 60)
            gage_resamp.append(cumu_gage)

    # Linear interpolation to upsample.
    else:
        # TODO: Write upsampling code.
        x = 1

    # Create a new dataframe with the resampled gage data.
    resamp_df_data = {
        "agency": np.repeat(gage_df.loc[0, "agency"], len(gage_resamp)),
        "site_no": np.repeat(gage_df.loc[0, "site_no"], len(gage_resamp)),
        "datetime": datetime_ar[:len(gage_resamp)],
        "tz_cd": np.repeat(gage_df.loc[0, "tz_cd"], len(gage_resamp)),
        gage_data_col_name: gage_resamp,
        "provis_accept": np.repeat(gage_df.loc[0, "provis_accept"], len(gage_resamp)),
    }
    gage_resamp_df = pd.DataFrame(resamp_df_data)

    return gage_resamp_df


if __name__ == '__main__':
    gage_id = "12451000"
    start_date = "2022-05-24"
    start_time = "13:00:00.000"
    end_date = "2022-05-26"
    end_time = "11:00:00.000"
    gmt_offset = "-05:00"

    gage_df = download_usgs_gage_data(gage_id, start_date, start_time, end_date, end_time, gmt_offset)
    print(gage_df.head())
    print(gage_df.tail())

    resample_gage_df = resample_gage_data(gage_df, 60, 'discharge')
    print(resample_gage_df.head())
    print(resample_gage_df.tail())