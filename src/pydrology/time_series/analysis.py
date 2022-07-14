# ======================================================
# Time series analysis functions.
# ======================================================

# Library imports.
import numpy as np
import pandas as pd
from datetime import timedelta
import matplotlib.pyplot as plt
from scipy import interpolate

# Local imports.
from pydrology.data_request_scripts import usgs_data

# ======================================================
# Functions.

def resample_data(df, dt, resample_col_names, time_col_name):
    """
    Upsample or downsample the hydrology data.
    Upsampling uses linear interpolation to increase the number of samples.
    Downsampling uses a windowed mean to compute the new data value.
    :param df: Pandas Data Frame containing the hydrology data.
    :param dt: New time step in minutes.
    :param resample_col_names: Column names to resample.
        If multiple columns, pass a list of strings.
        If single column, pass a string.
    :param time_col_name: Name of datetime column.
    :return: Pandas DataFrame containing the resampled data.
    """

    # If there is a single resample column, make it a list.
    if type(resample_col_names) == str:
        resample_col_names = [resample_col_names]

    # Convert datetime column to datetime object.
    df[time_col_name] = pd.to_datetime(df[time_col_name])

    # Delta T of the data. Converted to minutes.
    dt_data = (df.loc[1, time_col_name] - df.loc[0, time_col_name]).total_seconds() / 60

    # Resampling factor.
    dt_factor = dt_data / dt

    # Create the new datetime array.
    cur_dt = df.loc[0, time_col_name]
    end_dt = df.loc[df.index[-1], time_col_name]
    datetime_ar = []
    while cur_dt <= end_dt:
        datetime_ar.append(cur_dt)
        cur_dt = cur_dt + timedelta(minutes=dt)

    resampled_columns_dict = {}
    for resample_col in resample_col_names:
        # Use a windowed average to downsample.
        data_resamp = []
        if dt_factor < 1:
            for i, cur_ts in enumerate(datetime_ar):
                # Windowed average.
                data_pts = df.loc[(df.datetime >= cur_ts - timedelta(minutes=dt / 2))
                                  & (df.datetime <= cur_ts + timedelta(minutes=dt / 2)), resample_col]
                mean_val = np.mean(data_pts)
                data_resamp.append(mean_val)

        # Linear interpolation to upsample.
        else:
            # Convert original data to Unix timesteps.
            orig_unix_ts = (df[time_col_name] - pd.Timestamp("1970-01-01")) // pd.Timedelta('1s')

            # Convert interpolation time steps to unix.
            new_unix_ts = (pd.to_datetime(datetime_ar) - pd.Timestamp("1970-01-01")) // pd.Timedelta('1s')

            # Linear interpolation function.
            f = interpolate.interp1d(orig_unix_ts, df[resample_col])

            # Interpolate.
            data_resamp = f(new_unix_ts)

        # Add resampled column data to dictionary.
        resampled_columns_dict[resample_col] = data_resamp


    # Create a new dataframe with the resampled gage data.
    dynamic_col_names = resample_col_names + [time_col_name]
    static_col_names = [col for col in df.columns if col not in dynamic_col_names]
    resample_df_data = {}

    # Add any columns that weren't resampled back into data frame with a constant, repeated value.
    for static_col in static_col_names:
        resample_df_data[static_col] = np.repeat(df.loc[0, static_col], len(datetime_ar))

    # Add resampled data columns.
    for resample_col in resample_col_names:
        resample_df_data[resample_col] = resampled_columns_dict[resample_col]

    # Add time column.
    resample_df_data[time_col_name] = datetime_ar

    # Make data frame.
    resamp_df = pd.DataFrame(resample_df_data)

    return resamp_df


def interpolate_time_series(df, data_column, method='linear'):
    """
    Interpolation of Data Frame columns using Pandas interpolation routine.
    :param df: Data Frame containing data to interpolate.
    :param data_column: Column name to interpolate.
    :param method: Interpolation method. Assumes equally-spaced data.
    :return: Data Frame with column interpolated.
    """

    df[data_column] = df[data_column].interpolate(method=method, limit_direction='both')

    return df


if __name__ == '__main__':
    gage_id = "12451000"
    start_date = "2022-05-24"
    start_time = "13:00:00.000"
    end_date = "2022-05-24"
    end_time = "17:00:00.000"
    gmt_offset = "-05:00"

    gage_df = usgs_data.download_usgs_gage_data(gage_id, start_date, start_time, end_date, end_time, gmt_offset)
    print(gage_df.head())
    print(gage_df.tail())

    resample_gage_df = resample_data(gage_df, 6, ['discharge'], 'datetime')
    print(resample_gage_df.head())
    print(resample_gage_df.tail())

    # Interpolation.
    resample_gage_df.loc[1, 'discharge'] = np.nan
    print(resample_gage_df.head())
    resample_gage_df = interpolate_time_series(resample_gage_df, 'discharge')
    print(resample_gage_df.head())

    fig, ax = plt.subplots()
    ax.plot(gage_df.datetime, gage_df.discharge, 'b-o')
    ax.plot(resample_gage_df.datetime, resample_gage_df.discharge, 'rx')

    plt.show()

