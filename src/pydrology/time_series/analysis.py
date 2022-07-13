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
from pydrology.data_requests import usgs_data

# ======================================================
# Functions.

def resample_data(df, dt, data_col_name, time_col_name):

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

    # Use a windowed average to downsample.
    data_resamp = []
    if dt_factor < 1:
        for i, cur_ts in enumerate(datetime_ar):
            # Windowed average.
            data_pts = df.loc[(df.datetime >= cur_ts - timedelta(minutes=dt / 2))
                              & (df.datetime <= cur_ts + timedelta(minutes=dt / 2)), data_col_name]
            mean_val = np.mean(data_pts)
            data_resamp.append(mean_val)

    # Linear interpolation to upsample.
    else:
        # Convert original data to Unix timesteps.
        orig_unix_ts = (df[time_col_name] - pd.Timestamp("1970-01-01")) // pd.Timedelta('1s')

        # Convert interpolation time steps to unix.
        new_unix_ts = (pd.to_datetime(datetime_ar) - pd.Timestamp("1970-01-01")) // pd.Timedelta('1s')

        # Linear interpolation function.
        f = interpolate.interp1d(orig_unix_ts, df[data_col_name])

        # Interpolate.
        data_resamp = f(new_unix_ts)


    # Create a new dataframe with the resampled gage data.
    resamp_df_data = {
        "agency": np.repeat(df.loc[0, "agency"], len(data_resamp)),
        "site_no": np.repeat(df.loc[0, "site_no"], len(data_resamp)),
        "datetime": datetime_ar[:len(data_resamp)],
        "tz_cd": np.repeat(df.loc[0, "tz_cd"], len(data_resamp)),
        data_col_name: data_resamp,
        "provis_accept": np.repeat(df.loc[0, "provis_accept"], len(data_resamp)),
    }
    resamp_df = pd.DataFrame(resamp_df_data)

    return resamp_df

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

    resample_gage_df = resample_data(gage_df, 6, 'discharge', 'datetime')
    print(resample_gage_df.head())
    print(resample_gage_df.tail())

    fig, ax = plt.subplots()
    ax.plot(gage_df.datetime, gage_df.discharge, 'b-o')
    ax.plot(resample_gage_df.datetime, resample_gage_df.discharge, 'rx')

    plt.show()

