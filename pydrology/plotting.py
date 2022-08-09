# ==================================================
# Plotting functions for data analysis.
# ==================================================

# Library imports.
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ==================================================
# Plotting functions.

def plot_missing_ratio(df, plotting_col, missing_value='M'):

    data = df.loc[:, plotting_col].to_numpy()

    # Total length of data.
    total = len(data)

    # Number of missing values.
    # missing_count = np.sum(data == 'M')
    missing_count = np.sum(np.isin(missing_value, data))

    # Number of non-float values.
    data_float = pd.to_numeric(data, errors='coerce')
    other_nan_count = np.sum(np.isnan(data_float)) - missing_count

    # Number of valid values.
    value_count = total - missing_count - other_nan_count

    # Plot the missing and value ratio.
    labels = [plotting_col]
    width = 1
    fig, ax = plt.subplots()
    # Plot missing value count.
    ax.bar(labels, [missing_count], width, color='red', label='Missing Data')
    ax.text(labels, missing_count / 2, str(missing_count))
    # Plot other nan count.
    ax.bar(labels, [other_nan_count], width, bottom=[missing_count], color='green', label='Other Non-Valid Data')
    ax.text(labels, missing_count + other_nan_count / 2, str(other_nan_count))
    # Plot valid value count.
    ax.bar(labels, [value_count], width, color='skyblue', bottom=[missing_count + other_nan_count], label='Valid Data')
    ax.text(labels, missing_count + other_nan_count + value_count / 2, str(value_count))
    ax.set_ylabel('Data Count')
    ax.legend()

    plt.show()


def plot_data_timeseries(df, data_col, time_col, missing_value='M'):

    # Extract valid data.
    data = df.loc[:, data_col].to_numpy()
    data_valid = pd.to_numeric(data, errors='coerce')
    t = pd.to_datetime(df.loc[:, time_col]) # Time values.

    # Values of t where non-valid data occurs that is not "missing".
    # t_nonvalid = t[(np.isnan(data_valid)) & (data != missing_value)]
    t_nonvalid = t[(np.isnan(data_valid)) & (~np.isin(missing_value, data))]

    fig, ax = plt.subplots(figsize=(8,5))
    ax.plot(t, data_valid, 'b-+', label='Valid Data')

    # Plot missing data as red vertical lines.
    if np.isin(missing_value, data):
        t_missing = t[np.squeeze(np.argwhere(data == missing_value))]

        for i, m in enumerate(t_missing):
            if i == 0:
                ax.axvline(m, c='red', label='Missing Data')
            else:
                ax.axvline(m, c='red')

    # Plot other non-valid data as vertical green lines.
    for i, nv in enumerate(t_nonvalid):
        if i == 0:
            ax.axvline(nv, c='green', label='Non-Valid Data')
        else:
            ax.axvline(nv, c='green')

    ax.legend()
    ax.set_xlabel('Time')
    ax.set_ylabel('Value')


    # # Plot the missing values as filled vertical areas.
    # miss_flag = False
    # start_miss_t = 0
    # end_miss_t = 0
    # for i, data_val in enumerate(data):
    #     # Start of missing values.
    #     if miss_flag is False and np.isnan(data_val):
    #         start_miss_t = t[i]
    #         miss_flag = True
    #
    #     # End of missing values.
    #     if miss_flag is True and ~np.isnan(data_val):
    #         end_miss_t = t[i]
    #         ax.axvspan(start_miss_t, end_miss_t, alpha=0.5, color='red')
    #         miss_flag = False

    plt.show()


# def plot_nexrad_sweep(nexrad_nc):


if __name__ == '__main__':
    # data = {
    #     'a': [1, 8, 7, 4, 'M', 'M', 5, 4, 'T'],
    #     'b': [
    #         '2022-01-01',
    #         '2022-01-02',
    #         '2022-01-03',
    #         '2022-01-04',
    #         '2022-01-05',
    #         '2022-01-06',
    #         '2022-01-07',
    #         '2022-01-08',
    #         '2022-01-09',
    #     ]
    # }
    # df = pd.DataFrame(data)
    #
    # # plot_missing_ratio(df, 'a')
    #
    # plot_data_timeseries(df, 'a', 'b', missing_value='M')

    from pydrology.data_request_scripts import usgs_data

    # Parameters for the request.
    # ---------------------------

    # Gage ID. Found on the USGS page for the specific monitoring location.
    gage_id = "04234000"  # Fall Creek, Ithaca, NY

    # Start date in format yyyy-mm-dd. "2022-06-24"
    start_date = "2021-01-24"

    # Local start time in format HH:MM:SS.mmm. "11:17:05.203"
    start_time = "13:00:00.000"

    # End date in format yyyy-mm-dd. "2022-06-24"
    end_date = "2022-05-26"

    # Local end time in format HH:MM:SS.mmm. "11:17:05.203"
    end_time = "11:00:00.000"

    # Number of hour offset from GMT (+ or -) in format +/-HH:MM. "-04:00"
    gmt_offset = "-05:00"

    # Request the gage data as a DataFrame.
    gage_df = usgs_data.download_usgs_gage_data(gage_id, start_date, start_time, end_date, end_time, gmt_offset)

    # Print the head and tail.
    gage_df.head()
    gage_df.tail()

    # Plotting column names and missing value.
    data_column_name = 'discharge'
    time_column_name = 'datetime'
    missing_value = 'M'

    # Plot the valid, missing, and non-valid data as a bar chart.
    plotting.plot_missing_ratio(gage_df, data_column_name)

    # Plot the data as a time series.
    plotting.plot_data_timeseries(gage_df, data_column_name, time_column_name, missing_value=missing_value)

