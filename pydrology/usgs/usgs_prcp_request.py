# =======================================================
# USGS Stream Gage Data.
# =======================================================

# Library imports.
import pandas as pd
import requests
from io import StringIO

def request_usgs_prcp_data(gage_id, start_date, start_time, end_date, end_time, gmt_offset, timeout=10):
    """

    :param gage_id [string]: ID of gage to request data for.
    :param start_date [string]: Start date in format yyyy-mm-dd. "2022-06-24"
    :param start_time [string]: Local start time in format HH:MM:SS.mmm. "11:17:05.203"
    :param end_date [string]: End date in format yyyy-mm-dd.
    :param end_time [string]: Local end time in format HH:MM:SS.mmm
    :param gmt_offset [string]: Number of hour offset from gmt (+ or -) in format +/-HH:MM. "-04:00"
    :return: Data frame of data from the rain gage.
    """

    # Metadata Request.
    meta_url_format = ("https://waterservices.usgs.gov/nwis/site/?format=rdb&sites={}&"
                       "siteOutput=expanded&siteStatus=all")
    meta_url = meta_url_format.format(gage_id)
    meta_r = requests.get(meta_url, timeout=timeout)

    meta_df = create_metadata_df(meta_r.text)

    # Populate the request URL with ID and dates.
    usgs_url_format = ("https://waterservices.usgs.gov/nwis/iv/?sites={}&parameterCd=00045&"
                       "startDT={}T{}{}&"
                       "endDT={}T{}{}&siteStatus=all&format=rdb")
    usgs_url = usgs_url_format.format(gage_id, start_date, start_time, gmt_offset, end_date, end_time, gmt_offset)
    print(usgs_url)

    r = requests.get(usgs_url, timeout=timeout)
    r_text = r.text

    gage_df = create_discharge_df(r_text)

    # Add metadata.
    gage_df["lat"] = meta_df.loc[0, "dec_lat_va"]
    gage_df["long"] = meta_df.loc[0, "dec_long_va"]
    gage_df["alt"] = meta_df.loc[0, "alt_va"]
    gage_df["alt_acy"] = meta_df.loc[0, "alt_acy_va"]

    return gage_df

def create_metadata_df(usgs_meta_text) -> pd.DataFrame:
    cols = ["agency_cd", "site_no", "station_nm", "site_tp_cd", "lat_va", "long_va", "dec_lat_va", "dec_long_va",
            "coord_meth_cd", "coord_acy_cd", "coord_datum_cd", "dec_coord_datum_cd", "district_cd", "state_cd",
            "county_cd", "country_cd", "land_net_ds", "map_nm", "map_scale_fc", "alt_va", "alt_meth_cd", "alt_acy_va",
            "alt_datum_cd", "huc_cd", "basin_cd", "topo_cd", "instruments_cd", "construction_dt", "inventory_dt",
            "drain_area_va", "contrib_drain_area_va", "tz_cd", "local_time_fg", "reliability_cd", "gw_file_cd",
            "nat_aqfr_cd", "aqfr_cd", "aqfr_type_cd", "well_depth_va", "hole_depth_va", "depth_src_cd", "project_no"]

    skiprows = 10

    # Keep incrementing the number of rows to skip until the correct number is reached.
    while True:
        try:
            df = pd.read_csv(StringIO(usgs_meta_text), header=None, sep='\t', skiprows=skiprows)
            df.columns = cols

            # There are two rows of headers that also need to be skipped.
            if df.loc[0, 'agency_cd'] != 'USGS':
                skiprows += 1
                continue

            break
        except:
            skiprows += 1
            pass

        # Break loop if number of rows is too high.
        if skiprows >= 1000:
            df = None
            break

    return df


def create_discharge_df(usgs_text: str) -> pd.DataFrame:
    cols = ['agency', 'site_no', 'datetime', 'tz_cd', 'precip', 'provis_accept']
    skiprows = 26

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

        # Break loop if number of rows is too high.
        if skiprows >= 1000:
            df = None
            break

    return df

if __name__ == "__main__":
    # Gage ID. Found on the USGS page for the specific monitoring location.
    gage_id = "350110080502045"  # Fall Creek, Ithaca, NY


    # Start date in format yyyy-mm-dd. "2022-06-24"
    start_date = "2022-01-01"

    # Local start time in format HH:MM:SS.mmm. "11:17:05.203"
    start_time = "00:00:00.000"

    # End date in format yyyy-mm-dd. "2022-06-24"
    end_date = "2022-02-01"

    # Local end time in format HH:MM:SS.mmm. "11:17:05.203"
    end_time = "00:00:00.000"

    # Number of hour offset from GMT (+ or -) in format +/-HH:MM. "-04:00"
    gmt_offset = "-05:00"

    df = request_usgs_prcp_data(gage_id, start_date, start_time, end_date, end_time, gmt_offset)

    print(df)