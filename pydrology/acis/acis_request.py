# =======================================================
# NOAA SC-ACIS Meteorology.
# =======================================================

# Library imports.
import pandas as pd
import numpy as np
import requests
import os

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
                        'snow': 'SnowFall', 'snwd': 'SnowDepth'}

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

