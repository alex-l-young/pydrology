# ====================================================================================================
# Test script for pydrology/data_request_scripts/usgs_data.py
# ====================================================================================================

# Library imports.
import pandas as pd

# Local imports.
from pydrology.acis.acis_request import request_acis_data

class TestClassAcis:
    def setup(self):
        self.site_id = '304174 2'  # Ithaca.
        self.met_elements = ['maxt', 'mint', 'snow', 'avgt', 'pcpn']
        self.start_date = '2022-01-01'  # YYYY-mm-dd.
        self.end_date = '2022-01-21'

    def test_request_acis_data(self):
        acis_df = request_acis_data(
            self.met_elements,
            self.site_id,
            self.start_date,
            self.end_date
        )

        # List of assertion errors.
        errors = []
        # Test if the gage_df is a Pandas DataFrame.
        if not type(acis_df) == pd.DataFrame:
            errors.append('Gage DataFrame not generated.')
        # Test if data frame contains data.
        if acis_df.empty:
            errors.append("Gage DataFrame is empty.")
        # Test if the correct columns are present.
        cols = ['Date', 'MaxTemp', 'MinTemp', 'SnowDepth', 'AvgTemp', 'Precipitation', 'site_id', 'name', 'latitude', 'longitude']
        if set(cols) != set(acis_df.columns):
            errors.append("Not all ACIS columns are present.")

        assert not errors, "Errors occured:\n{}".format("\n".join(errors))