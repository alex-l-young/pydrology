# ====================================================================================================
# Test script for pydrology/data_request_scripts/usgs_data.py
# ====================================================================================================

# Library imports.
import pandas as pd

# Local imports.
from pydrology.data_requests import request_usgs_data
from pydrology.time_series import resample_data

class TestClassUsgs:
    def setup(self):
        # Test case parameters.
        self.gage_id = "04234000"
        self.parameter = "discharge"
        self.start_date = "2022-05-24"
        self.start_time = "13:00:00.000"
        self.end_date = "2022-05-26"
        self.end_time = "11:00:00.000"
        self.gmt_offset = "-05:00"


    def test_usgs_request_to_dataframe(self):
        gage_df = request_usgs_data(
            self.gage_id,
            self.parameter,
            self.start_date,
            self.start_time,
            self.end_date,
            self.end_time,
            self.gmt_offset
        )

        # List of assertion errors.
        errors = []
        # Test if the gage_df is a Pandas DataFrame.
        if not type(gage_df) == pd.DataFrame:
            errors.append('Gage DataFrame not generated.')
        # Test if data frame contains data.
        if gage_df.empty:
            errors.append("Gage DataFrame is empty.")

        assert not errors, "Errors occured:\n{}".format("\n".join(errors))


    def test_resample_gage_data_downsample(self):
        # Create gage DataFrame.
        gage_df = request_usgs_data(
            self.gage_id,
            self.parameter,
            self.start_date,
            self.start_time,
            self.end_date,
            self.end_time,
            self.gmt_offset
        )

        # Downsample to hourly.
        resample_gage_df = resample_data(gage_df, 'datetime', self.parameter, 60)

        # List of assertion errors.
        errors = []
        # Test if the resample_gage_df is a Pandas DataFrame.
        if not type(resample_gage_df) == pd.DataFrame:
            errors.append('Resampled Gage DataFrame not generated.')
        # Test if data frame contains data.
        if resample_gage_df.empty:
            errors.append("Resampled Gage DataFrame is empty.")
        # Test if data is now hourly.
        sec_diff = (resample_gage_df.loc[1, 'datetime'] - resample_gage_df.loc[0, 'datetime']).total_seconds()
        if sec_diff != 3600:
            errors.append("Resampling was performed incorrectly.")

        assert not errors, "Errors occured:\n{}".format("\n".join(errors))

