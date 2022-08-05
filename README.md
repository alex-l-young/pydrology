# Pydrology
A Python package for accessing and processing hydrologic data.

## Installation
To install pydrology, run the following within a virtual environment of choice.

```
c:\> pip install pydrology
```

If you want to work on developing Pydrology, you can clone the repository directly from GitHub.
You can also download the example Jupyter Notebooks directly.

```
c:\>git clone https://github.com/alex-l-young/pydrology.git
```

## Currently supported data sources

[**USGS Streamflow Conditions**](https://waterdata.usgs.gov/nwis/rt)
- Stream gage measurements
- Discharge values from rating curves

[**NOAA SC-ACIS Daily Meteorology**](https://scacis.rcc-acis.org/)
- Maximum Temperature
- Minimum Temperature
- Average Temperature
- Precipitation
- Snow Depth

[**NEXRAD Level-II Doppler Radar Reflectivity**](https://www.ncei.noaa.gov/access/metadata/landing-page/bin/iso?id=gov.noaa.ncdc:C00345)
- Currently only base reflectivity.
- Map of stations and coverage available [here](https://www.roc.noaa.gov/WSR88D/Maps.aspx)

## Always Something To-Do

1. Build out a module for basic statistics and plotting that can be done on the hydrology data. 
2. Extract dual-band radar data. And maybe NEXRAD Level-III?
