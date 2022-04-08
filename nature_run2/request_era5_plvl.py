''' Script to request for ERA-5 reanalysis data within a limited area for a limited time '''

import cdsapi
import numpy as np
import datetime
import sys

date = datetime.datetime.strptime( sys.argv[1], '%Y%m%d%H%M')


c = cdsapi.Client()

c.retrieve(
    'reanalysis-era5-pressure-levels',
    {
        'product_type':'reanalysis',
        'format':'grib',
        'variable':[  'geopotential','fraction_of_cloud_cover','ozone_mass_mixing_ratio','relative_humidity',
                      'specific_cloud_ice_water_content','specific_cloud_liquid_water_content',\
                       'specific_humidity','specific_rain_water_content','specific_snow_water_content',
                       'temperature','u_component_of_wind','v_component_of_wind','vertical_velocity'
        ],
        'pressure_level':[ '1','2','3', '5','7','10','20','30','50', '70','100','125','150','175','200',
                           '225','250','300','350','400','450','500','550','600','650','700','750',
                           '775','800','825','850','875','900','925','950','975','1000' ],
        'year': date.strftime('%Y'),
        'month': date.strftime( '%m' ),
        'day': date.strftime( '%d' ),
        'area': [22.5, 47.5, -22.5, 102.5], 
        'time': date.strftime( '%H' )
        },
   'raw_era5/era5_reanalysis_' + date.strftime('%Y-%m-%d_%H')+'UTC.grib')


