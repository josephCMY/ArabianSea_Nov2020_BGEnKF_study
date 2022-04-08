''' Script to request for ERA-5 reanalysis data within a limited area for a limited time '''

import cdsapi
import numpy as np
import datetime
import sys

date = datetime.datetime.strptime( sys.argv[1], '%Y%m%d%H%M')

outname = 'era5_data/era5_reanalysis_' + date.strftime('%Y%m%d%H%M')+'.nc'

c = cdsapi.Client()

c.retrieve(
    'reanalysis-era5-pressure-levels',
    {
        'product_type':'reanalysis',
        'format':'netcdf',
        'variable':[ 'geopotential', 'specific_humidity',
                     'temperature','u_component_of_wind',
                     'v_component_of_wind' ],
        'pressure_level':[ '1','2','3', '5','7','10','20','30','50', '70','100','125','150','175','200',
                           '225','250','300','350','400','450','500','550','600','650','700','750',
                           '775','800','825','850','875','900','925','950','975','1000' ],
        'year': '2017',
        'month': date.strftime( '%m' ),
        'day': date.strftime( '%d' ),
        'area': [42.5, 57.5,-22.5, 142.5],
        'time': date.strftime( '%H' )
        },
   outname)
