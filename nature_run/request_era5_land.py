''' Script to request for ERA-5 reanalysis data within a limited area for a limited time '''

import cdsapi
import numpy as np
import datetime
import sys

date = datetime.datetime.strptime( sys.argv[1], '%Y%m%d%H%M')


c = cdsapi.Client()

c.retrieve(
    'reanalysis-era5-single-levels',
    {
        'product_type':'reanalysis',
        'format':'grib',
        'variable': ['10m_u_component_of_wind','10m_v_component_of_wind','2m_dewpoint_temperature',
                    '2m_temperature','land_sea_mask','mean_sea_level_pressure',
                    'orography','sea_surface_temperature','skin_temperature',
                    'soil_temperature_level_1','soil_temperature_level_2','soil_temperature_level_3',
                    'soil_temperature_level_4','soil_type','standard_deviation_of_filtered_subgrid_orography',
                    'surface_pressure','volumetric_soil_water_layer_1','volumetric_soil_water_layer_2',
                    'volumetric_soil_water_layer_3','volumetric_soil_water_layer_4' ],
        'area': [20.0, 50.0, -20.0, 100.],
        'year':date.strftime('%Y'),
        'month':date.strftime('%m'),
        'day': date.strftime('%d'),
        'time': date.strftime('%H')
    },
    'raw_era5/era5_land_' + date.strftime('%Y-%m-%d_%H')+'UTC.grib')
