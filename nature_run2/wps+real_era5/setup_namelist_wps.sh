#!/bin/bash

# Script to make namelist.wps for ungrib.exe and metgrid.exe

# Load config file
. ../config

# Construct string
cat << EOF
&share
 wrf_core                     = 'ARW',
 max_dom                      = 1,
 start_year                   = `echo $date_st |cut -c1-4`,
 start_month                  = `echo $date_st |cut -c5-6`,  
 start_day                    = `echo $date_st |cut -c7-8`,  
 start_hour                   = `echo $date_st |cut -c9-10`,  
 start_minute                 = `echo $date_st |cut -c11-12`,  
 start_second                 = 00,  
 end_year                     = `echo $date_ed |cut -c1-4`,
 end_month                    = `echo $date_ed |cut -c5-6`,  
 end_day                      = `echo $date_ed |cut -c7-8`,  
 end_hour                     = `echo $date_ed |cut -c9-10`,  
 end_minute                   = `echo $date_ed |cut -c11-12`,  
 end_second                   = 00,  
 interval_seconds             = $(($t_int*60)), 
 io_form_geogrid              = 2, 
 opt_output_from_geogrid_path = './',
 debug_level                  = 0,
/

&geogrid
 parent_id                    = 0,
 parent_grid_ratio            = 1, 
 i_parent_start               = 1,
 j_parent_start               = 1,
 geog_data_res                = 30s,
 e_we                         = 431,
 e_sn                         = 401,
 dx                           = 9000,
 dy                           = 9000,
 map_proj                     = 'mercator',
 ref_lat                      =  0.0,
 ref_lon                      = 73.0,
 stand_lon                    = 73.0,
 truelat1                     = 0.0,
 truelat2                     = 0.0, 
 geog_data_res                = '30s',
 geog_data_path               = '/global/cscratch1/sd/my_chan/geog',
 opt_geogrid_tbl_path         = './geogrid'
/

&ungrib
 out_format                   = 'WPS'
 prefix                       = './processed_era5/WPS',
/

&metgrid
 fg_name = './processed_era5/PLVL','./processed_era5/LAND',
 io_form_metgrid = 2, 
 opt_output_from_metgrid_path = './processed_era5',
 opt_metgrid_tbl_path         = './metgrid',
/
EOF
