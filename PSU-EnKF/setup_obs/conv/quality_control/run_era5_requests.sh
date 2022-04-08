#!/bin/bash
# Script to retrieve ERA5 data for quality control

# Configuration file
. config

# Time interval between ERA5 files
era5_tint=60

# Retrieve hourly ERA5 files
DATE=$date_st
while [[ $DATE -le $date_ed ]]; do

  python request_era5.py $DATE
  # Increment time
  echo Retrieved ERA5 for $DATE
  DATE=`advance_time $DATE $era5_tint`
done


## Generate ERA5 for off-hour times
#DATE0=$date_st
#DATE1=`advance_time $date_st $era5_tint`
#while [[ $DATE1 -le $date_ed ]]; do
#
#  DATE=`advance_time $DATE0 $((era5_tint/2))`
#  f1='era5_data/era5_reanalysis_'$DATE0'.nc'
#  f2='era5_data/era5_reanalysis_'$DATE1'.nc'
#  fmid='era5_data/era5_reanalysis_'$DATE'.nc'
#  
#  # Take average
#  ncea $f1 $f2 $fmid
#
#  # Increment time
#  echo Generated ERA5 file for $DATE
#  DATE0=`advance_time $DATE0 $era5_tint`
#  DATE1=`advance_time $DATE1 $era5_tint`
#
#done

