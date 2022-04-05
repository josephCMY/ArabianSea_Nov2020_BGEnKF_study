#!/bin/bash

# Bash script to call down data from TIGGE and ERA5

# Load configuration
. config

# Function to increment time
function advance_time {
  ccyymmdd=`echo $1 |cut -c1-8`
  hh=`echo $1 |cut -c9-10`
  mm=`echo $1 |cut -c11-12`
  inc=$2
  date -u -d $inc' minutes '$ccyymmdd' '$hh':'$mm +%Y%m%d%H%M
}

# Function to loop over times requested
function request_date_loop { 

  # Function inputs
  CONFIG_FILE=$1
  python_script=$2

  # Starting loop
  date=$date_st

  while [[ $date -le $date_ed ]]; do
    python3 $python_script $date
    date=`advance_time $date $t_int`  
  done
}


# MAIN BASH SCRIPT
request_date_loop config request_era5_plvl.py >& log.era5_plvl &
request_date_loop config request_era5_land.py >& log.era5_land &
#python3 request_tigge_plvl.py 201110150000
#python3 request_tigge_land.py 201110150000
#request_date_loop config request_tigge_plvl.py >& log.tigge_plvl &
#request_date_loop config request_tigge_land.py >& log.tigge_land &

# Wait for requests to terminate before quitting.
wait

