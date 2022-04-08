#!/bin/bash

# Script to process conventional observations obtained from ncar cisl's archive


# config file
. config

# Request for ERA5 reanalysis 
#./run_era5_requests.sh

# Perform QC
python qc_gts_obs.py  $date_st $date_ed
