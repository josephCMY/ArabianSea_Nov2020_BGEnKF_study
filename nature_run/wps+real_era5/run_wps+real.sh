#!/bin/bash

# Script to run the WPS and the real.exe for ERA5

cd processed_era5

lfs setstripe --stripe-count 8 .

cd ~/DOE_reanalysis/setup_ens/wps+real_era5

# Run ungrib
./run_ungrib.sh

# Run metgrid
./run_metgrid.sh

# Run real.exe
./run_real.sh
