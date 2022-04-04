#!/bin/bash

# Script to run ungrib.exe for the era5 files

## Make namelist
#./setup_namelist_wps.sh > namelist.wps


# Symlink in all ERA5 GRIB files
./link_grib.csh raw_era5/*grib


# Run ungrib
echo Running ungrib
salloc -C knl -N 1 -n 1 -t 04:00:00 -q interactive srun ungrib.exe >& ungrib_srun.log 


# Clean up directory
rm GRIB*




