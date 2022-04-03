#!/bin/bash

# Script to run metgrid.exe to turn ERA5 WPS intermediate files
# to met_em files

# Make namelist
./setup_namelist_wps.sh > namelist.wps

# Evoke metgrid.exe
echo " "
echo "Running metgrid.exe"
#salloc -C knl -N 1 -q interactive -t 04:00:00 --job-name='metgrid' srun -N 1 -n 68 metgrid.exe 
./metgrid.exe >& metgrid.log 

mv metgrid.log.0000 metgrid.log
rm metgrid.log.0*


