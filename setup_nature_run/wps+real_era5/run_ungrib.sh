#!/bin/bash

# Script to run ungrib.exe for the era5 files

# Make namelist
./setup_namelist_wps.sh > namelist.wps


echo ""

# Link grib files
echo Linking GRIB files
./link_grib.csh raw_era5/*grib

# Run ungrib
echo Running ungrib
./ungrib.exe >& ungrib.log 

# Delete the linked grib files
echo Deleting GRIB file links
rm GRIB*

echo "Finished running ungrib"



## Iterating over both types of GRIB files
#for prefix in era5_reanalysis era5_land; do
#  echo ""
#  
#  # Link grib files
#  echo Linking GRIB files with prefix $prefix
#  ./link_grib.csh raw_era5/"$prefix"*
#
#  # Run ungrib
#  echo Running ungrib
#  #salloc -C knl -N 1 -q interactive -t 01:00:00 --job-name="ungrib" srun -N 1 -n 1 ungrib.exe 
#  ./ungrib.exe >& ungrib.log 
#
#  # Delete the linked grib files
#  echo Deleting GRIB file links
#  rm GRIB*
#
#  # Move rename UNGRIBBED files
#  echo Renaming WPS intermediate files
#  cd processed_era5
#  for ff in `ls --color=never WPS*`; do
#    if [[ $prefix == 'era5_reanalysis' ]] ; then
#      mv $ff PLVL_$ff
#    elif [[ $prefix == 'era5_land' ]] ; then
#      mv $ff LAND_$ff
#    fi
#  done
#  cd ..
#
#  # Proceed to next type of grib file
#
#
#done # End of loop over grib file
#
#
