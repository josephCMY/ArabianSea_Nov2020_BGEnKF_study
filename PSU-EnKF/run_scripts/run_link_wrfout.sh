#!/bin/bash

# Script to run nccopy compressions on wrfinput files

# Useful settings
source ~/.bashrc

# Useful functions
. util.sh

# Date range
date_st=201705312100
date_ed=201706010000
t_int=30

# Number of wrf domains
ndom=1

# For each experiment, run CRTM
for expt in Thom_GTS+ch08_30min; do

  echo Running nccopy for $expt
  #load configuration files, functions, parameters
  export CONFIG_FILE=/global/homes/m/my_chan/sumatra_satellite_DA/PSU_EnKF/scripts/config/$expt

  . $CONFIG_FILE
#  NUM_ENS=3


  # Initiate date
  date=$date_st

  # Loop over dates in question
  while [[ $date -le $date_ed ]]; do


    echo "    "Processing $date

    # Dealing with wrfout files in wrf_ens directory
    echo "    "Linking files in wrf_ens ens directories

    # Identify wrfout file names
    cd $WORK_DIR/run/$date/wrf_ens/001
    wrfout_list=`ls wrfout*`

    # Enter fc dir
    cd $WORK_DIR/fc/$date

    # Iterate over members
    for ee in `seq -f "%03g" 1 $NUM_ENS`; do

      # Iterate over available wrfout files
      for f_in in $wrfout_list; do

        ln -sf $WORK_DIR/run/$date/wrf_ens/$ee/$f_in $f_in"_"$ee

      done

    done
   

    # Increment time
    date=`advance_time $date $t_int`

  done # End of loop over time

done # End of loop over expt

# Wait for all background processes to clear before exiting
wait


