#!/bin/bash

# Script to process conventional observations obtained from NCAR CISL's archive

# Useful library
. util.sh

# Config file
. ../config

# Link in directory to output stuff
if [[ ! -e raw_conv ]]; then
  ln -sf $SCRATCH/DOE_reanalysis/obs/conv/raw raw_conv
fi

#------------------------------------------------
# SECTION 1: Remove irrelevant obs 
#------------------------------------------------
date_now=$date_st
while [[ $date_now -lt $date_ed ]]; do
  date_now2=`advance_time $date_now 360` 

  # Remove null characters
  fdate1=`echo $date_now  |cut -c1-10`
  fdate2=`echo $date_now2 |cut -c1-10`
  for prefix in raw_conv/OBS raw_conv/SURFACE_OBS; do
    echo $prefix":"$fdate1
    echo $prefix":"$fdate2
    ./remove_null.sh $prefix":"$fdate1
    ./remove_null.sh $prefix":"$fdate2
  done



  python subset_obs.py $date_now $date_now2
  date_now=`advance_time $date_now 360`
done

#------------------------------------------------
# SECTION 2: Iteratively run obsproc
#------------------------------------------------

# Submit obsproc to knl
salloc -C knl -N 1 -t 48:00:00 -q regular --job-name="obsproc" job_submit.sh
#salloc -C knl -N 1 -t 04:00:00 -q interactive --job-name="obsproc" job_submit.sh
