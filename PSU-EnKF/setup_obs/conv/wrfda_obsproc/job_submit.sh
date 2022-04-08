#!/bin/bash
#####header for Cori######
#SBATCH -J DOE 
#SBATCH -q interactive
#SBATCH -C knl 
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time=04:00:00
#SBATCH -o log.obsproc

# Script to iteratively run obsproc.exe 

# Useful library
. util.sh

# Config file
. ~/DOE_reanalysis/obs/conv/config

# Iterative loop
date_now=$date_st
tint=$cyc_int
date_end=$date_ed

while [[ $date_now -le $date_end ]]; do

  echo Running obsproc for $date_now
   
  # link relevant file
  if [[ -e obs.raw ]]; then
    rm obs.raw
  fi
  ln -s raw_conv/subsetted_obs:$date_now obs.raw

  # Generate namelist
  ./namelist_obsproc.sh $date_now > namelist.obsproc

  # Run the obsproc.exe
  srun -N 1 -n 1 obsproc.exe >& obsproc.log 

  # Move the outputs to a SCRATCH directory
  mv obs_* raw_conv

  # Increment time
  date_now=`advance_time $date_now $tint`

done

