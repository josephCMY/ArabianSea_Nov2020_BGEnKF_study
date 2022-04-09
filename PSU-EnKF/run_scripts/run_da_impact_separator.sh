#!/bin/bash

#####header for stampede2######
#SBATCH -J retry
#SBATCH -p normal
#SBATCH -N 7
#SBATCH -n 228
#SBATCH --time=03:00:00
#SBATCH -o log.diagnose_IR


# Script to run EnKF on the BGEnKF priors
# Done to attribute whether BGEnKF's QVAPOR performance is due to weaker CR.

# Dates to do diagnostics for
date_st=201110171800
date_ed=201110180900

# DA cycling time interval in minutes
t_int=180

NUM_ENS=49

nprocs=$((12*19))

# Directory containin all expts
mother_dir='/scratch/04920/tg842199/nonlinear_IR_DA/indian_ocean_osse/PSU_EnKF'

# Expts to process
exptlist='BGEnKF_3hr'

# Load convenient functions
. util.sh


# Iterate thru all dates
date_nw=$date_st
while [[ $date_nw -le $date_ed ]]; do


  # Iterate thru all expts
  for expt in $exptlist; do 

    # Enter relevant directory
    rundir=$mother_dir"/"$expt/run/$date_nw
    cd $rundir

    # Make directory to try enkf instead of bgenkf
    rm -r enkf_retry
    mkdir enkf_retry
    cd enkf_retry

    # Setting up directory
    ln -s ../enkf/d01/* .
    rm enkf.log fort.* enkf_time.log namelist.enkf
    ln -s ../enkf/d01/fort.8* .
    cp ../enkf/d01/fort.9* .
    rm fort.$((80010+$NUM_ENS+1)) fort.$((90010+$NUM_ENS+1))
    cp fort.80011 fort.$((80010+$NUM_ENS+1))
    cp fort.80011 fort.$((90010+$NUM_ENS+1))

    # Copying over namelist from EnKF expt
    cp $mother_dir"/EnKF"/run/$date_nw/enkf/d01/namelist.enkf .

    # Run EnKF
    ibrun -n 228 enkf.mpi >& enkf.log 

  done # --------------------- End of loop over expts

  # Advance time
  date_nw=`advance_time $date_nw $t_int`

done # ----------------------- End of loop over dates
