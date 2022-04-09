#!/bin/bash

#####header for stampede2######
#SBATCH -J sepIR
#SBATCH -p normal
#SBATCH -N 7
#SBATCH -n 228
#SBATCH --time=24:00:00
#SBATCH -o log.diagnose_IR


# Script to call EnSRF_diagnostic, which will output the ensemble mean thrice:
# 1) True prior ensemble mean
# 2) Ensemble mean after assimilating non-IR obs
# 3) Ensemble mean after assimilating IR obs

# Dates to do diagnostics for
date_st=201110160000
date_ed=201110171200

# DA cycling time interval in minutes
t_int=60

NUM_ENS=49

nprocs=$((12*19))

# Directory containin all expts
mother_dir='/scratch/04920/tg842199/nonlinear_IR_DA/indian_ocean_osse/PSU_EnKF'

# Expts to process
exptlist='BGEnKF_NoQ EnKF_NoQ'

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

    # Make diagnostic directory
    rm -r enkf_diagnostic
    mkdir enkf_diagnostic
    cd enkf_diagnostic
    ln -s ../enkf/d01/* .
    rm enkf.log fort.* enkf_time.log enkf.mpi
    ln -s ../enkf/d01/fort.8* .
    ln -s $mother_dir/$expt/code/EnSRF_diagnose/src/enkf.mpi .
    rm fort.$((80010+$NUM_ENS+1))
    cp fort.80011 fort.$((80010+$NUM_ENS+1))
    cp fort.80011 fort.$((30010+$NUM_ENS+1))
    cp fort.80011 fort.$((40010+$NUM_ENS+1))
    ibrun -n 228 enkf.mpi >& enkf.log 

  done # --------------------- End of loop over expts

  # Advance time
  date_nw=`advance_time $date_nw $t_int`

done # ----------------------- End of loop over dates
