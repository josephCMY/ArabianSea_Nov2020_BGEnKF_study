#!/bin/bash

#####header for stampede2######
#SBATCH -J btJSIR
#SBATCH -p development
#SBATCH -N 4
#SBATCH -n 272
#SBATCH --time=02:00:00
#SBATCH -o log.track_bt_JSEnKF_IR
#SBATCH -e err.track_bt_JSEnKF_IR


# SCRIPT PURPOSE:
# Run diagnostic enkf program to detect the evolution of obs space biases
# with the number of obs batches assimilated 


# USER SPECIFICATIONS:
# ---------------------
# USER-SPECIFIED CONFIGURATION FILE
export CONFIG_FILE=/work2/04920/tg842199/stampede2/nonlinear_IR-DA/indian_ocean_osse/PSU_EnKF/scripts/config/JointSpace_EnKF_IR_only 

# USER-SPECIFIED DATE RANGE
date_st=201110151200
date_ed=201110151400 #201110170000

# USER-SPECIFIED CYCLING INTERVAL IN MINUTES
cycle_interval=60 

# Location of the diagnostic EnKF
DIAGNOSTIC_PATH=/work/04920/tg842199/stampede2/nonlinear_IR-DA/indian_ocean_osse/PSU_EnKF/code/EnSRF_track_bias/src



# ACTUAL PROGRAM:
# ----------------

# Load config and useful functions
. $CONFIG_FILE
cd /work/04920/tg842199/stampede2/nonlinear_IR-DA/indian_ocean_osse/PSU_EnKF/scripts
. util.sh


# For each time, perform the BT tracking calculation.
date_now=$date_st
while [[ $date_now -le $date_ed ]] ; do

  # Enter run/$date_now/enkf
  cd $WORK_DIR/run/$date_now/enkf

  # Make copy of d01 and enter copied directory
  cp -r d01 diagnose
  cd diagnose

  # Symlink in the diagnostic
  ln -s $DIAGNOSTIC_PATH/enkf.mpi  diagnostic_enkf.mpi

  # Run the diagnostic enkf
  ibrun -n $((NICPU*NJCPU)) diagnostic_enkf.mpi >& log.diagnose 

  echo Obs space diagnostic created for $date_now


  # Advance time
  date_now=`advance_time $date_now $cycle_interval`
done


