#!/bin/bash

# Script to run deterministic forecasts

#####header for Cori######
#SBATCH -J deterfcst
#SBATCH -q regular
#SBATCH -C knl
#SBATCH -N 13 
#SBATCH -n 1700
#SBATCH --time=12:00:00
#SBATCH -o log.Thom_GTS
#SBATCH -e err.Thom_GTS

# OpenMP settings (seems to run WRF nicely)
export OMP_NUM_THREADS=1
export OMP_PLACES=cores
export OMP_PROC_BIND=spread


# Controls
date_st=201705310300      # Date for the first deterministic forecast
date_ed=201706020000      # Date for the last deterministic forecast
fcst_len=720              # Length of deterministic forecast
fcst_int=60               # Minutes between each fcst output
t_int=360                 # Time interval between each deterministic forecast
flag_wrf=false            # True = run deterministic fcst. False = no WRF integration.
flag_crtm=true            # True = CRTM on the deterministic fcst. False = no CRTM
expt_list="Thom_GTS+AMV+ch08_3hrly" 
#for expt in 030min_infl 060min 180min_infl; do
#  expt_list=$expt_list" "Thom_GTS+Him8_LM_$expt
#done


# MPI controls (modify accordingly for interactive nodes)
export SLURM_NTASKS=$((26*68))   # Uncomment and modify this if running on interactive nodes
nnodes_per_wrf=1
nprocs_per_wrf=$((68*$nnodes_per_wrf))
nnodes_per_crtm=1
nprocs_per_crtm=$((68*$nnodes_per_crtm))
nnodes_used=0
nnodes_total=$((SLURM_NTASKS/68))

# Environment settings
source ~/.bashrc

# Useful utilities
function advance_time {
  ccyymmdd=`echo $1 |cut -c1-8`
  hh=`echo $1 |cut -c9-10`
  mm=`echo $1 |cut -c11-12`
  inc=$2
  date -u -d $inc' minutes '$ccyymmdd' '$hh':'$mm +%Y%m%d%H%M
}
function wrf_time_string {
  ccyy=`echo $1 |cut -c1-4`
  mm=`echo $1 |cut -c5-6`
  dd=`echo $1 |cut -c7-8`
  hh=`echo $1 |cut -c9-10`
  ii=`echo $1 |cut -c11-12`
  echo ${ccyy}-${mm}-${dd}_${hh}:${ii}:00
}

# Moving onto running crtm
if $flag_crtm; then

  echo "Running CRTM on deter. fcst."

  # For each experiment
  for expt in $expt_list; do
    echo "    Running CRTM on $expt deter. fcsts."
  
    # Identify and read config file
    CONFIG_FILE=/global/homes/m/my_chan/sumatra_satellite_DA/PSU_EnKF/scripts/config/$expt
    . $CONFIG_FILE
  
    # For every forecast initiation time
    date=$date_st
    while [[ $date -le $date_ed ]]; do

      echo "        "$date
  
      cd $WORK_DIR/run/$date/wrf

      # Symbolic link crtm.exe in
      ln -sf $CODE_DIR/CRTM/crtm_wrf/him8_sumatra/crtm.exe crtm.exe 

      # For every forecasted time
      date0=$date
      date1=`advance_time $date $fcst_len`
      while [[ $date0 -le $date1 ]]; do
        echo Dealing with $date0

        # dont run crtm if required wrfout file does not exist
        wrfout_name=wrfout_d01_`wrf_time_string $date0`
        if [[ ! -e $wrfout_name ]]; then
          echo "    " $wrfout_name does not exist. skipping over
          date0=`advance_time $date0 $fcst_int`
          continue
        fi

        # Check if nodes are available
        if [[ $(( $nnodes_used + $nnodes_per_crtm )) -gt $nnodes_total ]]; then
          echo "        "Insufficient nodes available.
          echo "        "Waiting for background processes to clear
          wait
          echo "        "Background processes cleared!
          nnodes_used=0
        fi

        # Run CRTM
        srun -N $nnodes_per_crtm -n $nprocs_per_crtm crtm.exe $wrfout_name him8_df_$date0.bin >& crtm_"$date0".log &
        nnodes_used=$(( $nnodes_used + $nnodes_per_crtm ))

        # Go to fc directory and set up symbolic links to relevant files
        cd $WORK_DIR/fc/$date
        ln -sf $WORK_DIR/run/$date/wrf/him8_df_$date0.bin .
        ln -sf $WORK_DIR/run/$date/wrf/$wrfout_name df_$wrfout_name
        cd $WORK_DIR/run/$date/wrf

        # Advance time
        date0=`advance_time $date0 $fcst_int`
        
      done # End of while loop over fcsted times

      # Advance to next forecast init date
      date=`advance_time $date $t_int`

    done  # End of while loop over fcst init dates
  done  # End of loop over expt
fi  # End of if-statment for crtm

wait


