#!/bin/bash

# Script to run ensemble forecast for a short time. 
# Will be used for validation purposes

# OpenMP settings (seems to run WRF nicely)
export OMP_NUM_THREADS=1
export OMP_PLACES=cores
export OMP_PROC_BIND=spread


# Controls
date_st=201706010000      # Date for the first ensemble forecast
date_ed=201706010000      # Date for the last ensemble forecast
fcst_len=720              # Length of ensemble forecast
fcst_int=30               # Minutes between each fcst output
t_int=180                 # Time interval between each ensemble forecast
expt_list="Thom_GTS+AMV+ch08_30min"

# MPI controls (modify accordingly for interactive nodes)
export SLURM_NTASKS=$((25*68))   # Uncomment and modify this if running on interactive nodes
nnodes_per_wrf=1
nprocs_per_wrf=$((68*$nnodes_per_wrf))
nnodes_total=$((SLURM_NTASKS/68))
tid=0
ntid=$(($nnodes_total/$nnodes_per_wrf))

# Environment settings
source ~/.bashrc

#. util.sh
# Defining advance time function
function advance_time {
  ccyymmdd=`echo $1 |cut -c1-8`
  hh=`echo $1 |cut -c9-10`
  mm=`echo $1 |cut -c11-12`
  inc=$2
  date -u -d $inc' minutes '$ccyymmdd' '$hh':'$mm +%Y%m%d%H%M
}
export -f advance_time

function wrf_time_string {
  ccyy=`echo $1 |cut -c1-4`
  mm=`echo $1 |cut -c5-6`
  dd=`echo $1 |cut -c7-8`
  hh=`echo $1 |cut -c9-10`
  ii=`echo $1 |cut -c11-12`
  echo ${ccyy}-${mm}-${dd}_${hh}:${ii}:00
}
export -f wrf_time_string

function min {
  j=$1
  for i in $@; do
    if [[ $i -lt $j ]]; then
      j=$i
    fi
  done
  echo $j
}
export -f min

# Due to the long run times of WRF, will run all instances of WRF first before running CRTM
for expt in $expt_list; do
  echo "    Running ensemble fcsts on $expt"

  # Identify and read config file
  export CONFIG_FILE=/global/homes/m/my_chan/sumatra_satellite_DA/PSU_EnKF/scripts/config/$expt
  . $CONFIG_FILE

  # For every forecast initiation time
  date=$date_st
  while [[ $date -le $date_ed ]]; do

    # Running wrf on initiation date
    echo "        $date" 

    # Running short ensemble member fcst
    for ee in `seq 1 $NUM_ENS`; do

      eee=`printf %03d $ee`

      cd $WORK_DIR/run/$date/wrf_ens/$eee


      # Check if all needed wrfout files exists
      flag_wrfout_missing=false        # If true, not all wrfout files exist
      date0=$date
      date1=`advance_time $date $fcst_len`
      # Iterate over forecast dates
      while [[ $date0 -le $date1 ]]; do
        wrfout_name=wrfout_d01_`wrf_time_string $date0`
        # If file does not exist, break out of while and for loop
        if [[ ! -e $wrfout_name ]]; then
          flag_wrfout_missing=true
          break 
        fi
        date0=`advance_time $date0 $fcst_int`
      done  # End of while loop.


      # Now, if there are missing wrfout files, run WRF 
      if $flag_wrfout_missing;  then
        
        # Delete namelist.input and namelist.input.backup
        if [[ -e namelist.input ]]; then
          mv namelist.input namelist.input.original
        fi
        if [[ -e namelist.input.backup ]]; then
          mv namelist.input.backup namelist.input.origin
        fi

        # Run namelist_wrf.sh here to construct desired namelist
        export start_date=$date
        export run_minutes=$fcst_len
        export inputout_interval=99999
        export inputout_begin=0
        export inputout_end=99999
        export GET_PHYS_FROM_FILE=false
        if [[ $ee == 1 ]]; then 
          $SCRIPT_DIR/namelist_wrf.sh ensfcst 1 > namelist.input
        else
          ln -s ../001/namelist.input .
        fi

        # Adjust the boundary conditions
        rm wrfbdy_d01
        rm wrflowinp_d01
        ln -s ../../wrf/wrflowinp_d01 .
        ln -s ../../wrf/wrfbdy_d01 .


        # Run wrf with available nodes
        srun -N $nnodes_per_wrf -n $nprocs_per_wrf --cpus-per-task=1 --cpu_bind=cores wrf.exe >& wrf.log &
        tid=$(($tid+1))
        if [[ $tid -ge $ntid ]]; then
          tid=0
          wait
        fi

      fi  # End of if-statement to run WRF

    done # End of ensemble loop

    # Advance date to do the deterministic forecast and all
    date=`advance_time $date $t_int`

  done  # End of fcst init date loop
done  # End of expt loop


wait



