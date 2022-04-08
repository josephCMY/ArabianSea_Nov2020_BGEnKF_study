# Script to evoke CRTM on nature run
# Meant to be called within KNL nodes

# Library of useful functions
. util.sh

# Date range
date_st=202011061200
date_ed=202011110300
date_interval=30          # 30 minute interval


# CRTM parallel run controls
crtm_procs=68
tot_procs=272


# Enter the nature run
cd nature_run


# Start iterating through dates
date_nw=$date_st
active_procs=0
while [[ $date_nw -le $date_ed ]]; do

  # Define file names
  wrffile='wrfout_d01_'`wrf_time_string $date_nw`
  btfile=met8_bt_"$date_nw".nc

  # Evoke the CRTM and count number of active processes
  srun -n $crtm_procs --ntasks-per-node=68 --cpus-per-task=1 --cpu_bind=cores crtm.exe $wrffile $btfile >& log.crtm_$date_nw &
  sleep 3
  echo Running CRTM on $date_nw
  active_procs=$(( $active_procs + $crtm_procs ))

  # Wait if all processes are running
  if [[ $active_procs -ge $tot_procs ]]; then
    wait
    active_procs=0
  fi
    
  # Move onto next date
  date_nw=`advance_time $date_nw $date_interval`

done
