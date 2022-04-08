#!/bin/bash

# Special script to submit an ensemble of WRF runs

# Load config 
. $CONFIG_FILE


# Iterate over all ensemble members
tid=0
ntid=25

echo Starting WRF ensemble on `date`
for NE in `seq -f "%03g" 1 $NUM_ENS`; do

  id=$NE  #`expr $NE + 1000 |cut -c2-`
  if [[ ! -d $id ]]; then mkdir $id; fi
  touch $id/rsl.error.0000
  if [[ `tail -n 20 $id/rsl.error.0000 |grep SUCCESS` ]]; then 
    echo " id $id, `tail -n 20 $id/rsl.error.0000 |grep SUCCESS`"
    continue
  fi
  echo "  id $id, RUNNING WRF"

  cd $id
  #lfs setstripe -c 1 .


#  # Check if wrflowinp has been constructed
#  # The wrflowinp construction step was thrown into the background earlier. 
#  while [[ ! -f wrflowinp_d01 ]]; do
#    sleep 5
#  done

  # Run WRF
#  $SCRIPT_DIR/job_submit.sh $wrf_ntasks $((tid*$wrf_ntasks)) $HOSTPPN ./wrf.exe >& wrf.log &
  srun -N 1 -n 68 --ntasks-per-node=68 --cpus-per-task=1 --cpu_bind=cores wrf.exe >& wrf.log &
  tid=$(($tid+1))
  if [[ $tid -ge $ntid ]]; then
    wait
    tid=0
  fi
  sleep 3

  cd ..

done
wait

echo Finished running WRF ensemble on `date`
