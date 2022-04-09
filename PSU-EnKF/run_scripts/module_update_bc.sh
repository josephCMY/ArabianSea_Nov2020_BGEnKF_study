#!/bin/bash
. $CONFIG_FILE

rundir=$WORK_DIR/run/$DATE/update_bc
if [[ ! -d $rundir ]]; then mkdir -p $rundir; echo running > $rundir/stat; fi
cd $rundir
if [[ `cat stat` == "complete" ]]; then exit; fi

#Check dependency
#----------------
if $RUN_ENKF; then
  if [[ $DATE == $DATE_START ]]; then
    wait_for_module ../perturb_ic ../icbc
  else
    wait_for_module ../enkf
  fi
fi

echo running > stat
echo "  Running UpdateBC..."


# Construct "namelist" needed to run da_update_bc.exe
# -----------------------------------------------------
cat > parame.in << EOF
&control_param
 da_file               = 'wrfinput_d01_update'
 wrf_bdy_file          = 'wrfbdy_d01_update'
 domain_id             = 1
 debug                 = .true. 
 update_lateral_bdy    = .true. 
 update_low_bdy        = .false.
 update_lsm            = .false.
 iswater               = 17
 var4d_lbc             = .false.
/
EOF



# Construct individual directories to run update_bc for each member
# ------------------------------------------------------------------
for id in `seq -f "%03g" 1 $NUM_ENS`; do

  # Construct the directory itself
  if [[ ! -d $id ]]; then mkdir $id; fi
  touch $id/update_bc.log

  ## Skip over directory if the update_bc already ran
  #if [[ `tail -n1 $id/update_bc.log |grep successfully` ]]; then continue; fi

  # Enter member's update_bc directory and link in executable and "namelist"
  cd $id
  ln -fs ../parame.in .
  ln -fs $WRFDA_DIR/var/da/da_update_bc.exe .

  # Use Joseph's interpolate_wrfbdy to interpolate wrfbdy to desired time
  ln -fs $SCRIPT_DIR/interpolate_wrfbdy.py .
  ln -fs $BDY_DIR/$id/wrfbdy_d01 wrfbdy_src
  python interpolate_wrfbdy.py wrfbdy_src $DATE $CYCLE_PERIOD wrfbdy_d01_update

  # Link in the wrfinput made by the EnKF
  ln -fs $WORK_DIR/fc/$DATE/wrfinput_d01_$id wrfinput_d01_update

  # Link in job submission script
  ln -fs $SCRIPT_DIR/job_submit.sh

  cd ..
done
#wait



# Run the da_update_bc.exe for all ensemble members
# -------------------------------------------------
# Outer rerun loop in event of random code failure (HPCs can be unstable)
for rerun in `seq 1 3`; do

  # Counter for number of actively running processes
  active_procs=0

  # Iterate over all ensemble members
  for id in `seq -f "%03g" 1 $NUM_ENS`; do

    # Skip over if da_update_bc already ran successfully
    if [[ `tail -n 10 $id/update_bc.log |grep successfully` ]]; then continue; fi

    echo Rerun $rerun -- running UpdateBC for member $id

    cd $id

    # If da_update_bc hasn't been successfully run, run da_update_bc using 
    # HPC-appropriate commands (job_submit.sh)
    ./job_submit.sh 1 $active_procs $HOSTPPN ./da_update_bc.exe >& update_bc.log &

    cd ..

    # Increment number of processes just to keep track
    active_procs=$(( $active_procs + 1 ))

  done # ----- End of loop over members

  # Wait for all background processes (ie, all instances of da_update_bc) to clear.
  wait

done # --- End of rerun loop




# Check if all update_bcs ran successfully and also symlink the updated bdy files 
# to the FC directory
flag=1
for id in `seq -f "%03g" 1 $NUM_ENS`; do

   if [[ `tail -n 10 $id/update_bc.log |grep successfully` ]]; then
     # Symlink in appropriate directory
     cd $WORK_DIR/fc/$DATE
     ln -sf $rundir/$id/wrfbdy_d01_update wrfbdy_d01_update_$id
     cd $rundir
   else
     # If missing, issue print statement
     echo "  "da_update_bc failed for member $id ! Please check.
     flag=-1
   fi

done

# If all update_bcs ran, move along to next item.
if [[ $flag -eq 1 ]]; then
  echo complete > stat
fi

