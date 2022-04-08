#!/bin/bash
. $CONFIG_FILE
rundir=$WORK_DIR/run/$DATE/wrf_ens
if [[ ! -d $rundir ]]; then mkdir -p $rundir; echo waiting > $rundir/stat; fi

cd $rundir
if [[ `cat stat` == "complete" ]]; then exit; fi

#Check dependency
#wait_for_module ../update_bc ../icbc 
#if [[ $DATE -gt $DATE_START ]]; then
#  wait_for_module ../enkf
#fi

#Setup for wrf run
echo running > stat
date
echo "  Preparing to run WRF ensemble..."

# Generate the things needed to run the wrf ens
tid=1  #does not start from 0, because the wrf forecast runs with ens at the same time.
if [[ `cat ../wrf/stat` == "complete" ]]; then tid=0; fi
#if [[ `tail -n2 ../wrf/rsl.error.0000 |grep SUCCESS` ]]; then tid=0; fi
nt=$((total_ntasks/$wrf_ntasks))
tid=1
if [[ `cat ../wrf/stat` == "complete" ]]; then tid=0; fi
#  if [[ `tail -n2 ../wrf/rsl.error.0000 |grep SUCCESS` ]]; then tid=0; fi
. $CONFIG_FILE
for NE in `seq -f "%03g" 1 $NUM_ENS`; do
    id=$NE    #`expr $NE + 1000 |cut -c2-`
    if [[ ! -d $id ]]; then mkdir $id; fi
    touch $id/rsl.error.0000
    if [[ `tail -n2 $id/rsl.error.0000 |grep SUCCESS` ]]; then continue; fi

    cd $id
    lfs setstripe -c 1 $rundir/$id

    for n in `seq 1 $MAX_DOM`; do
      dm=d`expr $n + 100 |cut -c2-`
      ln -fs $WORK_DIR/fc/$DATE/wrfinput_${dm}_$id wrfinput_$dm
    done
    ln -sf $BDY_DIR/$id/wrfbdy_d01 .

    if [[ $SST_UPDATE == 1 ]]; then
      if [[ $NE == 001 ]]; then
        cp $WORK_DIR/DA/wrflowinp_interpolate.py .
        python3 wrflowinp_interpolate.py $BDY_DIR/$id/wrflowinp_d01 $DATE wrflowinp_d01
      else
        ln -s ../001/wrflowinp_d01
      fi
    fi

    if $FOLLOW_STORM; then
      cp $WORK_DIR/rc/$DATE/ij_parent_start .
      cp $WORK_DIR/rc/$DATE/domain_moves .
      ln -fs $WRF_PRESET_DIR/run/* .
    else
      ln -fs $WRF_DIR/run/* .
    fi
    rm -f namelist.*

    export start_date=$start_date_cycle
    export run_minutes=$run_minutes_cycle 
    export inputout_interval=$run_minutes
    export inputout_begin=0
    export inputout_end=$run_minutes
    export GET_PHYS_FROM_FILE=false
    if [[ $NE == 001 ]]; then
      $SCRIPT_DIR/namelist_wrf.sh wrf $RUN_DOMAIN > namelist.input
    else
      ln -sf ../001/namelist.input namelist.input
    fi

    cd ..
done
#wait

date

# RUN WRF
echo "  Proceeding to run WRF ensemble"
tid=1  #does not start from 0, because the wrf forecast runs with ens at the same time.
if [[ `cat ../wrf/stat` == "complete" ]]; then tid=0; fi
#if [[ `tail -n2 ../wrf/rsl.error.0000 |grep SUCCESS` ]]; then tid=0; fi
nt=$((total_ntasks/$wrf_ntasks))
. $CONFIG_FILE
for rerun in `seq 1 3`; do

  # Special job submission script to run all WRF members in one shot
  if [[ $HOSTTYPE == 'cori-login' ]]; then

    echo Rerun $rerun Submitting $NUM_ENS ensemble runs to interactive node
    echo Submission date: `date`

    wrf_runtime=$(($CYCLE_PERIOD*2))
    if [[ $wrf_runtime -gt 120 ]]; then
      wrf_runtime=120
    fi

    salloc -C knl -N 25 -q interactive -t $wrf_runtime --job-name='wrf-ens' $SCRIPT_DIR/job_submit_wrf_ens.sh

  # Special job submission script to run all WRF members in one shot
  elif [[ $HOSTTYPE == 'cori-regular' ]]; then

    echo Rerun $rerun Submitting $NUM_ENS ensemble runs to regular node
    echo Submission date: `date`

    wrf_runtime=$(($CYCLE_PERIOD*2))
    if [[ $wrf_runtime -gt 120 ]]; then
      wrf_runtime=120
    fi

    salloc -C knl -N 25 -q regular -t $wrf_runtime --job-name='wrf-ens' $SCRIPT_DIR/job_submit_wrf_ens.sh



  # If not in login node
  else

    for NE in `seq -f "%03g" 1 $NUM_ENS`; do
      id=$NE  #`expr $NE + 1000 |cut -c2-`
      if [[ ! -d $id ]]; then mkdir $id; fi
      touch $id/rsl.error.0000
      if [[ `tail -n 20 $id/rsl.error.0000 |grep SUCCESS` ]]; then 
        echo "  rerun $rerun, id $id, `tail -n 20 $id/rsl.error.0000 |grep SUCCESS`"
        continue
      fi
      echo "  rerun $rerun, id $id, RUNNING WRF"

      cd $id
      lfs setstripe -c 1 $rundir/$id


      # Check if wrflowinp has been constructed
      # The wrflowinp construction step was thrown into the background earlier. 
      while [[ ! -f wrflowinp_d01 ]]; do
        sleep 5
      done

      # Run WRF
      $SCRIPT_DIR/job_submit.sh $wrf_ntasks $((tid*$wrf_ntasks)) $HOSTPPN ./wrf.exe >& wrf.log &
      sleep 3

      tid=$((tid+1))
      if [[ $tid == $nt ]]; then

         
        wait
        # Reset tid
        if [[ `cat ../../wrf/stat` == "complete" ]]; then
          tid=0
        else
          tid=1
        fi

        
      fi
      cd ..
    done
  wait
  fi
done
wait

cd $WORK_DIR/fc/$DATE

for NE in `seq -f "%03g" 1 $NUM_ENS`; do
  id=$NE   #`expr $NE + 1000 |cut -c2-`
#  watch_log $id/rsl.error.0000 SUCCESS 1 $rundir
  outfile=$id/wrfinput_d01_`wrf_time_string $NEXTDATE`
#  watch_file $outfile 1 $rundir
  ln -fs $rundir/$outfile wrfinput_d01_`wrf_time_string $NEXTDATE`_$id
#  ln -s $outfile $WORK_DIR/fc/$DATE/wrfinput_d01_`wrf_time_string $NEXTDATE`_$id
  if [ $MAX_DOM -gt 1 ]; then
    for n in `seq 2 $MAX_DOM`; do
      dm=d`expr $n + 100 |cut -c2-`
      outfile=$id/wrfout_${dm}_`wrf_time_string $NEXTDATE`
      watch_file $outfile 1 $rundir
      ln -fs $rundir/outfile wrfinput_${dm}_`wrf_time_string $NEXTDATE`_$id
#      ln -fs $rundir/$outfile $WORK_DIR/fc/$DATE/wrfinput_${dm}_`wrf_time_string $NEXTDATE`_$id
    done
  fi
done

cd $rundir
date
#if $CLEAN; then rm $rundir/$id/wrfout* ; fi
echo complete > stat
