#!/bin/bash

source ~/.bashrc
## setting
ENS_LIST=`seq 1 30`

# reading config file
export CONFIG_FILE=PATH_TO_YOUR_CONFIG
. $CONFIG_FILE
cd $SCRIPT_DIR
. util.sh

if [[ ! -d $WORK_DIR ]]; then mkdir -p $WORK_DIR; fi
cd $WORK_DIR

####total_ntasks####
if [[ $HOSTTYPE == "stampede" ]]; then
  export total_ntasks=$SLURM_NTASKS
fi
if [[ $HOSTTYPE == "jet" ]]; then
  export total_ntasks=$PBS_NP
fi

#################################
#date

export DATE=$DATE_START
export PREVDATE=$DATE_START
export NEXTDATE=$DATE_START

while [[ $NEXTDATE -le $DATE_CYCLE_END ]]; do

  export OWMAX=`min ${OBS_WIN_MAX[@]}`
  export OWMIN=`min ${OBS_WIN_MIN[@]}`
  export DT=$TIME_STEP
  export MPS=`min ${MINUTE_PER_SLOT[@]}`
  export FCSTM=`min ${FORECAST_MINUTES[@]}`
  export CP=`min ${CYCLE_PERIOD[@]}`
  if [[ $DATE == $DATE_START ]]; then
    export CP=`diff_time $DATE $DATE_CYCLE_START`
  fi

#calculate start_date and run_minutes, used by namelist_wrf.sh to generate correct
#time in namelist.input
  if $RUN_4DVAR; then
    export start_date_cycle=$DATE
    export run_minutes_cycle=`echo $CP+$OWMAX |bc`
    if [[ $DATE -ge $DATE_CYCLE_START ]]; then
      export start_date_cycle=`advance_time $start_date_cycle $OWMIN`
      export run_minutes_cycle=`echo $run_minutes+$OWMIN |bc`
    fi
  else
    export start_date_cycle=$DATE
    export run_minutes_cycle=$CP
  fi
  if $RUN_DETERMINISTIC; then
    export run_minutes_forecast=`diff_time $DATE $DATE_END`
  else
    export run_minutes_forecast=`max $CP $FCSTM`
  fi

#LBDATE
  export minute_off=`echo "(${start_date_cycle:8:2}*60+${start_date_cycle:10:2})%$LBC_INTERVAL" |bc`
  if [[ $DATE == `advance_time $start_date_cycle -$minute_off` ]]; then
    export LBDATE=`advance_time $start_date_cycle -$minute_off`
    #export LBDATE=$DATE_START
  fi
  export NEXTDATE=`advance_time $DATE $CP`
  echo "----------------------------------------------------------------------"
  echo "CURRENT CYCLE: `wrf_time_string $DATE` => `wrf_time_string $NEXTDATE`"

  mkdir -p {run,rc,fc,output,obs}/$DATE

#clear error tags
  for d in `ls run/$DATE/`; do
    if [[ `cat run/$DATE/$d/stat` != "complete" ]]; then
      echo waiting > run/$DATE/$d/stat
    fi
  done

  ######
  ## run components
  rundir=$WORK_DIR/run/$DATE/wrf_ens
  if [[ ! -d $rundir ]]; then mkdir -p $rundir; echo waiting > $rundir/stat; fi
  
  cd $rundir
  if [[ `cat stat` == "complete" ]]; then 
     echo 'going to next step'
  else
 
     #Check dependency
     wait_for_module ../update_bc ../icbc
     if [[ $DATE -gt $DATE_START ]]; then
       wait_for_module ../enkf
     fi
   
     #Setup for wrf run
     echo running > stat
     echo "  Submitting WRF ensemble jobs..."
     
     for NE in $ENS_LIST; do
       id=`expr $NE + 1000 |cut -c2-`
       if [[ ! -d $id ]]; then mkdir $id; fi
       touch $id/rsl.error.0000
       if [[ `tail -n15 $id/rsl.error.0000 |grep SUCCESS` ]]; then continue; fi
     
       cd $id
       lfs setstripe -c 1 $rundir/$id
     
       for n in `seq 1 $MAX_DOM`; do
         dm=d`expr $n + 100 |cut -c2-`
         ln -fs $WORK_DIR/fc/$DATE/wrfinput_${dm}_$id wrfinput_$dm
       done
       ln -fs $WORK_DIR/fc/$DATE/wrfbdy_d01_$id wrfbdy_d01
     
       if [[ $SST_UPDATE == 1 ]]; then
         ln -fs $WORK_DIR/rc/$LBDATE/wrflowinp_d?? .
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
       $SCRIPT_DIR/namelist_wrf.sh wrfw $RUN_DOMAIN > namelist.input
       $SCRIPT_DIR/prepare_job_script.sh $id > run_wrf
     
       echo "    submit ens #$id"
       sbatch run_wrf
       cd ..
     done
   
     date
   
   ##check errors  
   #  for d in `ls -t run/$DATE/`; do
   #    if [[ `cat run/$DATE/$d/stat` == "error" ]]; then
   #      echo CYCLING STOP DUE TO FAILED COMPONENT: $d
   #      exit 1
   #    fi
   #  done

  fi

#advance to next cycle
  cd $WORK_DIR
  export PREVDATE=$DATE
  export DATE=$NEXTDATE

done
echo CYCLING COMPLETE

