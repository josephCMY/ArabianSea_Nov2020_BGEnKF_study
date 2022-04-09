#!/bin/bash

#####header for Cori######
#SBATCH -J SBATCH_EXP_NAME 
#SBATCH -q SBATCH_QUEUE
#SBATCH -C knl
#SBATCH -N SBATCH_NNODES
#SBATCH -n SBATCH_NTASKS
#SBATCH --time=SBATCH_TIME
#SBATCH -o log.SBATCH_LOG
#SBATCH -e err.SBATCH_ERR

# OpenMP settings
export OMP_NUM_THREADS=1
export OMP_PLACES=threads
export OMP_PROC_BIND=spread


source ~/.bashrc


#load configuration files, functions, parameters
export CONFIG_FILE=CONFIG_FILE_PATH

. $CONFIG_FILE
cd $SCRIPT_DIR
. util.sh

if [[ ! -d $WORK_DIR ]]; then mkdir -p $WORK_DIR; fi

cd $WORK_DIR
# Forcing system to skip first day stuff of the run is completely fresh
if [[ ! -d run/201705300000 ]]; then
  mkdir -p run/201705300000
  cd run/201705300000
  mkdir icbc perturb_ic wrf wrf_ens
  for f in `ls`; do
    echo complete > $f/stat
  done
fi

cd $WORK_DIR


####total_ntasks####
if [[ $HOSTTYPE == "stampede" ]]; then
  export total_ntasks=$SLURM_NTASKS
fi
if [[ $HOSTTYPE == "jet" ]]; then
  export total_ntasks=$PBS_NP
fi
if [[ $HOSTTYPE == "cori-knl" ]]; then
  export total_ntasks=$SLURM_NTASKS
fi

echo total_ntasks is $total_ntasks
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
#run components
  #$SCRIPT_DIR/module_wps.sh &
  #$SCRIPT_DIR/module_real.sh &
  $SCRIPT_DIR/module_icbc.sh &
  if $RUN_ENKF && [ $DATE == $DATE_START ]; then
    $SCRIPT_DIR/module_perturb_ic.sh &
    $SCRIPT_DIR/module_update_bc.sh &
    $SCRIPT_DIR/module_wrf_ens.sh &
    wait
  fi
  if [ $DATE == $LBDATE ]; then
    $SCRIPT_DIR/module_perturb_ic.sh &
  fi
  if [ $DATE -ge $DATE_CYCLE_START ] && [ $DATE -le $DATE_CYCLE_END ]; then
    if $RUN_ENKF || $RUN_4DVAR; then
      $SCRIPT_DIR/module_obsproc.sh &
    fi
#    if $RUN_ENKF; then
    $SCRIPT_DIR/module_enkf.sh &
    wait
#    fi
    if $RUN_4DVAR; then
      $SCRIPT_DIR/module_4dvar.sh &
    fi
    if $RUN_ENKF || $RUN_4DVAR; then
      $SCRIPT_DIR/module_update_bc.sh &
    fi
#    if $RUN_ENKF; then
    $SCRIPT_DIR/module_wrf.sh &
    $SCRIPT_DIR/module_wrf_ens.sh &
    wait
#    fi
  fi
#  $SCRIPT_DIR/module_wrf.sh &

  wait
  date

#check errors  
  for d in `ls -t run/$DATE/`; do
    if [[ `cat run/$DATE/$d/stat` == "error" ]]; then
      echo CYCLING STOP DUE TO FAILED COMPONENT: $d
      exit 1
    fi
  done

#advance to next cycle
#  if [[ $DATE == $DATE_START ]]; then
#    export CP=`min ${CYCLE_PERIOD[@]}`
#    export NEXTDATE=`advance_time $DATE $CP`
#    export PREVDATE=$DATE
#    export DATE=$DATE
#  else
#    export PREVDATE=$DATE
#    export DATE=$NEXTDATE
#  fi

  export PREVDATE=$DATE
  export DATE=$NEXTDATE

done
echo CYCLING COMPLETE



#########################################
#echo Proceed with computing the simulated brightness temperatures
#rm -r $WORK_DIR/prior_Tb
#mkdir $WORK_DIR/prior_Tb
#cd $WORK_DIR/prior_Tb
#
#fc_dir=$WORK_DIR'/fc'
#ntasks=16
#final_id=300
#
##rm -r 2011*
#
#tid=0
#for date in `ls $fc_dir`; do
#  # Checking if CRTM has been run before
#  if [[ ! -d $date ]]; then
#    # If CRTM hasnt been run yet, folder should not exist
#    mkdir $date
#    cd $date
#    # Make links
#    ln -s /work/04920/tg842199/stampede2/mjo_dynamo_real_enkf/crtm_src/crtm.exe .
#    for ens in `seq -f "%03g" 1 50`; do
#      fname=`readlink $fc_dir/$date/wrf_enkf_input_d01_$ens`
#      echo $fname
#      ln -s $fname enkf_input_d01_$ens
#    done
#    # Run the CRTM
#    for ens in `seq -f "%03g" 1 50`; do
#      ibrun -n $ntasks -o $((tid*ntasks)) crtm.exe enkf_input_d01_$ens met7_ir_wv_"$ens".bin &
#      tid=$(($tid+1))
#      if [[ $tid == $final_id ]]; then
#        tid=0
#        wait
#      fi
#    done
#    cd ..
#  fi
#done
#wait
#
#echo Finished running CRTM
#
