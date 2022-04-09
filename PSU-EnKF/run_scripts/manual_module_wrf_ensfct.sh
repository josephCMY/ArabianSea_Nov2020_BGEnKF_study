#!/bin/bash

## setting
export CONFIG_FILE=PATH_TO_CONFIG
export DATE_FCT_START=
export DATE_FCT_END=
export DATE_DIR=
export WRF_STORM="vortex" or "preset" or "none"
export RESTART=true or false
export NUM_ENS_FCT=ENSEMBLE_SIZE

# reading config file
. $CONFIG_FILE
cd $SCRIPT_DIR
. util.sh
export DT=$TIME_STEP


rundir=$WORK_DIR/run/$DATE_DIR/wrf_ens_fct
if [[ ! -d $rundir ]]; then mkdir -p $rundir; echo waiting > $rundir/stat; fi

cd $rundir

#Setup for wrf run
echo running > stat
echo "  Submitting WRF ensemble jobs..."

for NE in `seq 1 $NUM_ENS_FCT`; do
  id=`expr $NE + 1000 |cut -c2-`
  if [[ ! -d $id ]]; then mkdir $id; fi
  touch $id/rsl.error.0000
  if [[ `tail -n2 $id/rsl.error.0000 |grep SUCCESS` ]]; then continue; fi

  cd $id
  lfs setstripe -c 1 $rundir/$id

  for n in `seq 1 $MAX_DOM`; do
    dm=d`expr $n + 100 |cut -c2-`
    ln -fs $WORK_DIR/fc/$DATE_DIR/wrfinput_${dm}_$id wrfinput_$dm
  done
  ln -fs $WORK_DIR/fc/$DATE_DIR/wrfbdy_d01_$id wrfbdy_d01

  if $FOLLOW_STORM; then
      i_parent_start="1 "
      j_parent_start="1 "
      for n in `seq 2 $MAX_DOM`; do
        dm=d`expr $n + 100 |cut -c2-`
        outfile=wrfinput_${dm}
        i_parent_start="$i_parent_start $(ncdump -h $outfile |grep :I_PARENT_START |awk '{print $3}')"
        j_parent_start="$j_parent_start $(ncdump -h $outfile |grep :J_PARENT_START |awk '{print $3}')"
      done
      echo $i_parent_start > ij_parent_start
      echo $j_parent_start >> ij_parent_start
  fi

  rm -fr domain_moves
  if [[ $WRF_STORM == "vortex" ]]; then
    ln -fs $WRF_VORTEX_DIR/run/* .
  elif [[ $WRF_STORM == "preset" ]]; then
    ln -fs $WRF_PRESET_DIR/run/* .
    $SCRIPT_DIR/calc_domain_moves.sh $DATE_FCT_START $DATE_FCT_END domain_moves >& follow_storm.log
  else
    ln -fs $WRF_DIR/run/* .
  fi
  rm -f namelist.*

  if [[ $SST_UPDATE == 1 ]]; then
    ln -fs $WORK_DIR/rc/$DATE_DIR/wrflowinp_d?? .
  fi

  export start_date=$DATE_FCT_START
  export run_minutes=`diff_time $DATE_FCT_START $DATE_FCT_END`
  #export inputout_interval=$run_minutes
  #export inputout_begin=0
  #export inputout_end=$run_minutes
  export GET_PHYS_FROM_FILE=false
  $SCRIPT_DIR/namelist_wrf.sh wrf_fct $RUN_DOMAIN > namelist.input
  $SCRIPT_DIR/prepare_job_script.sh $id > run_wrf

  echo "    submit ens #$id"
  sbatch run_wrf
  cd ..
done

