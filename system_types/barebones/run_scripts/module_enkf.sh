#!/bin/bash
. $CONFIG_FILE
rundir=$WORK_DIR/run/$DATE/enkf
if [[ ! -d $rundir ]]; then mkdir -p $rundir; echo waiting > $rundir/stat; fi

if [[ $RUN_ENKF == 'false' ]]; then 
  echo complete > $rundir/stat
  cd $WORK_DIR/fc/$DATE
  domlist=`seq 1 $MAX_DOM`
  for n in $domlist; do
    dm=d`expr $n + 100 |cut -c2-`
    ln -sf $WORK_DIR/fc/$PREVDATE/wrfinput_${dm}_`wrf_time_string $DATE` .
    for NE in `seq 1 $NUM_ENS`; do
      id=`expr $NE + 1000 |cut -c2-`
      ln -sf $WORK_DIR/fc/$PREVDATE/wrfinput_${dm}_`wrf_time_string $DATE`_$id wrfinput_${dm}_$id
    done
  done
fi


cd $rundir
if [[ `cat stat` == "complete" ]]; then exit; fi

#Check dependency
wait_for_module ../../$PREVDATE/wrf_ens ../obsproc ../icbc
if [[ $DATE == $LBDATE ]]; then
  wait_for_module ../perturb_ic
fi

#Run EnKF
echo running > stat
echo "  Running EnKF..."

domlist=`seq 1 $MAX_DOM`
#link files
for n in $domlist; do
  dm=d`expr $n + 100 |cut -c2-`
  if [[ ! -d $dm ]]; then mkdir -p $dm; fi
  if [ -f $dm/${DATE}.finish_flag ]; then continue; fi
  cd $dm
  lfs setstripe -c 8 $rundir/$dm

  echo "    Linking files for domain $dm"


  for NE in `seq 1 $NUM_ENS`; do
    id=`expr $NE + 1000 |cut -c2-`
    ln -fs $WORK_DIR/fc/$PREVDATE/wrfinput_${dm}_`wrf_time_string $DATE`_$id $WORK_DIR/fc/$DATE/wrf_enkf_input_${dm}_$id

    # If USE_ESTIMATE_INF is true
    if $USE_ESTIMATE_INF; then
      cp -L $WORK_DIR/fc/$PREVDATE/wrfinput_${dm}_`wrf_time_string $DATE`_$id fort.`expr 80010 + $NE` 

    ## If using ABEI (fast MPI doesnt need copying)
    #elif $USE_ABEI; then
    #  cp -L $WORK_DIR/fc/$PREVDATE/wrfinput_${dm}_`wrf_time_string $DATE`_$id fort.`expr 80010 + $NE` 

    else
      ln -fs $WORK_DIR/fc/$PREVDATE/wrfinput_${dm}_`wrf_time_string $DATE`_$id fort.`expr 80010 + $NE`
    fi

    # Preparing posterior file
    cp -L fort.`expr 80010 + $NE` fort.`expr 90010 + $NE` >> link.log 2>&1 

  done
  wait


#  ncea `ls fort.800*` fort.`expr 80011 + $NE`
  cp -L fort.`expr 80010 + $NUM_ENS` fort.`expr 80011 + $NUM_ENS` 
  cp -L fort.`expr 80011 + $NUM_ENS` fort.`expr 90011 + $NUM_ENS`
  cp -L fort.`expr 80011 + $NUM_ENS` fort.`expr 60011 + $NUM_ENS`
  wait


#  echo $PREVDATE
#  echo $DATE
#  echo $WORK_DIR/fc/$PREVDATE/wrfinput_${dm}_`wrf_time_string $DATE`
#
#  cp -L $WORK_DIR/fc/$PREVDATE/wrfinput_${dm}_`wrf_time_string $DATE` fort.`expr 80011 + $NUM_ENS`
#  cp -L fort.`expr 80011 + $NUM_ENS` fort.`expr 90011 + $NUM_ENS` 
#  cp -L fort.`expr 80011 + $NUM_ENS` fort.`expr 60011 + $NUM_ENS`

  ln -fs $WRF_DIR/run/LANDUSE.TBL .
  # coefficients for CRTM
  ln -fs $CRTM_DIR/coefficients .
  # Empirical Localization Function
  #ln -fs $ELF_DIR/elf .
  #Observations
  #LITTLE_R format from obsproc
  ln -fs $OBS_DIR/obs_gts_`wrf_time_string $DATE`.3DVAR obs_3dvar_${DATE}00
#  if [ -f $WORK_DIR/obs/${DATE}/obs_gts_`wrf_time_string $DATE`.3DVAR ]; then
#     ln -fs $WORK_DIR/obs/${DATE}/obs_gts_`wrf_time_string $DATE`.3DVAR obs_3dvar_${DATE}00
#  fi
  #airborne radar superobs
  ln -fs $DATA_DIR/so/${DATE:0:4}/${DATE}_all.so_ass airborne_${DATE}_so

  ln -fs $RADIANCE_DIR/Tb_${dm}_${DATE}_so radiance_${DATE}_so

#  # If radiance bias correction enabled
#  if $RADIANCE_NBC; then
#    raw_dir=$WORK_DIR/run/$PREVDATE/enkf/$dm
#    out_dir=$WORK_DIR/run/$DATE/enkf/$dm
#    python $SCRIPT_DIR/radiance_bc_regress_only.py $raw_dir
#    cp $raw_dir/radiance_bc_coeffs.txt $out_dir
#  fi


  ln -fs $ENKF_DIR/enkf.mpi .
#  ln -fs $ENKF_DIR/enkf_abei.mpi .
#  ln -fs $ENKF_DIR/enkf_after_abei.mpi .

#  # updating non-Q variables only every 1-hour
#  if [[ ${DATE:10:2} == '00'  ]]; then
#    ln -fs $ENKF_DIR/enkf.mpi .
#  else
#    #ln -fs $ENKF_DIR/enkf.mpi .
#    ln -fs $ENKF_DIR/enkf_hydro.mpi enkf.mpi
#  fi

  # multiplicative inflation
  if $USE_ESTIMATE_INF; then
    ln -fs $ENKF_DIR/cal_inflate.mpi .
    if [[ $PREVDATE -ge $DATE_CYCLE_START ]]; then
      cp $WORK_DIR/run/$PREVDATE/enkf/$dm/parameters_update${PREVDATE} parameters_update
    fi
  fi

  if $USE_ABEI && [[ $PREVDATE -ge $DATE_CYCLE_START ]]; then
    ln -sf $WORK_DIR/run/$PREVDATE/enkf/$dm/parameters_update${PREVDATE} parameters_update
  fi
#  if $USE_ABEI && [[ $PREVDATE -lt $DATE_CYCLE_START ]]; then
#    USE_ABEI=false
#    USE_RADIANCE=false
#    abei_temp_flag=true
#  fi

  
#  echo $USE_ABEI

#  $SCRIPT_DIR/namelist_enkf.sh $n $USE_ABEI $USE_RADIANCE > namelist.enkf
  $SCRIPT_DIR/namelist_enkf.sh $n  > namelist.enkf

#  if $abei_temp_flag; then
#    USE_ABEI=true
#    USE_RADIANCE=true
#    abei_temp_flag=false
#  fi

  cd ..
done


# replacing prior mean with determinstic forecast
if [ $DATE != $DATE_CYCLE_START ]; then
  if $REPLACE_MEAN; then
   if [[ $REPLACE_MEAN_WITH == "prior_forecast" ]]; then
    tid=0
    nn=$((($enkf_ntasks+$HOSTPPN-$enkf_ntasks%$HOSTPPN)/$HOSTPPN))
    nt=$(($total_ntasks/$HOSTPPN/$nn))
    for n in $domlist; do
     dm=d`expr $n + 100 |cut -c2-`
     cd $dm
     if [[ ! -d replace_mean ]]; then mkdir -p replace_mean; fi
      cd replace_mean
      echo "  Replacing ens mean with $REPLACE_MEAN_WITH for domain $dm"
      for NE in `seq 1 $((NUM_ENS+1))`; do
        mv ../fort.`expr 80010 + $NE` fort.`expr 80010 + $NE`
        cp -L fort.`expr 80010 + $NE` fort.`expr 90010 + $NE`
      done
      if [[ $REPLACE_MEAN_WITH == "prior_forecast" ]]; then
        ln -fs $WORK_DIR/fc/$PREVDATE/wrfinput_${dm}_`wrf_time_string $DATE` fort.70010
      fi
      ln -fs $ENKF_DIR/replace_mean.exe .
      export SLURM_TASKS_PER_NODE="$HOSTPPN(x$SLURM_NNODES)"
      ibrun -n $enkf_ntasks -o $((tid*$enkf_ntasks)) ./replace_mean.exe $NUM_ENS >& replace_mean.log &
      export SLURM_TASKS_PER_NODE="$((SLURM_NTASKS/$SLURM_NNODES))(x$SLURM_NNODES)"
      tid=$((tid+1))
      if [[ $tid == $nt ]]; then
        tid=0
        wait
      fi
      cd ../..
    done
  
    for n in $domlist; do
      dm=d`expr $n + 100 |cut -c2-`
      cd $dm/replace_mean/
      watch_log replace_mean.log Successful 1 $rundir
      for NE in `seq 1 $((NUM_ENS+1))`; do
        mv fort.`expr 90010 + $NE` ../fort.`expr 80010 + $NE`
        cp ../fort.`expr 80010 + $NE` ../fort.`expr 90010 + $NE`
      done
      cd ../..
    done
   fi
  fi
fi

wait

# inflate prior members
if $USE_ESTIMATE_INF; then
  if [[ ${DATE:10:2} == '00'  ]]; then
    tid=0
    nn=$((($enkf_ntasks+$enkf_ppn-$enkf_ntasks%$enkf_ppn)/$enkf_ppn))
    nt=$(($total_ntasks/$HOSTPPN/$nn))
    for n in $domlist; do
      dm=d`expr $n + 100 |cut -c2-`
      if [ -f $dm/${DATE}.inffin_flag ]; then continue; fi
      cd $dm
      echo "    Running inflation for domain $dm"
      wait
      $SCRIPT_DIR/job_submit.sh $enkf_ntasks $((tid*$enkf_ntasks)) $enkf_ppn ./cal_inflate.mpi >& inflate.log &
      tid=$((tid+1))
      if [[ $tid == $nt ]]; then
        tid=0
        wait
      fi
      cd ..
    done
    wait
    #Check output
    for n in $domlist; do
      dm=d`expr $n + 100 |cut -c2-`
    #  watch_log $dm/enkf.log Successful 5 $rundir
      watch_log $dm/${DATE}.inffin_flag _ 5 $rundir
    done
  fi
fi

#run enkf.mpi
tid=0
nn=$((($enkf_ntasks+$enkf_ppn-$enkf_ntasks%$enkf_ppn)/$enkf_ppn))
nt=$(($total_ntasks/$HOSTPPN/$nn))
for n in $domlist; do
  dm=d`expr $n + 100 |cut -c2-`
  if [ -f $dm/${DATE}.finish_flag ]; then continue; fi
  cd $dm
  echo "    Running enkf.mpi for domain $dm"
  $SCRIPT_DIR/job_submit.sh $enkf_ntasks $((tid*$enkf_ntasks)) $enkf_ppn ./enkf.mpi >& enkf.log #&

#  if $USE_ABEI && [[ $PREVDATE -ge $DATE_CYCLE_START ]]; then
#    $SCRIPT_DIR/job_submit.sh $enkf_ntasks $((tid*$enkf_ntasks)) $enkf_ppn ./enkf_abei.mpi >& enkf_abei.log
#    $SCRIPT_DIR/job_submit.sh $enkf_ntasks $((tid*$enkf_ntasks)) $enkf_ppn ./enkf_after_abei.mpi >& enkf.log &
#  else
#    $SCRIPT_DIR/job_submit.sh $enkf_ntasks $((tid*$enkf_ntasks)) $enkf_ppn ./enkf.mpi >& enkf.log &
#  fi
#  tid=$((tid+1))
#  if [[ $tid == $nt ]]; then
#    wait
#    tid=0
#  fi
  cd ..
done

#Check output
for n in $domlist; do
  dm=d`expr $n + 100 |cut -c2-`
#  watch_log $dm/enkf.log Successful 5 $rundir
  watch_log $dm/${DATE}.finish_flag _ 5 $rundir
done

# replacing mean with first guess (GFS/FNL) reanalysis
if $REPLACE_MEAN; then
 if [[ $REPLACE_MEAN_WITH != "prior_forecast" ]]; then
  tid=0
  nn=$((($enkf_ntasks+$HOSTPPN-$enkf_ntasks%$HOSTPPN)/$HOSTPPN))
  nt=$(($total_ntasks/$HOSTPPN/$nn))
  for n in $domlist; do
    dm=d`expr $n + 100 |cut -c2-`
    cd $dm
    if [[ ! -d replace_mean ]]; then mkdir -p replace_mean; fi
    cd replace_mean
    echo "  Replacing ens mean with $REPLACE_MEAN_WITH for domain $dm"
    for NE in `seq 1 $((NUM_ENS+1))`; do
      mv ../fort.`expr 90010 + $NE` fort.`expr 80010 + $NE`
      cp fort.`expr 80010 + $NE` fort.`expr 90010 + $NE`
    done
    if [[ $REPLACE_MEAN_WITH == "forecast" ]]; then
      ln -fs $WORK_DIR/fc/$PREVDATE/wrfinput_${dm}_`wrf_time_string $DATE` fort.70010
    elif [[ $REPLACE_MEAN_WITH == "gfs" ]]; then
      ln -fs $WORK_DIR/rc/$DATE/wrfinput_$dm fort.70010
    fi
    ln -fs $ENKF_DIR/replace_mean.exe .
    export SLURM_TASKS_PER_NODE="$HOSTPPN(x$SLURM_NNODES)"
    ibrun -n $enkf_ntasks -o $((tid*$enkf_ntasks)) ./replace_mean.exe $NUM_ENS >& replace_mean.log &
    export SLURM_TASKS_PER_NODE="$((SLURM_NTASKS/$SLURM_NNODES))(x$SLURM_NNODES)"
    tid=$((tid+1))
    if [[ $tid == $nt ]]; then
      tid=0
      wait
    fi
    cd ../..
  done
  for n in $domlist; do
    dm=d`expr $n + 100 |cut -c2-`
    cd $dm/replace_mean/
    watch_log replace_mean.log Successful 1 $rundir
    for NE in `seq 1 $((NUM_ENS+1))`; do
      mv fort.`expr 90010 + $NE` ../
    done
    cd ..
    cd ..
  done
 fi
fi


### replacing outside with first guess (GFS/FNL) reanalysis
if $REPLACE_ENVIRONMENT; then
 nt=$((total_ntasks/$HOSTPPN))
 if [ $DATE == $LBDATE ]; then
  tid=0
  for n in $domlist; do
#  for n in `seq 1 $((MAX_DOM-1))`; do
    dm=d`expr $n + 100 |cut -c2-`
    cd $dm
    if [[ ! -d replace_environment ]]; then mkdir -p replace_environment; fi
    cd replace_environment
    echo "  Replacing environment with GFS for domain $dm"
    for NE in `seq 1 $((NUM_ENS+1))`; do
      id=`expr $NE + 1000 |cut -c2-`
      if [[ ! -d $id ]]; then mkdir $id; fi
      if [[ `tail -n2 ${id}/replace_environment.log |grep Successful` ]]; then continue; fi
      cd $id

      ln -sf $ENKF_DIR/replace_environment_by_gfs.exe .
      ln -sf $TCVITALS_DIR/${DATE:0:4}/${DATE}.${STORM_ID}-tcvitals.dat tcvitals.dat
      mv ../../fort.`expr 90010 + $NE` wrfinput
      if  [[ $NE == `expr $((NUM_ENS+1))` ]]; then
        ln -sf $WORK_DIR/rc/$DATE/wrfinput_${dm} wrfinput_gfs
      else
        ln -sf $WORK_DIR/rc/$DATE/wrfinput_${dm}_$id wrfinput_gfs
      fi
      ./replace_environment_by_gfs.exe >& replace_environment.log
      tid=$((tid+1))
      if [[ $tid == $nt ]]; then
        tid=0
        wait
      fi
      cd ..
    done
    cd ../..
  done

  for n in $domlist; do
#  for n in `seq 1 $((MAX_DOM-1))`; do
    dm=d`expr $n + 100 |cut -c2-`
    cd $dm
    for NE in `seq 1 $((NUM_ENS+1))`; do
      id=`expr $NE + 1000 |cut -c2-`
      watch_log replace_environment/$id/replace_environment.log Successful 1 $rundir
      mv replace_environment/$id/wrfinput fort.`expr 90010 + $NE` 
    done
    cd ..
  done
 fi
fi
###

for n in $domlist; do
  dm=d`expr $n + 100 |cut -c2-`
  for NE in `seq 1 $NUM_ENS`; do
    id=`expr $NE + 1000 |cut -c2-`
    ln -sf $rundir/$dm/fort.`expr 90010 + $NE` $WORK_DIR/fc/$DATE/wrf_enkf_output_${dm}_$id
#    cp $dm/fort.`expr 90010 + $NE` $WORK_DIR/fc/$DATE/wrf_enkf_output_${dm}_$id
    ln -fs $WORK_DIR/fc/$DATE/wrf_enkf_output_${dm}_$id $WORK_DIR/fc/$DATE/wrfinput_${dm}_$id
    if [ ! -f $WORK_DIR/fc/$DATE/wrf_enkf_input_${dm}_$id ]; then
      ln -fs $WORK_DIR/fc/$PREVDATE/wrfinput_${dm}_`wrf_time_string $DATE`_$id $WORK_DIR/fc/$DATE/wrf_enkf_input_${dm}_$id
    fi
  done
  cp $dm/fort.`expr 80011 + $NUM_ENS` $WORK_DIR/fc/$DATE/wrf_enkf_input_${dm}_mean
  cp $dm/fort.`expr 90011 + $NUM_ENS` $WORK_DIR/fc/$DATE/wrf_enkf_output_${dm}_mean
  ln -fs $WORK_DIR/fc/$DATE/wrf_enkf_output_${dm}_mean $WORK_DIR/fc/$DATE/wrfinput_${dm}
done


echo complete > stat

