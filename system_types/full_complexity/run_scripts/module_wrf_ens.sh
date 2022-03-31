#!/bin/bash --login
. $CONFIG_FILE
rundir=$WORK_DIR/run/$DATE/wrf_ens
if [[ ! -d $rundir ]]; then mkdir -p $rundir; echo waiting > $rundir/stat; fi

cd $rundir
if [[ `cat stat` == "complete" ]]; then exit; fi

#Check dependency
wait_for_module ../update_bc ../icbc 

#if [[ $DATE == $LBDATE ]]; then 
#  wait_for_module ../perturb_ic; 
#fi

#Setup for wrf run
echo running > stat
echo "  Running WRF ensemble..."

tid=0 
nt=$((total_ntasks/$wrf_ntasks))

if $RECENTER; then NUM_ENS=`expr $NUM_ENS + 1`; fi


for rerun in 1 1 1 1 1 1; do
   for NE in `seq 1 $NUM_ENS`; do
   #for NE in `seq 1 1`; do
     id=`expr $NE + 1000 |cut -c2-`
     if [[ ! -d $id ]]; then mkdir $id; fi
     touch $id/rsl.error.0000
     if [[ `tail -n1 $id/rsl.error.0000 |grep SUCCESS` ]]; then continue; fi
     if [[ `tail -n50 $id/rsl.error.0000 |grep segmentation` ]]; then
        echo Segmentation fault occurred during WRF ENS due to problem with posterior of ensemble member $NE - exiting WRF ENS now ...
        touch $WORK_DIR/run/$DATE/segmentation
        exit
     fi
     cd $id
     #ln -fs $WRF_DIR/run/* .
     ln -fs $WRF_PRESET_DIR/run/* .
     rm -f namelist.*
   
     for n in `seq 1 $MAX_DOM`; do
       dm=d`expr $n + 100 |cut -c2-`
       #if [[ $NE -eq $NUM_ENS ]]; then
       #  ln -fs $WORK_DIR/fc/$DATE/wrfinput_${dm} wrfinput_$dm
       #  ln -fs $WORK_DIR/fc/$DATE/wrfbdy_d01 wrfbdy_d01
       #else
       ln -fs $WORK_DIR/fc/$DATE/wrfinput_${dm}_$id wrfinput_$dm
       ln -fs $WORK_DIR/fc/$DATE/wrfbdy_d01_$id wrfbdy_d01
       #fi
       ncl -Q $SCRIPT_DIR/util_change_nc_time.ncl 'ncfile="wrfinput_d01"' 'time="'`wrf_time_string $DATE`'"'
     done
     echo Member $NE 

     if $FOLLOW_STORM; then
       cp $WORK_DIR/rc/$DATE/ij_parent_start .
       cp $WORK_DIR/rc/$DATE/domain_moves .
     fi
    
     if [[ $SST_UPDATE == 1 ]]; then
       ln -fs $WORK_DIR/rc/$LBDATE/wrflowinp_d?? .
     fi
   
     export start_date=$start_date_cycle
     export run_minutes=$run_minutes_cycle
     export inputout_interval=$run_minutes
     export inputout_begin=0
     export inputout_end=$run_minutes
     export GET_PHYS_FROM_FILE=false
   
     if [[ $DATE == $DATE_START ]]; then
       export wrf_for=forecast
       $SCRIPT_DIR/namelist_wrf_realtime.sh wrfw $RUN_DOMAIN > namelist.input
     else
       export wrf_for=cycle
       $SCRIPT_DIR/namelist_wrf_realtime.sh wrfw $RUN_DOMAIN > namelist.input
     fi
     
     $SCRIPT_DIR/job_submit.sh $wrf_ntasks $((tid*$wrf_ntasks)) $HOSTPPN ./wrf.exe >& wrf.log &
     tid=$((tid+1))
     if [[ $tid -ge $nt ]]; then
       tid=0
       wait
     fi
     cd ..
   done
done 
wait
	

for NE in `seq 1 $NUM_ENS`; do
  id=`expr $NE + 1000 |cut -c2-`
  watch_log $id/rsl.error.0000 SUCCESS 5 $rundir
  rm $id/rsl.*
  echo SUCCESS > $id/rsl.error.0000
  #outfile=$id/wrfinput_d01_`wrf_time_string $NEXTDATE`
  outfile=$id/wrfout_d01_`wrf_time_string $NEXTDATE`
  #outfile=$id/wrfrst_d01_`wrf_time_string $NEXTDATE`
  watch_file $outfile 5 $rundir
  if $RECENTER && [[ $NE -eq $NUM_ENS ]]; then
    mv $outfile $WORK_DIR/fc/$DATE/wrfinput_d01_`wrf_time_string $NEXTDATE`
  else
    cp $outfile $WORK_DIR/fc/$DATE/wrfinput_d01_`wrf_time_string $NEXTDATE`_$id
  fi
  if [ $MAX_DOM -gt 1 ]; then
    for n in `seq 2 $MAX_DOM`; do
      dm=d`expr $n + 100 |cut -c2-`
      outfile=$id/wrfout_${dm}_`wrf_time_string $NEXTDATE`
      #outfile=$id/wrfrst_${dm}_`wrf_time_string $NEXTDATE`
      watch_file $outfile 5 $rundir
      if $RECENTER && [[ $NE -eq $NUM_ENS ]]; then
        mv $outfile $WORK_DIR/fc/$DATE/wrfinput_${dm}_`wrf_time_string $NEXTDATE`
      else
        cp $outfile $WORK_DIR/fc/$DATE/wrfinput_${dm}_`wrf_time_string $NEXTDATE`_$id
      fi 
    done
  fi
done

if $CLEAN; then rm $rundir/$id/wrfout* ; fi
echo complete > stat
