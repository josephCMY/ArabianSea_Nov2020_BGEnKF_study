#!/bin/bash
. $CONFIG_FILE
rundir=$WORK_DIR/run/$DATE/msfat
if [[ ! -d $rundir ]]; then mkdir -p $rundir; echo waiting > $rundir/stat; fi

if [[ $RUN_MSFAT == 'false' ]]; then 
  echo complete > $rundir/stat
fi

cd $rundir
if [[ `cat stat` == "complete" ]]; then exit; fi


# Check dependency
wait_for_module ../../$PREVDATE/wrf_ens ../obsproc ../icbc
if [[ $DATE == $LBDATE ]]; then
  wait_for_module ../perturb_ic
fi

# Tell stat file that we are running FAT
echo running > stat
echo "  Running multismooth FAT..."



###############################################################
# SECTION 1: SETTING UP DIRECTORY FOR MULTISMOOTH FAT
###############################################################

domlist=`seq 1 $MAX_DOM`

# Iterating over all domains
for n in $domlist; do


  dm=d`expr $n + 100 |cut -c2-`
  if [[ ! -d $dm ]]; then mkdir -p $dm; fi
  if [ -f $dm/${DATE}.finish_flag ]; then continue; fi
  cd $dm
  #lfs setstripe -c 1 $rundir/$dm


  # Copying over ensemble files
  echo "    Linking files for domain $dm"
  for id in `seq -f "%03g" 1 $NUM_ENS`; do

    # Link in prior files
    ln -fs $WORK_DIR/fc/$PREVDATE/wrfinput_${dm}_`wrf_time_string $DATE`_$id wrfinput_xf_${dm}_$id

    # Copy in posterior files
    cp -L wrfinput_xf_${dm}_$id wrfinput_fat_${dm}_$id

  done

  # Symbolic link in CRTM program and necessary items
  ln -s $CRTM_SAT_DIR/crtm.exe
  ln -fs $WRF_DIR/run/LANDUSE.TBL .
  ln -fs $CRTM_DIR/crtm_wrf/coefficients .

  # Symbolic link in FAT code
  ln -s $SCRIPT_DIR/funclib_FAT.py .
  ln -s $SCRIPT_DIR/funclib_penalty_smoothness.py .
  ln -s $SCRIPT_DIR/funclib_penalty_residual.py .


  # Radiance observations
  ln -fs $RADIANCE_DIR/Tb_${dm}_${DATE}_so radiance_${DATE}_so

  cd ..
done



###############################################################
# SECTION 2: RUN CRTM ON THE ENSEMBLE FILE
###############################################################
ntid=$(($total_ntasks/$crtm_ntasks))
tid=0
for id in `seq -f "%03g" 1 $NUM_ENS`; do

  # Run CRTM on this ensemble member
  $SCRIPT_DIR/job_submit.sh $crtm_ntasks $(($tid*$wrf_ntasks))  \
                            $HOSTPPN ./crtm.exe wrfinput_xf_$id \
                            met7_xf_$id.bin >& log.crtm_$id &

  # Check the number of active processes
  tid=$(($tid+1))

  # If needed, wait for all processes to clear
  if [[ $tid -le $ntid ]]; then
    wait
    tid=0
  fi

done

wait # Wait for all CRTM runs to clear


##############################################################
# SECTION 3: EVOKE MULTISMOOTH FAT CODE
##############################################################

for id in `seq -f "%03f" 1 %NUM_ENS`; do

  # Special exception for stampede
  if [[ $HOSTTYPE == 'stampede' ]]; then
    echo "python3 ens_msfat.py $fname_xf_crtm 


done
