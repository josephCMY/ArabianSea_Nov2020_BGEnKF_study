#!/bin/bash

# Script to run nccopy compressions on wrf files in the experiment

# Useful settings
source ~/.bashrc

# Useful functions
. util.sh



#-------------------------------------------------------------------
# SECTION 0: USER CONTROLS
#-------------------------------------------------------------------

# Date range
date_st=201705301200
date_ed=201706020000

# Date of first cycle
date_firstcycle=201705301200

# Number of processes to perform nccopy and compression
ntid=1
tid=0

# Number of wrf domains
ndom=1


# For each experiment, run CRTM
for expt in Thom_GTS+AMV+ch08_30min; do

  echo Running nccopy for $expt
  #load configuration files, functions, parameters
  export CONFIG_FILE=/global/homes/m/my_chan/sumatra_satellite_DA/PSU_EnKF/scripts/config/$expt

  . $CONFIG_FILE


  # Initiate date
  date=$date_st

  # Loop over dates in question
  while [[ $date -le $date_ed ]]; do

    echo "  "Processing $date

    mkdir -p $WORK_DIR/storage/$date

    cd $WORK_DIR/fc/$date

    ## Compressing ensemble files
    ## ----------------------------

    ## For each domain and each member
    #for dd in `seq -f "%02g" 1 $ndom`; do
    #  for ee in `seq -f "%03g" 1 $NUM_ENS`; do

    #    # Compress prior
    #    date2=`advance_time $date $CYCLE_PERIOD`
    #    f_in='wrfinput_d'$dd'_'`wrf_time_string $date2`_$ee
    #    f_out=$WORK_DIR/storage/$date/$f_in
    #    nccopy -d 1 $f_in $f_out &
    #    tid=$(($tid+1))
    #    if [[ $tid -ge $ntid ]]; then
    #      wait
    #      tid=0
    #    fi


    #    # Compress posterior
    #    f_in='wrfinput_d'$dd'_'$ee
    #    f_out=$WORK_DIR/storage/$date/$f_in
    #    nccopy -d 1 $f_in $f_out &
    #    tid=$(($tid+1))
    #    if [[ $tid -ge $ntid ]]; then
    #      wait
    #      tid=0
    #    fi

    # done
    #done


    ## Compress mean fields
    ## ---------------------
    #for dd in `seq -f "%02g" 1 $ndom`; do

    #  # Dealing with posterior
    #  f_in='wrfinput_d'$dd
    #  f_out=$WORK_DIR/storage/$date/$f_in
    #  nccopy -d 1 $f_in $f_out &
    #  tid=$(($tid+1))
    #  if [[ $tid -ge $ntid ]]; then
    #    wait
    #    tid=0
    #  fi

    #  # Dealing with prior
    #  f_in='wrf_enkf_input_d'$dd'_mean'
    #  f_out=$WORK_DIR/storage/$date/$f_in
    #  nccopy -d 1 $f_in $f_out &
    #  tid=$(($tid+1))
    #  if [[ $tid -ge $ntid ]]; then
    #    wait
    #    tid=0
    #  fi
    #done
  

    ## Copy over CRTM outputs
    ## ----------------------
    #cp him8_x*bin $WORK_DIR/storage/$date 


    # Compress deterministic forecasts
    # ---------------------------------
    for dd in `seq -f "%02g" 1 $ndom`; do
      for ff in `ls --color=never df_wrfout_d"$dd"*`; do
        # Dealing with forecast files
        f_in=$ff
        f_out=$WORK_DIR/storage/$date/$f_in
        nccopy -d 1 $f_in $f_out &
        tid=$(($tid+1))
        if [[ $tid -ge $ntid ]]; then
          wait
          tid=0
        fi
      done
    done
    cp him8_df*bin $WORK_DIR/storage/$date 


    # Prepare symbolic links for ensemble members
    # -------------------------------------------
    prevdate=`advance_time $date -$CYCLE_PERIOD`
    cd $WORK_DIR/storage/$date
    for dd in `seq -f "%02g" 1 $ndom`; do

      # Linking members
      for ee in `seq -f "%03g" 1 $NUM_ENS`; do

        # Dealing with posterior memebrs
        ln -s wrfinput_d"$dd"_$ee wrf_enkf_output_d"$dd"_$ee

        # Dealing with priors
        ln -s ../$prevdate/wrfinput_d"$dd"_`wrf_time_string $date`_$ee wrf_enkf_input_d"$dd"_$ee

      done

      # Linking means
      ln -s wrfinput_d"$dd" wrf_enkf_output_d"$dd"_mean 

    done 
    cd $WORK_DIR/fc/$date


    
    # Store the WRFOUT ensemble files
    # -------------------------------
    prevdate=`advance_time $date -$CYCLE_PERIOD`
    cd $WORK_DIR/storage/$date
    for dd in `seq -f "%02g" 1 $ndom`; do
      for ee in `seq -f "%03g" 1 $NUM_ENS`; do

        # Dealing with prior
        if [[ $date -gt 201705301200 ]]; then
          f_in=$WORK_DIR/run/$prevdate/wrf_ens/$ee/wrfout_d"$dd"_`wrf_time_string $date`
          f_out=wrfout_d"$dd"_xf_$ee
          nccopy -d 1 $f_in $f_out &
          tid=$(($tid+1))
          if [[ $tid -ge $ntid ]]; then
            wait
            tid=0
          fi
        fi

        # Dealing with posterior
        f_in=$WORK_DIR/run/$date/wrf_ens/$ee/wrfout_d"$dd"_`wrf_time_string $date`
        f_out=wrfout_d"$dd"_xa_$ee
        nccopy -d 1 $f_in $f_out &
        tid=$(($tid+1))
        if [[ $tid -ge $ntid ]]; then
          wait
          tid=0
        fi

      done
    done
    cd $WORK_DIR/fc/$date

    
    # Store the enkf logs
    # --------------------
    cd $WORK_DIR/storage/$date
    for dd in `seq -f "%02g" 1 $ndom`; do
      for ff in fort.10000 fort.10001 enkf.log; do
        cp $WORK_DIR/run/$date/enkf/d"$dd"/$ff $WORK_DIR/storage/$date
      done
    done
    cd $WORK_DIR/fc/$date


    # Increment time
    date=`advance_time $date $CYCLE_PERIOD`

  done # End of loop over time

done # End of loop over expt

# Wait for all background processes to clear before exiting
wait


