#!/bin/bash
# Script to compress and archive a specified run directory into HPSS

. util.sh

maindir=/scratch/04920/tg842199/nonlinear_IR_DA/indian_ocean_osse/PSU_EnKF

date_st=201110151200
date_ed=201110180000

# Dealing with run directory files
for expt in JointSpace_EnKF_IR_only JointSpace_EnKF_WV_only NoDA  StateSpace_EnKF_IR_only  StateSpace_EnKF_WV_only; do

  homedir=$maindir/$expt
  rdir=run

  # For each specified date in run dir
  date=$date_st
  while [[ $date -le $date_ed ]]; do    #--------------- date loop start

    ##########################################
    # Compressing the wrf_ens files
    ##########################################
    for ee in `seq -f "%03g" 1 50`; do  #--------------- WRF ens loop start
      cd $homedir/$rdir/$date/wrf_ens/$ee

      # Dealing with lower boundary conditions
      if [ "$ee" == "001" ] ; then
        deflate_ncfile wrflowinp_d01  &
      fi

      # Dealing with wrfinput_d01 and wrfout files
      for fname in `ls --color=none wrfinput_d01_* wrfout_d01*`; do
        echo $fname
        deflate_ncfile $fname &
      done

      # Wait for this member to be comrpessed
      wait
      echo Compressed mem $ee at `date`

    done                                #--------------- WRF ens loop end

    #########################################
    # Compressing enkf/d01 mean and posterior files
    #########################################
    cd $homedir/$rdir/$date/enkf/d01
    for ee in `seq 1 51`; do
      # Posterior file
      deflate_ncfile fort.$((90010+$ee))
    done
    deflate_ncfile fort.90061 
    echo Deflated $homedir/$rdir/$date/enkf/d01 files



    #########################################
    # Advance time for 1 hour
    #########################################
    date=`advance_time $date 60`

  done                                  #--------------- date loop end

done


