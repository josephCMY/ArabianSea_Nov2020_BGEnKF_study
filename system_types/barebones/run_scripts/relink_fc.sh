#!/bin/bash

# Recently reorganised the EnKF experiments' directories. This script is to perform relinking in fc

. util.sh

tint=180

dir0=$SCRATCH/sumatra_satellite_DA_ens/enkf_expts
config_dir=~/sumatra_satellite_DA/PSU_EnKF/scripts/config
cd $dir0

# Iterate over available expts
for expt in aeThom_NoDA aeThom_GTS WDM6_NoDA WDM6_GTS WDM6_GTS+Him8_LM WDM6_GTS+Him8_NLM; do

  echo Patching up symbolic links in $expt

  # Calling on the relevant config file
  . $config_dir/$expt

  # Redoing the links in the fc directory
  cd $WORK_DIR/fc
  for date in `ls`; do
    echo $date
    cd $WORK_DIR/fc/$date

    # Determinin the previous date
    if [[ $date == 201705310000 ]]; then
      prev_date=201705300000
    else
      prev_date=`advance_time $date -$tint`
    fi

    # Ignore spinup
    if [[ $date == 201705300000 ]]; then
      continue
    fi


  
    # For each ens member
    for mem in `seq 1 $NUM_ENS`; do

      eee=`printf "%03g" $mem`
      
      # If NoDA, wrfinput_d01 is from the previous date
      if [[ $expt == 'WDM6_NoDA' ]] || [[ $expt == 'aeThom_NoDA' ]] || [[ $expt == 'Thom_NoDA' ]]; then

        # Relink wrfinput_d01 files
        ln -sf ../$prev_date/wrfinput_d01_`wrf_time_string $date`_$eee wrfinput_d01_$eee

      # If not NoDA, we need to link a load of stuff
      else
        # Link the wrf enkf files
        ln -sf ../$prev_date/wrfinput_d01_`wrf_time_string $date`_$eee wrf_enkf_input_d01_$eee
        nn=$((90010+$mem))
        ln -sf $WORK_DIR/run/$date/enkf/d01/fort.$nn wrf_enkf_output_d01_$eee
        # Link the wrfinput file
        ln -sf wrf_enkf_output_d01_$eee wrfinput_d01_$eee
      fi
    done 
  done

done


