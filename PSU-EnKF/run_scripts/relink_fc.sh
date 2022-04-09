#!/bin/bash

# Recently reorganised the EnKF experiments' directories. This script is to perform relinking in fc

. util.sh

tint=60

dir0=$SCRATCH/nonlinear_IR_DA/indian_ocean_osse/PSU_EnKF
config_dir=$SCRATCH/nonlinear_IR_DA/indian_ocean_osse/PSU_EnKF/NoDA/DA/config
cd $dir0

# Iterate over available expts
for expt in NoDA; do

  echo Patching up symbolic links in $expt

  # Calling on the relevant config file
  . $config_dir/$expt

  # Redoing the links in the fc directory
  cd $WORK_DIR/fc
  for date in `ls`; do
    echo $date
    cd $WORK_DIR/fc/$date

    prev_date=`advance_time $date -$tint`
  
    # For each ens member
    for mem in `seq 1 $NUM_ENS`; do

      eee=`printf "%03g" $mem`
      
      # If NoDA, wrfinput_d01 is from the previous date
      if [[ $expt == 'NoDA' ]] ; then

        # Relink wrfinput_d01 files
        prev_date_file=$WORK_DIR/run/$prev_date/wrf_ens/$eee/wrfinput_d01_`wrf_time_string $date`
        ln -sf $prev_date_file wrfinput_d01_`wrf_time_string $date`"_"$eee
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


