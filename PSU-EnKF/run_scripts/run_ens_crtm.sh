#!/bin/bash

#####header for stampede2######
#SBATCH -J ensCRTM
#SBATCH -p development
#SBATCH -N 10
#SBATCH -n 476
#SBATCH --time=02:00:00
#SBATCH -o log.ens_crtm
#SBATCH -e err.ens_crtm

# Parallelization control
tid=0
ntid=10


# Experiment list control
expt_list="NoDA"


# Date controls
date_st=201110160000
date_ed=201110160100
t_int=60

NUM_ENS_CRTM=5


# Defining advance time function
function advance_time {
  ccyymmdd=`echo $1 |cut -c1-8`
  hh=`echo $1 |cut -c9-10`
  mm=`echo $1 |cut -c11-12`
  inc=$2
  date -u -d $inc' minutes '$ccyymmdd' '$hh':'$mm +%Y%m%d%H%M
}



# Loop over experiments
for expt in $expt_list; do
 
  # Load relevant config file
  export CONFIG=/work2/04920/tg842199/stampede2/nonlinear_IR-DA/indian_ocean_osse/PSU_EnKF/scripts/config/$expt
  . $CONFIG



  # Iterative date loop
  date=$date_st
  while [[ $date -le $date_ed ]]; do

    echo Dealing with $expt $date

    # Entering date directory
    cd $WORK_DIR/fc/$date

    # Setting up crtm symbolic link
    rm crtm.exe
    ln -sf $CODE_DIR/CRTM/crtm_wrf/met7_iodc/crtm.exe .
    ln -sf $CODE_DIR/CRTM/crtm_wrf/met7_iodc/crtm.cloudless.exe .
    ln -sf $CODE_DIR/CRTM/crtm_wrf/met7_iodc/coefficients .


    # Running CRTM on wrfinput_d01
    # If DA was performed, wrfinput_d01 is the posterior
    # If no DA, wrfinput_d01 is the prior.
    echo "    Running CRTM on posterior ensemble"
    for ee in `seq -f "%03g" 1 $NUM_ENS_CRTM`; do
      if [[ ! -f him8_xa_"$ee".bin ]]; then

        # Perform cloudy CRTM
        ibrun -n 68 -o $((tid*68)) crtm.exe wrfinput_d01_$ee met7_xa_"$ee".bin >& log.crtm_xa_$ee &
        tid=$(($tid+1))
        if [[ $tid -ge $ntid ]]; then
          tid=0
          wait
        fi

        # Perform cloudless CRTM
        ibrun -n 68 -o $((tid*68)) crtm.cloudless.exe wrfinput_d01_$ee met7_cloudless_xa_"$ee".bin >& log.crtm_cloudless_xa_$ee &
        tid=$(($tid+1))
        if [[ $tid -ge $ntid ]]; then
          tid=0
          wait
        fi


      fi
    done
  

    ## Running CRTM on wrf_enkf_input_d01 (if they exist)
    ## wrf_enkf_input_d01 only exists in the situation where DA was performed
    #echo "    Running CRTM for prior ensemble"
    #for ee in `seq -f "%03g" 1 $NUM_ENS`; do

    #  # Situation where DA happened (wrf_enkf_input exists)
    #  if [[ -e wrf_enkf_input_d01_$ee ]]; then
    #    if [[ ! -f him8_xf_"$ee".bin ]]; then
    #      srun -N 1 -n 68 crtm.exe wrf_enkf_input_d01_$ee him8_xf_"$ee".bin >& log.crtm_xf_$ee &
    #      tid=$(($tid+1))
    #      if [[ $tid -ge $ntid ]]; then
    #        tid=0
    #        wait
    #      fi
    #    fi

    #  # Situation where DA did not happen
    #  else
    #    ln -sf him8_xa_"$ee".bin him8_xf_"$ee".bin
    #  fi

    #done # End of iteration over ens


    # Advancing time
    date=`advance_time $date $t_int`

  done # End of date loop
done # End of expt loop

wait
