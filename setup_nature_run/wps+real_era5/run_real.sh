#!/bin/bash

# Script to run real.exe to turn met_em and geo_em files to 
# wrfbdy_d01, wrflowinp_d01 and wrfinput_d01

# Load config file
. ../config

# Make namelist
./setup_namelist_real.sh > processed_era5/namelist.input

# Enter into scratch directory
cd ~/DOE_reanalysis/setup_ens/wps+real_era5/processed_era5

# Symbolic link in real.exe
ln -sf ~/DOE_reanalysis/PSU_EnKF/code/WRFV3/run/real.exe

# Evoke real.exe 
echo " "
echo "Running real.exe"
#salloc -C knl -N 4 -q interactive -t 04:00:00 --job-name='real' srun -N 4 -n 272 real.exe >& real.log
ibrun -n

mv rsl.out.0000 rsl.log.$date_st"_"$date_ed
rm rsl.out.* rsl.error*

# Rename files
for fname in wrfbdy_d01 wrfinput_d01 wrflowinp_d01; do
  mv $fname "$fname"_"$date_st"_"$date_ed"
done

# Return to HOME directory
cd ~/DOE_reanalysis/setup_ens/wps+real_era5

