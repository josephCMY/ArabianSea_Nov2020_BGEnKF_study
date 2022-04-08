#!/bin/bash

# Script to extract 10.8 and 6.2 micron BTs within study domain from SEVIRI native data files.

# Important notes:
# -----------------
# 1) study domain geog limits are specified within convert_seviri_nat2nc.py
# 2) the SEVIRI native files must live within seviri_native_file directory
# 3) the outputted NetCDF4 files will be saved within extracted_seviri_ncfiles
# 4) all available native files will be converted
# 5) parallelization will be employed to convert multiple native files at a time
# 6) to turn off parallelization, just set max_active_procs=1


# Parallelization controls
max_active_procs=3


# Count number of files to process
tot_cnt=0
for ff in `ls --color=none seviri_native_file/*.nat`; do
  tot_cnt=$(( $tot_cnt + 1 ))
done
echo Total number of native files to extract: $tot_cnt


# Perform data extraction
num_active_procs=0
cnt=0
for ff in `ls --color=none seviri_native_file/*.nat`; do

  # Initiate data extraction for the desired file
  python convert_seviri_nat2nc.py $ff >> log.proc_$num_active_procs &
  cnt=$(( $cnt + 1 ))
  echo Processing $cnt out of $tot_cnt files
  num_active_procs=$(( $num_active_procs + 1 ))

  # Wait for all procs to clear if num_active_procs >= max_active_procs
  if [[ $num_active_procs -ge $max_active_procs ]]; then
    wait
    num_active_procs=0
  fi

done
wait

echo Finished processing all $tot_cnt files
