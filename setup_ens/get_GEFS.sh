#!/bin/bash
# Script to download the GEFS files

# Read config file 
. config_get_GEFS

# Load utility functions
. util.sh

 
# Echo date and some msg
date
echo "Starting to download GEFS data"


# ===============================================================
# PART A: For each PSU-EnKF ensemble member, construct a list of 
# GEFS GRIB filenames on the AWS server
# ===============================================================
#
# General idea: 
# -------------
# Suppose the GEFS has 30 members and initiates every 6 hrs, 
# and we want 91 members.
#
# 1) First 30 members use analysis states at 00, 06, 12, 18 UTC
# 2) Next 30 members use 6-fcst states valid at 00, 06, 12, 18 UTC
# 3) Next 30 members use 12-fcst states valid at 00, 06, 12, 18 UTC
# 4) Final member uses 18-fcst state 
#
# -------------------------------------------------------------

# Iterating thru each PSU-EnKF member
for id in `seq 1 $ENSEMBLE_SIZE`; do

  # Initialize variable to store all AWS paths for this member
  fname_list=""

  # Determine which GEFS lead time this member uses
  lead_time=$(( ( ($id-1)/$GEFS_SIZE )*$GEFS_INIT_INTERVAL ))
  lead_time_str=`printf "%03g" $lead_time`

  # Pretty print ensemble ID
  id_str=`printf "%03d" $id`
 
  # Determine the GEFS ensemble id corresponding to this member
  gefs_id=$(( $id - ( ($id-1)/$GEFS_SIZE )*$GEFS_SIZE  ))
  echo PSU-EnKF ID $id
  echo GEFS ID $gefs_id
  
  # Pretty print GEFS ensemble id
  gefs_id_str=`printf "%02d" $gefs_id`

  # Construct all the file names desired
  date_nw=$DATE_START
  while [[ $date_nw -le $DATE_END ]]; do

    # Compute init date corresponding to date_nw and lead_hours
    date_init=`advance_time $date_nw -$(( $lead_time*60 ))`

    # Extract relevant elements of initiation date
    ccyymmdd=`echo $date_init |cut -c1-8`
    HH=`echo $date_init |cut -c9-10`

    # Iterating thru path types
    for raw_pgrb_fmt in $PGRB2A_PATH_FORMAT $PGRB2B_PATH_FORMAT; do
      # Load path format specified in config file config file
      aws_path=$raw_pgrb_fmt

      # Overwrite date and time-related items in file path
      aws_path=`echo $aws_path | sed 's|ccyymmdd|'$ccyymmdd'|g' `
      aws_path=`echo $aws_path | sed 's|HH|'$HH'|g'`
      aws_path=`echo $aws_path | sed 's|LL|'$lead_time_str'|g'`
      aws_path=`echo $aws_path | sed 's|ID|'$gefs_id_str'|g'`

      # Store the AWS path into the fname_list
      fname_list="$fname_list"$'\n'"$aws_path"
    
    done # --- Iterated over both pgrb2 types 

    # Advance date
    date_nw=`advance_time $date_nw $(( $GEFS_INIT_INTERVAL * 60 ))`

  done # --- Iterated through all dates

  # Saving the list of filenames produced
  mkdir -p $RAW_GEFS_DIR/$id_str
  echo "$fname_list" >& $RAW_GEFS_DIR/$id_str/aws_path_list.txt

done # --- Iterated over all members

  

# ============================================================================
# PART B: Download GEFS GRIB files from the AWS
# ----------------------------------------------------------------------------
# In part A, we constructed text files containing lists of AWS paths that
# point to the desired GEFS grib files. 
# In this part, we will use rclone to download these files.
# ----------------------------------------------------------------------------
active_procs=0
# Iterate over all ensemble members
for id_str in `seq -f "%03g" 1 $ENSEMBLE_SIZE`; do

  # Enter the appropriate directory
  cd $RAW_GEFS_DIR/$id_str

  # Load the list of file names
  fname_list=`cat aws_path_list.txt`


  # Download GRIB files from AWS server
  for ff in $fname_list; do
    echo Rcloning $ff
    ~/rclone-v1.53.2-linux-amd64/rclone copy AWS_us-east-1:$ff . 
  done


done

date
echo "Finished downloading GEFS data"
