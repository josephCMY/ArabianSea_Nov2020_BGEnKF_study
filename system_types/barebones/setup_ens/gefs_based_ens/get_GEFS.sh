#!/bin/bash
# Script to download the GEFS files

# Read config file 
. config_ens_setup

# Load utility functions
. util.sh

 
# Echo date and some msg
date
echo "Starting to download GEFS data"


# ===============================================================
# PART A: For each PSU-EnKF ensemble member, construct a list of 
# GEFS GRIB filenames on the AWS server
# ===============================================================
# General idea: use multiple sets of GEFS forecasts to form the
# desired ensemble. I.e., use members from time-staggered forecasts
# 
# E.g., to generate a 60-member ensemble:
# The first 20 members will use initial and boundary conditions
# based on the GEFS forecasts initiated at time $DATE_START
# The second 20 members will use initial and boundary conditions
# based on the GEFS forecasts initiated at $DATE_START - 6 hours
# The third 20 members will use initial and boundary conditions
# based on the GEFS forecasts initiated at $DATE_START - 12 hours
# These three sets of 20-member forecasts are all valid at the 
# same date.
# -------------------------------------------------------------

# Iterating thru each PSU-EnKF ensemble member
for id in `seq 1 $ENSEMBLE_SIZE`; do
  
  # Initialize variable to hold GRIB file names
  fname_list=''

  # Figure out which GEFS initiation date is used by this ensemble member
  lag_hours=$(( ( ($id-1)/$GEFS_SIZE )*$GEFS_TIME_INTERVAL ))
  date_init=`advance_time $DATE_START -$(( $lag_hours*60 ))`

  # Pretty print ensemble ID
  id_str=`printf "%03d" $id`
 
  # Determine the GEFS ensemble id corresponding to this member
  gefs_id=$(( $id - ( ($id-1)/$GEFS_SIZE )*$GEFS_SIZE  ))
  echo PSU-EnKF ID $id
  echo GEFS ID $gefs_id

  # Pretty print GEFS ensemble id
  gefs_id_str=`printf "%02d" $gefs_id`

  # Now figure out all the forecast lead times needed
  lead_time=$lag_hours
  lead_time_list=$lead_time
  date_nw=$date_init
  while [[ `advance_time $date_nw $(( $lead_time*60 ))` -le $DATE_END ]]; do
   
    # Increment lead_time
    lead_time=$(( $lead_time + $GEFS_TIME_INTERVAL ))

    # Store new lead time
    if [[ `advance_time $date_nw $(( $lead_time*60 ))` -le $DATE_END ]]; then
      lead_time_list=$lead_time_list" "$lead_time
    fi
    
  done
  echo Lead times generated $lead_time_list


  # Extract relevant elements of initiation date
  ccyymmdd=`echo $date_init |cut -c1-8`
  HH=`echo $date_init |cut -c9-10`


  # Construct list of AWS paths to appropriate pgrb files and store in fname_list
  for raw_pgrb_fmt in $PGRB2A_PATH_FORMAT $PGRB2B_PATH_FORMAT; do
    for lead_time in $lead_time_list; do

      # Make pretty string
      lead_time_str=`printf "%03g" $lead_time`

      # Load path format specified in config file config file
      aws_path=$raw_pgrb_fmt

      # Overwrite date and time-related items in file path
      aws_path=`echo $aws_path | sed 's|ccyymmdd|'$ccyymmdd'|g' `
      aws_path=`echo $aws_path | sed 's|HH|'$HH'|g'`
      aws_path=`echo $aws_path | sed 's|LL|'$lead_time_str'|g'`
      aws_path=`echo $aws_path | sed 's|ID|'$gefs_id_str'|g'`

      # Store the AWS path into the fname_list
      fname_list="$fname_list"$'\n'"$aws_path"
    
    done # --- Iterated over lead times
  done # --- Iterated over both types of pgrb files


  # Store file name list into a textfile in RAW_GEFS_DIR
  mkdir -p $RAW_GEFS_DIR/$id_str
  echo "$fname_list" >& $RAW_GEFS_DIR/$id_str/aws_path_list.txt


done # --- Iterated over ensemble members



# ============================================================================
# PART B: Download GEFS GRIB files from the AWS
# ----------------------------------------------------------------------------
# In part A, we constructed text files containing lists of AWS paths that
# point to the desired GEFS grib files. 
# In this part, we will use rclone to download these files.
# ----------------------------------------------------------------------------

# Iterate over all ensemble members
for id_str in `seq -f "%03g" 1 $ENSEMBLE_SIZE`; do

  # Enter the appropriate directory
  cd $RAW_GEFS_DIR/$id_str

  # Load the list of file names
  fname_list=`cat aws_path_list.txt`


  # Download GRIB files from AWS server
  for ff in $fname_list; do
    echo Rcloning $ff
    rclone copy AWS:$ff . 
  done

done

date
echo "Finished downloading GEFS data"
