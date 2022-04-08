#!/bin/bash
#SBATCH --job-name="setup_ens"
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=127
#SBATCH --cpus-per-task=1
#SBATCH --export=ALL
#SBATCH -t 02:00:00
#SBATCH --account=pen116
#SBATCH --output="setup_ens.log"

# ======================================================================== #
# Script to apply the WRF Pre-processing System onto the GEFS GRIB files   #
# Author: Man-Yau (Joseph) Chan                                            #
#									   #
# IMPORTANT NOTE: total number of processes must be > ENSEMBLE_SIZE        #
# ======================================================================== #


# Loading configuration and utility functions 
# -------------------------------------------
. config_ens_setup
. util.sh


# Make directories to process stuff
# ---------------------------------
mkdir -p $RAW_GEFS_DIR
mkdir -p $WPS_GEFS_DIR


echo ' '
echo ========================================================
echo Starting to construct WRF ensemble from GEFS
echo `date`
echo --------------------------------------------------------


# Download data from AWS
# ----------------------
echo ' '
echo "Downloading data from AWS"
./get_GEFS.sh >& log.download_GEFS 
echo "Finished downloading data from AWS"


# Construct geographical domain
# -----------------------------
echo ' '
echo Executing geogrid.exe
# Enter the designated directory
mkdir -p $WPS_GEFS_DIR/setup_geog_domain
cd $WPS_GEFS_DIR/setup_geog_domain

# Symbolic link in stuff needed for geogrid
ln -sf $WPS_NAMELIST .
ln -sf $WPS_DIR/geogrid/geogrid.exe .
ln -sf $WPS_DIR/geogrid/GEOGRID.TBL .

# Run geogrid
module restore intel
ibrun geogrid.exe >& geogrid.log 
module restore default

echo ' '
echo "Geogrid finished running"
echo `date`
echo --------------------------------------------------------


# Convert every member's GRIB file into WPS intermediate files via ungrib
# -----------------------------------------------------------------------
echo ' '
echo "Executing ungrib.exe"
for ee in `seq -f "%03g" 1 $ENSEMBLE_SIZE`; do
  
  # Construct directory and enter
  mkdir -p $WPS_GEFS_DIR/ungrib_member_$ee
  cd $WPS_GEFS_DIR/ungrib_member_$ee 

  # Symlink in stuff from WPS' ungrib directory
  ln -s $WPS_DIR/ungrib.exe .
  ln -s $WPS_DIR/ungrib/Variable_Tables/Vtable.GFSENS Vtable
  ln -s $WPS_DIR/link_grib.csh .
  ln -s $WPS_GEFS_DIR/setup_geog_domain/geo_em*nc .
  ln -s $WPS_NAMELIST .

  # Construct links to the appropriate grib files
  ./link_grib.csh $RAW_GEFS_DIR/$ee/gep*

  # Call ungrib on GRIB files and move on
  module restore intel
  ibrun -n 1 ungrib.exe >& ungrib.log &
  module restore default
  
done

# Wait for all ungrib executions to finish
wait

echo ' '
echo "Ungrib finished running"
echo `date`
echo --------------------------------------------------------



# Construct met_em files from every member's WPS intermediate files
# ------------------------------------------------------------------
echo ' '
echo Executing metgrid on the GEFS ensemble
for ee in `seq -f "%03g" 1 $ENSEMBLE_SIZE`; do

  # Construct directory and enter
  mkdir -p $WPS_GEFS_DIR/metgrid_member_$ee
  cd $WPS_GEFS_DIR/metgrid_member_$ee

  # Symlink in the metgrid stuff
  ln -s $WPS_DIR/metgrid/metgrid.exe .
  ln -s $WPS_DIR/metgrid/METGRID.TBL.ARW METGRID.TBL
  ln -s $WPS_DIR/metgrid/gribmap.txt .
  ln -s $WPS_NAMELIST namelist.wps
  ln -s $WPS_GEFS_DIR/setup_geog_domain/geo_em*nc .
  ln -s $WPS_GEFS_DIR/ungrib_member_$ee/* .


  # Run metgrid and move on to next member
  module restore intel
  ibrun -n 1 metgrid.exe >& metgrid.log &
  module restore default

done

# Wait for all metgrid executions to finish
wait
echo ' '
echo "Metgrid finished running"
echo `date`
echo --------------------------------------------------------


# Construct convert met_em files into wrfinput and wrfbdy files
# --------------------------------------------------------------
echo ' '
echo 'Converting met_em files to wrfinput and wrfbdy files'
for ee in `seq -f "%03g" 1 $ENSEMBLE_SIZE`; do

  # Construct directory and enter
  mkdir -p $WPS_GEFS_DIR/real_member_$ee
  cd $WPS_GEFS_DIR/real_member_$ee

  # Symlink in things needed to use real.exe
  ln -s $WRF_DIR/main/real.exe 
  ln -s $WRF_REAL_NAMELIST namelist.input
  ln -s $WPS_GEFS_DIR/metgrid_member_$ee/met_em*.nc .

  # Run real.exe and move onto the next member
  module restore intel
  ibrun real.exe >& real.log 
  module restore default
  echo Finished running real.exe on $ee

done

echo " "
echo "WRF real.exe finished running"
echo `date`
#echo --------------------------------------------------------



# Executing ensemble spin-up
# --------------------------
echo ' '
echo Executing ensemble spin-up
for ee in `seq -f "%03g" 1 $ENSEMBLE_SIZE`; do

  # Construct directory and enter
  mkdir -p $WPS_GEFS_DIR/spinup_member_$ee
  cd $WPS_GEFS_DIR/spinup_member_$ee

  # Symlink in things needed to use real.exe
  ln -s $WRF_DIR/main/wrf.exe 
  ln -s $WRF_SPINUP_NAMELIST namelist.input
  ln -s $WPS_GEFS_DIR/real_member_$ee/wrfinput_d01 .
  ln -s $WPS_GEFS_DIR/real_member_$ee/wrfbdy_d01 .
  ln -s $WPS_GEFS_DIR/real_member_$ee/wrflowinp_d01 .
  ln -s $WRF_DIR/run/* . 


  # Run real.exe and move onto the next member
  module restore intel
  ibrun wrf.exe >& wrf.log 
  module restore default

  echo Finished running spinup on member $ee

done


echo " "
echo "Finished spinnup up ensemble"
echo `date`
echo --------------------------------------------------------



echo ' '
echo ========================================================
echo Finished constructing the WRF ensemble from the GEFS
echo `date`
echo ========================================================
echo ' '
