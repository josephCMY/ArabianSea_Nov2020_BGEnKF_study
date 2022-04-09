#========================================================================
# Script to generate appropriate wrfbdy for rapid cycling experiments.
# -----------------------------------------------------------------------
# Author: Man-Yau Chan
#
# Description:
# ------------
# This script will read in the original wrfbdy file, wrf date string, the 
# cycling interval (minutes), and the output wrfbdy file name.
# Then, the script will construct a wrfbdy file containing the boundary 
# values and tendencies for the time in question.
#
# This script uses linear interpolation to interpolate boundary values and
# tendencies between the times available in the original wrfbdy file.
# 
# Need to carefully redesign the netCDF4 attributes and all that.
# ======================================================================


# Basic libraries
import netCDF4 as nc4
import numpy as np
import math as m
import copy
import sys
import datetime
import warnings

warnings.filterwarnings("ignore")

# Read in desired arguments
# -------------------------
fname_src = sys.argv[1]
date_out = datetime.datetime.strptime( sys.argv[2], '%Y%m%d%H%M' )
time_interval = int( sys.argv[3] )
fname_out = sys.argv[4]


# Open netcdf files for reading and writing
# ------------------------------------------
f_src = nc4.Dataset(fname_src, 'r')
f_out = nc4.Dataset(fname_out, 'w')



# ==================================================================================
# SECTION 1: SET UP TIME INTERPOLATION
# ==================================================================================

# Load the time variable and convert to minutes since first date on bdy file
# ---------------------------------------------------------------------------
src_times = nc4.chartostring( f_src.variables['Times'][:]  )
src_times = [ datetime.datetime.strptime( tt, '%Y-%m-%d_%H:%M:%S') for tt in src_times ]
time_coord_list = [ ( tt - src_times[0] ).total_seconds()/60 for tt in src_times ]

# Determine time index of the f_src date that immediately preceeds date_out
# --------------------------------------------------------------------------
time_coord_out0 = ( date_out - src_times[0]).total_seconds()/60
t_ind_left0 = np.searchsorted( time_coord_list, time_coord_out0, side='left' )
if ( time_coord_out0 < time_coord_list[t_ind_left0] ):
  t_ind_left0 -= 1

# Determine time index of f_src date that immediately preceeds date_out + time_interval
# -------------------------------------------------------------------------------------
time_coord_out1 = ( date_out + datetime.timedelta( minutes=time_interval ) 
                    - src_times[0]
                  ).total_seconds()/60
t_ind_left1 = np.searchsorted( time_coord_list, time_coord_out1, side='left' )
if ( time_coord_out1 < time_coord_list[t_ind_left1] ):
  t_ind_left1 -= 1


# Determine time interpolation weights (needed for linear time interpolation later)
# ----------------------------------------------------------------------------------
# For date_out
t0 = time_coord_list[ t_ind_left0 ]
t1 = time_coord_list[ t_ind_left0 + 1 ]
t_weights0 = [ (t1 - time_coord_out0) / (t1 - t0), \
               (time_coord_out0 - t0) / (t1 - t0) ]
# For date_out + time_interval
t0 = time_coord_list[ t_ind_left1 ]
t1 = time_coord_list[ t_ind_left1 + 1 ]
t_weights1 = [ (t1 - time_coord_out1) / (t1 - t0), \
               (time_coord_out1 - t0) / (t1 - t0) ]




# =================================================================================
# SECTION 2: SET UP NON-VARIABLE QUANTITIES IN OUTPUT FILE
# =================================================================================

# Initialize output file dimensions and attributes
# ------------------------------------------------
# Iterate thru each source file attribute, and copy
# over to the output file.
for attr_name in f_src.ncattrs():
  attr = f_src.getncattr( attr_name )
  f_out.setncattr( attr_name, attr )


# Iterate through each source file dimension, and
# copy dimensions over to the output file
for dimname in f_src.dimensions.keys():
  # Load dimension size 
  dim_size = f_src.dimensions[dimname].size

  # Special handling for time dimension since the output 
  # file will only contain bdy values and tendencies at
  # 2 time points (start and end of fcst cycle)
  if dimname == 'Time':
    dim_size=2
  # End of special time dimension handling

  # Construct dimension in output file
  f_out.createDimension( dimname, dim_size )




# Construct output file time variables (there are three of them)
# --------------------------------------------------------------
# Construct time strings corresponding to start and end times of forecast cycle
time0 = date_out.strftime('%Y-%m-%d_%H:%M:%S')
time1 = ( date_out + datetime.timedelta( minutes=time_interval ) ).strftime('%Y-%m-%d_%H:%M:%S')

# Special time string for meta data
time2 = ( date_out + datetime.timedelta( minutes=time_interval*2 ) ).strftime('%Y-%m-%d_%H:%M:%S')

# Converting Python strings to character byte arrays for netCDF3 files
time0 = nc4.stringtochar( np.array( time0, dtype=("S%d" % len(time0)) ) )
time1 = nc4.stringtochar( np.array( time1, dtype=("S%d" % len(time1)) ) )
time2 = nc4.stringtochar( np.array( time2, dtype=("S%d" % len(time2)) ) )

# Initialize Times variable and store stuff
var = f_out.createVariable('Times','S1', dimensions=('Time','DateStrLen',) )
var[0] = time0
var[1] = time1

# Initialize meta data variable for current time
var = f_out.createVariable( 'md___thisbdytimee_x_t_d_o_m_a_i_n_m_e_t_a_data_', 'S1',
                             dimensions=('Time','DateStrLen',) )
var[0] = time0
var[1] = time1

# Initialize meta data variable for the next time
var = f_out.createVariable('md___nextbdytimee_x_t_d_o_m_a_i_n_m_e_t_a_data_', 'S1',
                            dimensions=('Time','DateStrLen',) )
var[0] = time1
var[1] = time2



# ===================================================================================
# SECTION 3: TIME-INTERPOLATE BOUNDARY CONDITIONS TO DESIRED TIMES
# ===================================================================================

# Iterate through all float variables
for var_key in f_src.variables.keys():

  # Load variable
  src_var = f_src.variables[var_key]

  # Skip over non-numeric variables
  if ( not np.issubdtype( src_var.dtype, np.number ) ):
    continue

  # Now interpolate to first time point
  out_var0 = src_var[ t_ind_left0 ] * t_weights0[0] \
             + src_var[ t_ind_left0 + 1 ] * t_weights0[1]

  # Now interpolate to second time point
  out_var1 = src_var[ t_ind_left1 ] * t_weights1[0] \
             + src_var[ t_ind_left1 + 1 ] * t_weights1[1]


  # Store output into f_out
  out_var = f_out.createVariable( var_key, src_var.dtype, dimensions=src_var.dimensions )
  out_var[0] = out_var0
  out_var[1] = out_var1

  # Add variable attributes
  for attr_name in src_var.ncattrs():
    attr =  src_var.getncattr( attr_name )
    out_var.setncattr( attr_name, attr )


# Close the stuff to save
f_out.close()
f_src.close()

