# Script to generate appropriate wrflowinp for rapid cycling experiments.
# From eyeballing the technical note for WRFV3, it appears that the bottom boundary is
# read in at the netCDF frequency. I.e, no temporal interpolation for the btm bdy.

# This script will read in the original wrflowinp file, wrf date string and output 
# wrflowinp filename. It will then seek the entry in the wrflowinp files that immediately
# precedes the desired time.

# Creates a wrflowinp with values only at the specified time

# Need to carefully redesign the netCDF4 attributes and all that.

# Basic libraries
import netCDF4 as nc4
import numpy as np
import math as m
import copy
import sys
import datetime

# Read in desired arguments
fname_src = sys.argv[1]
time_out = datetime.datetime.strptime( sys.argv[2], '%Y%m%d%H%M' )
#print( time_out)
fname_out = sys.argv[3]



# Functions to handle date string manipulations
# ---------------------------------------------
# Read in wrf date string
def wrf_time_str_read( t_str ):
  out = datetime.datetime.strptime( t_str, '%Y-%m-%d_%H:%M:%S')
  return out

# Increment time in terms of seconds
def wrf_time_change( now, dt ):
  now += datetime.timedelta( seconds = dt )
  return now

# Return wrf datestring
def wrf_time_str( now ):
  return now.strftime('%Y-%m-%d_%H:%M:%S')


# Interpolate wrflowinp to desired date
# -------------------------------------

# Open netcdf files for reading and writing
f_src = nc4.Dataset(fname_src, 'r')
f_out = nc4.Dataset(fname_out, 'w')


# Load and compute array dimensions
dims = f_src.dimensions
dims_key = ['Time','DateStrLen','west_east','south_north']
dims_len = [ dims[i].size for i in dims_key]
dims_len[0] = 1
dims_input = [None, dims_len[1], dims_len[2], dims_len[3] ]
var_shp = [ dims_len[0], dims_len[2], dims_len[3] ]

# Insert array dimensions
[times, datestr, west_east, south_north] \
= [ f_out.createDimension(dims_key[i], dims_input[i]) for i in range(4) ]

# Read in all global attributes and insert
g_attr = f_src.ncattrs()
for attr_name in g_attr:
  attr = f_src.getncattr( attr_name )
  f_out.setncattr( attr_name, attr)

# Read in all the times in the source file, and then iteratively seek out the time index
# immediately preceding the desired time
src_times = f_src.variables['Times'][:]
#print((np.array(src_times[0,:]).tostring()).decode('UTF-8'))
src_times = [ wrf_time_str_read( (np.array(src_times[i]).tostring()).decode('UTF-8') ) for i in range( len(src_times) ) ]
for t_ind in range( len(src_times) ):
  if time_out < src_times[t_ind]:
    break
# preceding index is the one right before breaking
t_ind -= 1
#print( t_ind)

# Generating the desired wrflowinp
var_keys = f_src.variables.keys()
#print( ['desired time' src_times[t_ind] )

for var_key in var_keys:
  # Read in variable
  var = f_src.variables[var_key]
  var_arr = (var)[:]
  var_shp = var_arr.shape

  # Special handling for Times variable
  if var_key == 'Times':
    # Generate the time value desired
    new_arr = [ list(wrf_time_str(time_out)) ]
    # Store new variable
    new_var = f_out.createVariable('Times','S1',dimensions=('Time','DateStrLen',))
    new_var[0] = new_arr[:]

  # Other variables
  else:
    # Set up new array
    new_arr = np.zeros( [1, var_shp[1], var_shp[2] ])
    new_arr[0,:,:] = var_arr[ t_ind, :,: ]
    new_var = f_out.createVariable( var_key, 'f4', dimensions = ('Time','south_north','west_east',) )
    # Insert array data
    new_var[0] = new_arr[:,:,:]
    # Include other attributes
    attrs_list = var.ncattrs()
    for attr_name in attrs_list:
      attr = var.getncattr( attr_name )
      new_var.setncattr( attr_name, attr )


# Close the stuff to save
f_out.close()
f_src.close()
