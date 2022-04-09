'''
Script to check if two netcdf files are identical

Takes names of the two netcdf files as inputs

Will check all variables that are not strings
'''

# Standard libraries
import numpy as np
from netCDF4 import Dataset as ncopen
import sys

# Names of netcdf files being compared
ncfnames = [ sys.argv[1], sys.argv[2] ]

# Open the files
ncfiles = [ ncopen( ncfnames[i], 'r') for i in [0,1] ]

# Quit flag (triggers whenever something is mismatched)
quit_flag = False


''' Check if both systems have identical variable names '''
keylists = [ ncfiles[i].variables.keys() for i in [0,1] ]

# Performing two-way check
for key in keylists[0]:
    if key in keylists[1]:
        continue
    else:
        print( "Variable %s present in %s, absent in %s" % (key, ncfnames[0], ncfnames[1]))
        quit_flag = True
for key in keylists[1]:
    if key in keylists[0]:
        continue
    else:
        print( "Variable %s absent in %s, present in %s" % (key, ncfnames[0], ncfnames[1]))
        quit_flag = True

# If missing variable, quit.
if quit_flag:
    quit()


''' Check if variables that are integers and floats are the same '''
for key in keylists[0]:

    # Load the variable
    arr1 = ncfiles[0].variables[key]
    arr2 = ncfiles[1].variables[key]

    # Do not check non-integer variables
    if not np.issubdtype( arr1.dtype, np.number):
        continue

    # Compare the dimensions
    arr1 = np.array(arr1) #+1e-16
    arr2 = np.array(arr2) #+1e-16
    if arr1.shape != arr2.shape:
        quit_flag = True
        print("%s and %s have different %s dimensions" % (ncfnames[0], ncfnames[1], key))
    
    # Compare the variable fields
    diff = np.fabs(arr2 - arr1)
    print("variable %s diff: %e " 
           % (  key, diff.max()))

    diff /= arr1
    diff[np.isnan(diff)] = 0.
    if diff.max() > 1e-5:
        print("%s and %s variable %s diff/arr1: %f " 
               % ( ncfnames[0], ncfnames[1], key, np.nanmax(diff)))
        quit_flag = True

# If files contain variables with different values or different dimensions, quit
if quit_flag:
    quit()


# If both files have the same variables, then print out approval.
print('match')
