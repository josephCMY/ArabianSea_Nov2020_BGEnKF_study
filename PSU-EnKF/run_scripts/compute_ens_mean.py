''' Computation of ensemble mean for deterministic forecast
    NetCDF operator NCO's ensemble mean computation is too slow.
    Will use Python to do the job.
    Requires the presence of a template file '''

import numpy as np
import sys
from netCDF4 import Dataset as ncopen
import glob

# Inputs
fname_avg = sys.argv[1]
fname_ens_prefix = sys.argv[2]
n_ens = int(sys.argv[3])

# Open up the average file and zero all non-character variables
favg = ncopen( fname_avg, 'a' )
vlist = favg.variables.keys()
for vname in vlist:
    if favg.variables[vname].dtype == 'float32':
        favg.variables[vname][:] *= 0.

# Now put in some ensemble means
fname_ens = [ "%s_%03d" % ( fname_ens_prefix, ee+1 ) for ee in range(n_ens) ]
for ee in range( n_ens):
    print('  Including mem %02d' % (ee+1))
    fens = ncopen( fname_ens[ee], 'r')
    for vname in vlist:
        if favg.variables[vname].dtype == 'float32':
            favg.variables[vname][:] += fens.variables[vname][:]/n_ens
    fens.close() 

# Close file
favg.close()
