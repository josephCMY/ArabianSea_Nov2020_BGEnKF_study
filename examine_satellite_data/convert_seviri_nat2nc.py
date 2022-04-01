''' Script to convert SEVIRI native file into netCDF '''
# IMPORTANT: designed to work with a conda environment that has satpy, numpy and matplotlib installed


# import modules
from satpy import Scene
import sys
import numpy as np
import datetime
import xarray
from netCDF4 import Dataset

# Read in SEVIRIR native file name from command line
fname = 'seviri_native_file/MSG1-SEVI-MSG15-0100-NA-20200315152741.253000000Z-NA.nat'


# Open SEVIRI file
scn = Scene( reader='seviri_l1b_native', filenames=[fname])


# Load all 3km observations an
id_list= scn.all_dataset_ids()
for datid in id_list:
  print(datid)

# Load window bt DataArray from file
scn.load(['IR_108'])
seviri_windowbt_obj = scn['IR_108']
print( seviri_windowbt_obj.values )
seviri_windowbt_obj.to_netcdf(path='trial.nc', format='NETCDF4')


quit()
# Loading bts, lat and lon into numpy arrays
seviri_windowbt_obj = scn['IR_108']
windowbt = seviri_windowbt_obj.values
lon, lat = seviri_windowbt_obj.attrs['area'].get_lonlats()
date = seviri_windowbt_obj.attrs['start_time']



# Plot out bts
fig, axs = plt.subplots(nrows=1,ncols=1, figsize=(8,6))
cnf_clr, cnf_cld = plot_windowbt( lon, lat, windowbt, axs )
axs.set_ylim([-30,30])
axs.set_xlim([50,130])
axs.set_aspect(1)
axs.set_title( 'SEVIRI 10.8 $\mu$ on %s' % date.strftime('%d-%m-%Y %H:%M UTC'))
plt.savefig('trial_seviri_bt_plot.png')
