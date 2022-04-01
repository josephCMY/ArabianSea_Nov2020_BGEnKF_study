''' Experimental script to plot out seviri bt from native data format '''
# IMPORTANT: designed to work with a conda environment that has satpy, numpy and matplotlib installed


# import modules
from satpy import Scene
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import sys
import numpy as np
import matplotlib.colors as mcolors
import datetime


# Read in SEVIRIR native file name from command line
fname = 'MSG1-SEVI-MSG15-0100-NA-20200318065740.443000000Z-NA.nat'




''' Special things to plot out OLR-BT'''
cmap_data = [(0.0, 0.7529411911964417, 0.0),
             (0.501960813999176, 0.8784313797950745, 0.0),
             (1.0, 1.0, 0.0),
             (1.0, 0.6274510025978088, 0.0),
             (1.0, 0.0, 0.0),
             (1.0, 0.125490203499794, 0.501960813999176),
             (0.9411764740943909, 0.250980406999588, 1.0),
             (0.501960813999176, 0.125490203499794, 1.0)][::-1]
cloud_cmap = mcolors.LinearSegmentedColormap.from_list('ppt',cmap_data)


# Special function to plot Window-BT
# ----------------------------------
def plot_windowbt( lon2d, lat2d, bt2d, ax ):

  # Plot out clr BTs
  clr_bts = bt2d*1.
  clr_bts[ bt2d < 280 ] = np.nan
  cnf_clr = ax.contourf( lon2d, lat2d, clr_bts, np.linspace(280,300,11),
                         cmap = 'binary', extend='max')

  # Plot out cld BTs
  cld_bts = bt2d*1.
  cld_bts[ bt2d > 264 ] = np.nan
  cnf_cld = ax.contourf( lon2d, lat2d, cld_bts, np.linspace(200,280,11),
                         cmap = cloud_cmap, extend='min')


  return cnf_clr, cnf_cld 




# Open SEVIRI file
scn = Scene( reader='seviri_l1b_native', filenames=[fname])

# Load window bt DataArray from file
scn.load(['IR_108'])

# Loading bts, lat and lon into numpy arrays
seviri_bt_obj = scn['IR_108']
bt = seviri_bt_obj.values
lon, lat = seviri_bt_obj.attrs['area'].get_lonlats()
date = seviri_bt_obj.attrs['start_time']
print( lon.shape)
print(date)

# Plot out bts
fig, axs = plt.subplots(nrows=1,ncols=1, figsize=(8,6))
cnf_clr, cnf_cld = plot_windowbt( lon, lat, bt, axs )
axs.set_ylim([-30,30])
axs.set_xlim([50,100])
axs.set_aspect(1)
plt.savefig('trial_seviri_bt_plot.png')
