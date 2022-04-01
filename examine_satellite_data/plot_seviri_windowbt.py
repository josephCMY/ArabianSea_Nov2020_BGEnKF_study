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
fname = sys.argv[1] 



''' Special things to plot out OLR-BT'''
cmap_data = [ (0.3137255012989044, 0.8156862854957581, 0.8156862854957581),
             (0.0, 1.0, 1.0),
             (0.0, 0.8784313797950745, 0.501960813999176),
             (0.0, 0.7529411911964417, 0.0),
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
  clr_bts[ bt2d < 285 ] = np.nan
  cnf_clr = ax.contourf( lon2d, lat2d, clr_bts, np.linspace(285,305,11),
                         cmap = 'binary', extend='max')

  # Plot out cld BTs
  cld_bts = bt2d*1.
  cld_bts[ bt2d > 256 ] = np.nan
  cnf_cld = ax.contourf( lon2d, lat2d, cld_bts, np.linspace(200,280,11),
                         cmap = cloud_cmap, extend='min')


  return cnf_clr, cnf_cld 


# Open SEVIRI file
scn = Scene( reader='seviri_l1b_native', filenames=[fname])
#id_list= scn.all_dataset_ids()
#for datid in id_list:
#    print( datid)

# Load window bt DataArray from file
scn.load(['IR_108'])


# Loading bts, lat and lon into numpy arrays
seviri_windowbt_obj = scn['IR_108']
windowbt = seviri_windowbt_obj.values
lon, lat = seviri_windowbt_obj.attrs['area'].get_lonlats()
date = seviri_windowbt_obj.attrs['start_time']


# Plot out bts
fig, axs = plt.subplots(nrows=1,ncols=1, figsize=(8,6))
cnf_clr, cnf_cld = plot_windowbt( lon, lat, windowbt, axs )
axs.set_ylim([-30,30])
axs.set_xlim([50,115])
axs.set_aspect(1)
axs.set_title( 'SEVIRI 10.8 $\mu$m on %s' % date.strftime('%d-%b-%Y %H:%M UTC'))
plt.savefig('ir10.8_%s.png' % date.strftime( "%Y%m%d%H%M") )
plt.close()
