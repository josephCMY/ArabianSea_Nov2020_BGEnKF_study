''' 
Script to plot out the AMV density and frequency
'''

import numpy as np
import datetime
from mpl_toolkits.basemap import Basemap
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import sys

# User controls
raw_dir = sys.argv[1]
date_st=datetime.datetime.strptime( "2017053000", "%Y%m%d%H")
date_ed=datetime.datetime.strptime( "2017053001", "%Y%m%d%H")
date_diff = date_ed - date_st
date_diff = date_diff.total_seconds()
date_diff /= (60*60*24)
zlvl_vlist = ['pres','spd','dir','hgt','temp','dwpt', 'rh']
obs_interval = 60

# Make list of all GTS files
flist = []
datelist=[]
date = date_st + datetime.timedelta(minutes=0)
while date < date_ed:
  fname = date.strftime( raw_dir+"/obs_gts_%Y-%m-%d_%H:%M:%S.3DVAR"  )
  flist.append(fname)
  datelist.append(date)
  date += datetime.timedelta( minutes=obs_interval)
print( flist)

# Setup container variable to hold density results
amv_frequency ={}
amv_frequency['lat edges'] = np.linspace( -20,  40, 60*4+1)
amv_frequency['lon edges'] = np.linspace(  60, 140, 80*4+1 )
tmp_lat = np.zeros( [80000], dtype='float') +np.nan
tmp_lon = np.zeros( [80000], dtype='float') +np.nan

# For each GTS file and its corresponding date,
for ff in range(len(flist)):

#    print( "\n\n\nDealing with %s" % (raw_flist[ff]))
    # Read in all the lines in gts file
    f = open( flist[ff], 'r')
    raw_lines = [ line for line in f ]
    nL = len(raw_lines)
    f.close()

    # Call observation file time
    hh = datelist[ff].hour

    # Useful lines to preserve
    header_invariants = raw_lines[5:21]

    # Define container for the raw and qc'ed obs
    gts_obs = {}            # Raw obs
    out_obs = {}            # QC'ed obs
   
    # Seek out line numbers where individual AMV obs entries start
    obsflags = np.zeros( nL, dtype=bool )
    for ll in range( 21, nL ):
        if raw_lines[ll].find( 'FM-88 SATOB') != -1:
            obsflags[ll] = True

    # Compute number of entries in obs file
    nObs = len( np.arange(nL)[obsflags] )
    print( "# of AMV found: %d" % nObs )

    # For each obs entry, read in all data
    tmp_lat[:] = np.nan; tmp_lon[:] = np.nan
    for oo in np.arange(nObs):
        
        ll =  np.arange(nL)[obsflags][oo]

        # Interpret obs header to extract lat lon info
        info_line = raw_lines[ll]
        info_line.strip()
        tmp_lat[oo] = float(info_line[80:92])
        tmp_lon[oo] = float(info_line[103:115])

    # Count AMV obs
    if ff == 0:
      amv_frequency['count'], xedges, yedges \
              = np.histogram2d( tmp_lat[:nObs], tmp_lon[:nObs], 
                          bins = ( amv_frequency['lat edges']
                                   , amv_frequency['lon edges'] ),
#                          range= [ [amv_frequency['lat edges'].min(), 
#                                    amv_frequency['lat edges'].max()],
#                                   [amv_frequency['lon edges'].min(),
#                                    amv_frequency['lon edges'].max()]
#                                 ],
                          density = False )
    else:
      tmp, xedges, yedges \
              = np.histogram2d( tmp_lat[:nObs], tmp_lon[:nObs], 
                          bins = ( amv_frequency['lat edges']
                                   , amv_frequency['lon edges'] ),
#                          range= [ [amv_frequency['lat edges'].min(), 
#                                    amv_frequency['lat edges'].max()],
#                                   [amv_frequency['lon edges'].min(),
#                                    amv_frequency['lon edges'].max()]
#                                 ],
                          density = False )
      amv_frequency['count'] += tmp


            



# Set up figure
fig = plt.figure( figsize=(5,4))
bm = Basemap( projection ='cyl', lat_0=0, lon_0=90, resolution='l', 
               llcrnrlat = -20., urcrnrlat=40.,
               llcrnrlon =  60, urcrnrlon =140.)
bm.drawcoastlines()
bm.drawparallels( np.arange( -20, 40+1, 10), labels=[True,False,False,False])
bm.drawmeridians( np.arange(60, 140+1, 10), labels=[False,False,False,True])

pc = bm.pcolormesh( amv_frequency['lon edges'], amv_frequency['lat edges'],
                    amv_frequency['count']*1./len(flist), vmin = 0.0, vmax = 1.,
                    cmap = 'jet' )
#bm.ax.set_title('Hourly AMV from %s to %s'
#             % (date_st.strftime('%HUTC %d-%b'), date_ed.strftime('%HUTC %d-%b')))

fig.subplots_adjust( right=0.85,left=0.05, wspace=0.15, hspace=0.3, bottom=0.05, top=0.95)
cb_ax = fig.add_axes( [0.88,0.2, 0.05, 0.6])
cb = fig.colorbar( pc, cax=cb_ax, orientation='vertical')
cb.ax.set_ylabel('Obs / hr')
plt.savefig('amv_obs_%s.png' % (raw_dir), dpi=200 )
plt.close()


