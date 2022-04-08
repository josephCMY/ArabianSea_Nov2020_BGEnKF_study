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
date_st=datetime.datetime.strptime( "2017053012", "%Y%m%d%H")
date_ed=datetime.datetime.strptime( "2017060103", "%Y%m%d%H")
date_diff = date_ed - date_st
date_diff = date_diff.total_seconds()
date_diff /= (60*60*24)
zlvl_vlist = ['pres','spd','dir','hgt','temp','dwpt', 'rh']
obs_interval = 3*60

# Make list of all GTS files
flist = []
datelist=[]
date = date_st + datetime.timedelta(minutes=0)
while date < date_ed:
  fname = date.strftime( raw_dir+"/obs_gts_%Y-%m-%d_%H:%M:%S.3DVAR"  )
  flist.append(fname)
  datelist.append(date)
  date += datetime.timedelta( minutes=obs_interval)


# Setup container variable to hold density results
amv_density ={}
amv_density['lat edges'] = np.linspace( -22.5,  17.5, int((17.5+22.5)/0.5 + 1))
amv_density['lon edges'] = np.linspace(  82.5, 132.5, int((132.5-82.5)/0.5 + 1))
amv_density['p edges'] = np.linspace( 0,1000,6)*100
nP = len( amv_density['p edges'] )
nLat = len( amv_density['lat edges'] )
nLon = len( amv_density['lon edges'] )
tlist = ['all']
for hh in np.arange(0,24, int( obs_interval/60) ):
    amv_density["%02dUTC" % hh] = np.zeros( [nP, nLat-1, nLon-1], dtype= np.float )
    tlist.append( "%02dUTC" %hh )
amv_density['all'] = np.zeros( [nP, nLat-1, nLon-1], dtype= np.float )


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
    for oo in np.arange(nObs):
        
        ll =  np.arange(nL)[obsflags][oo]

        # Interpret obs header
        info_line = raw_lines[ll]
        info_line.strip()
        platform = info_line[:12].strip()
        date = info_line[13:32].strip()
        name = info_line[33:73].strip()
        lvl = int(info_line[74:80])
        lat = float(info_line[80:92])
        lon = float(info_line[103:115])
        ele = float(info_line[126:138])
        ID = info_line[155:195].strip()

        # Load into obs header into container
        gts_obs[oo] = {}
        gts_obs[oo]['fm'] = platform
        gts_obs[oo]['keep_flag'] = True
        gts_obs[oo]['date'] = date
        gts_obs[oo]['name'] = name
        gts_obs[oo]['zlvl'] = lvl
        gts_obs[oo]['lat'] = lat
        gts_obs[oo]['lon'] = lon
        gts_obs[oo]['ele'] = ele
        gts_obs[oo]['id'] = ID
        gts_obs[oo]['reject'] = False

        # Interpret sfc values (all gts obs have a line for sfc)
        sfc_line = raw_lines[ll+1]
        sfc_line.strip()
        sfc_line = sfc_line.split()
        gts_obs[oo]['sfc'] = {}
        gts_obs[oo]['sfc']['slp'] = {'data': float( sfc_line[0]), 'qc': int( sfc_line[1]), \
                                     'err': float( sfc_line[2]) }
        gts_obs[oo]['sfc']['pw'] = {'data': float( sfc_line[3]), 'qc': int( sfc_line[4]), \
                                     'err': float( sfc_line[5]) }

        # Move on to reading next entry if no vertical levels
        if gts_obs[oo]['zlvl'] == 0:
            continue

        # If obs has vertical levels, we will interpret vertical levels too
        gts_obs[oo]['vert'] = {}
        # Allocate memory
        for vname in zlvl_vlist:
            gts_obs[oo]['vert'][vname] = { 'data': np.zeros( lvl, dtype=np.float ),\
                                           'qc': np.zeros( lvl, dtype=np.int ),  \
                                           'err': np.zeros( lvl, dtype=np.float ) }
        for vname in ['q','u','v']:
            gts_obs[oo]['vert'][vname] = { 'data': np.zeros( lvl, dtype=np.float ),\
                                           'qc': np.zeros( lvl, dtype=np.int ),  \
                                           'err': np.zeros( lvl, dtype=np.float ) }

        # Read in the values
        vert_lines = raw_lines[(ll+2):(ll+2+lvl)]
        for kk in range(lvl):
            subline = vert_lines[kk]
            subline.strip()
            subline=subline.split()
            for ii in range(len(zlvl_vlist)):
                vname = zlvl_vlist[ii]
                gts_obs[oo]['vert'][vname]['data'][kk] = float(subline[ii*3])
                gts_obs[oo]['vert'][vname]['qc'][kk] = int(subline[ii*3+1])
                gts_obs[oo]['vert'][vname]['err'][kk] = float(subline[ii*3+2])

    # Now iterate through each obs and check where the obs belongs to
    for oo in np.arange(nObs):

        # Assess geolocation of obs
        inds = {'lat': -999, 'lon': -999, 'p':-999}
        for coord in ['lat','lon']:
          mask = amv_density[coord+' edges'][:-1] <= gts_obs[oo][coord]
          mask *= ( amv_density[coord+' edges'][1:] > gts_obs[oo][coord] )
          mask = np.array( mask, dtype='bool')
          inds[coord] = int(np.arange( len(amv_density[coord+' edges'])-1)[mask])

        # For each valid wind observation, increment histogram counter
        for kk in range(gts_obs[oo]['zlvl']):

          inds['p'] = -999

          # Exclude invalid obs
          if( gts_obs[oo]['vert']['spd']['data'][kk] == -888888 
                  or gts_obs[oo]['vert']['dir']['data'][kk] == -888888 ):
              continue

          # Determine the pressure index
          mask = amv_density['p edges'][:-1] <= gts_obs[oo]['vert']['pres']['data'][kk]
          mask *= (amv_density['p edges'][1:] > gts_obs[oo]['vert']['pres']['data'][kk] )
          inds['p'] = int(np.arange( len( amv_density['p edges']-1 )-1)[mask])

          # Increment the appropriate 3D histogram
          for utc in [ '%02dUTC' % hh, 'all']:
            amv_density[utc][inds['p'], inds['lat'], inds['lon']] += 1 
            

# Plot spatial distribution of obs by time windows and pressure layer
xc = (amv_density['lon edges'][:-1] + amv_density['lon edges'][1:])/2.
yc = (amv_density['lat edges'][:-1] + amv_density['lat edges'][1:])/2.
xx, yy =np.meshgrid(xc,yc)

for utc in ['all']:

  # Set up figure
  fig, axs = plt.subplots( nrows=3, ncols=2, figsize=(9,9))
  axs = [val for sub in axs for val in sub ]
  axs[-1].remove()

  # For each pressure layer
  for kk in range( nP-1 ):
    ax = axs[kk]

    # Set up map
    bm = Basemap( projection='merc', llcrnrlat=-20., urcrnrlat=10., 
                                      llcrnrlon= 85., urcrnrlon=130., 
                                      lat_ts=0.0, resolution='c', ax = ax)
    bm.drawcoastlines()
    bm.drawparallels( np.arange( -20, 10+1, 5), labels=[True,False,False,False])
    bm.drawmeridians( np.arange(85, 130+1, 5), labels=[False,False,False,True])

    hist = amv_density[utc][kk,:,:]/(date_diff*8)
    hist = np.ma.masked_where(hist < 1, hist)

    pc = bm.pcolormesh( xx,yy, hist, latlon=True, vmin=1., vmax=10, cmap="jet")
    ax.set_title( 'P $\in$ [%d, %d) hPa' % ( int(amv_density['p edges'][kk]/100), 
                                             int(amv_density['p edges'][kk+1]/100) ))

  fig.subplots_adjust( right=0.85,left=0.05, wspace=0.15, hspace=0.3, bottom=0.05, top=0.95)
  cb_ax = fig.add_axes( [0.88,0.2, 0.05, 0.6])
  cb = fig.colorbar( pc, cax=cb_ax, orientation='vertical')
  cb.ax.set_ylabel('Obs / 3 hrs')
  plt.savefig('amv_obs_%s_%s.png' % (raw_dir, utc), dpi=200 )
  plt.close()


