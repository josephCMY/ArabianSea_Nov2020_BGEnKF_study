''' 
Script to plot out soundings before and after QC
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
date_st=datetime.datetime.strptime( "2017052000", "%Y%m%d%H")
date_ed=datetime.datetime.strptime( "2017060500", "%Y%m%d%H")
date_diff = date_ed - date_st
date_diff = date_diff.total_seconds()
date_diff /= (60*60*24)
zlvl_vlist = ['pres','spd','dir','hgt','temp','dwpt', 'rh']
obs_interval = 60

# Make list of all GTS files
flist = []
datelist=[]
date = date_st + datetime.timedelta(minutes=0)
while date <= date_ed:
  fname = date.strftime( raw_dir+"/obs_gts_%Y-%m-%d_%H:%M:%S.3DVAR"  )
  flist.append(fname)
  datelist.append(date)
  date += datetime.timedelta( minutes=obs_interval)


# Setup container variabe to hold location results
loc_contain ={}
for hh in np.arange(0,24, int( obs_interval/60) ):
    loc_contain["%02dUTC" % hh] = {}
    loc_contain["%02dUTC" % hh]['latlon'] = []
    loc_contain["%02dUTC" % hh]['oo'] = []
loc_contain['all']={}
loc_contain['all']['latlon'] = []
loc_contain['all']['oo'] = []


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
   
    # Seek out line numbers where individual obs entries start
    obsflags = np.zeros( nL, dtype=bool )
    for ll in range( 21, nL ):
        if raw_lines[ll].find( 'FM') != -1:
            obsflags[ll] = True

    # Compute number of entries in obs file
    nObs = len( np.arange(nL)[obsflags] )

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

    # Now iterate through each obs and check for desireable properties
    for oo in np.arange(nObs):

        # Ignore all obs with insufficient zlvls
        if gts_obs[oo]['zlvl'] < 8:
            continue

        # For scenarios where obs is in P coordinates
        # Check if there is at least one obs per 100 hPa interval from 900 hPa to 100 hPa
        cntP, xe = np.histogram( gts_obs[oo]['vert']['pres']['data'],
                               bins=np.linspace( 100*100,900*100, 9)  )


        # For scenarios where obs is in Z coordinates
        # Check if there is at least one temperature obs for height intervals corresponding to 100 hPa intervals
        zbins = [1e3,2e3,3e3,4e3,6e3,8e3,10e3,12e3,16e3]
        cntZ, xe = np.histogram( gts_obs[oo]['vert']['hgt']['data'],
                             bins=zbins )


        # Combined binning effect
        cnt = cntP + cntZ
        # Seek out empty spots
        flags = (cnt==0)
        # Ignore obs with insufficient density
        if np.sum(flags):
            continue



        # All obs that survive to this point will be used.
        for utc in [ '%02dUTC' % hh, 'all']:
            loc_contain[utc]['latlon'].append( [gts_obs[oo]['lon'], gts_obs[oo]['lat']] )
            loc_contain[utc]['oo'].append( oo )



# Plot spatial distribution of obs by time windows.
for utc in loc_contain.keys():

    # Prepare to sort everything
    loc_contain[utc]['latlon']= np.array( loc_contain[utc]['latlon'])

    # Set up map
    bm = Basemap( projection='merc', llcrnrlat=-20., urcrnrlat=10., 
                                      llcrnrlon= 85., urcrnrlon=130., 
                                      lat_ts=0.0, resolution='c')
    bm.drawcoastlines()
    bm.drawparallels( np.arange( -20, 10+1, 10), labels=[True,False,False,False])
    bm.drawmeridians( np.arange(85, 130+1, 10), labels=[False,False,False,True])

    # Aggregate obs
    if loc_contain[utc]['latlon'].shape[0] > 0:
        # Plot out .5 degree aggregated obs counts
        latbins= np.linspace( -22.5,17.5, int((17.5+22.5)/0.5 + 1))
        lonbins= np.linspace( 82.5,132.5, int((132.5-82.5)/0.5 + 1)) 
        hist, xe,ye = np.histogram2d( loc_contain[utc]['latlon'][:,0], loc_contain[utc]['latlon'][:,1], bins = [lonbins, latbins])
        hist = hist/date_diff
        xc = (xe[:-1] + xe[1:])/2.
        yc = (ye[:-1] + ye[1:])/2.
        xx, yy =np.meshgrid(xc,yc)

        if utc != 'all':
            hist = np.ma.masked_where(hist <1., hist)
            bm.pcolormesh( xx,yy, hist.T, latlon=True, vmin=1., vmax=1.01, cmap=plt.cm.bwr)
        else:
            hist = np.ma.masked_where(hist < 1, hist)
            bm.bluemarble()
            bm.pcolormesh( xx,yy, hist.T, latlon=True, vmin=0.5, vmax=3.0, cmap='Reds', zorder=100)
        cb = plt.colorbar(orientation='horizontal')
        cb.ax.set_xlabel('Obs / day')
        plt.savefig('vert_obs_%s_%s.png' % (raw_dir, utc), dpi=200 )
    plt.close()


# Print out locations with at least 2 obs per day
for ii in range( len(xc)):
    for jj in range(len(yc)):
        if hist[ii,jj] >= 2.:
            print( [xc[ii],yc[jj]])


