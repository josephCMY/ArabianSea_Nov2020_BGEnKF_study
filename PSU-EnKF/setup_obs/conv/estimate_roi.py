''' 
Script to read in obs_gts files, then exclude the rows that contain validation obs and construct a new obs_gts file
'''

from metpy import calc as met
from metpy.units import units
import numpy as np
import datetime
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from netCDF4 import Dataset as ncopen
from math import sin, cos, atan2, sqrt

# User controls
gts_dir = "gts_assimilate"
era5_dir = "era5_data"
date_st=datetime.datetime.strptime( "2017053100", "%Y%m%d%H")
date_ed=datetime.datetime.strptime( "2017060900", "%Y%m%d%H")
date_diff = date_ed - date_st
date_diff = date_diff.total_seconds()
date_diff /= (60*60*24)
typestr_list = ['SYNOP','METAR','SHIP','BUOY','BOGUS','TEMP','AMDAR','AIREP','TAMDAR',\
                'PILOT','SATEM','SATOB','GPSPW','GPSZD','GPSRF','GPSEP','SSMT1','SSMT2',\
                'TOVS', 'QSCAT','PROFL','AIRSR','OTHER']
zlvl_vlist = ['pres','spd','dir','hgt','temp','dwpt', 'rh']


# Haversine formula (only put in lat lon in terms of radians!)
def haversine_distance( lat1, lon1, lat2, lon2 ):
    R = 6371.
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = sin(dlat/2)*sin(dlat/2) + cos(lat1)*cos(lat2) * sin( dlon/2) * sin( dlon/2)
    c = 2 * atan2( sqrt(a), sqrt(1.-a))
    return R * c



# Make list of all GTS files
gts_flist = []
date_list = []
date = date_st + datetime.timedelta(minutes=0)
while date <= date_ed:
    gts_fname = date.strftime( gts_dir+"/obs_gts_%Y-%m-%d_%H:%M:%S.3DVAR"  )
    gts_flist.append(gts_fname)
    date_list.append(date)
    date += datetime.timedelta( minutes=180)


# Nearest neighbour distance container
nndist = {}
for tt in typestr_list:
    nndist[tt] = []


# For each GTS file and its corresponding date,
for ff in range(len(gts_flist)):

    print( "\nDealing with %s" % (gts_flist[ff]))
    # Read in all the lines in raw file
    f = open( gts_flist[ff], 'r')
    raw_lines = [ line for line in f ]
    nL = len(raw_lines)
    f.close()

    # Define container for the gts obs
    gts_obs = {}
   
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
        gts_obs[oo]['paired'] = False              # Flag to help with nearest distance computation

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


    # For each observation type, compute nearest neighbour distances
    for tt in typestr_list:

        # Iterate over obs
        for o1 in np.arange(nObs):

            # Ignore irrelevant obs
            if gts_obs[o1]['fm'].find( tt ) == -1:
                continue

            # Ignore obs that were alrdy been paired in nearest distance approach
            if gts_obs[o1]['paired']:
                continue

            # First distance
            dist = 1e10          

            # Pair starting index
            o2_nn = -999

            # iterate ovr Obs
            for o2 in np.arange(o1+1, nObs):
                
                # Ignore irrelevant obs
                if gts_obs[o2]['fm'].find(tt) == -1:
                    continue

                # Ignore paired obs
                if gts_obs[o2]['paired']:
                    continue

                # Compute distance between the 2 obs.
                dist_test = haversine_distance( gts_obs[o1]['lat'] * 3.1415/180.,
                                                gts_obs[o1]['lon'] * 3.1415/180.,
                                                gts_obs[o2]['lat'] * 3.1415/180.,
                                                gts_obs[o2]['lon'] * 3.1415/180.)

                # Check if distance is smaller
                if dist > dist_test:
                    dist = dist_test
                    o2_nn = o2

            # Flip the flag!
            if dist < 1e9:
                gts_obs[o1]['paired'] = True 
                gts_obs[o2_nn]['paired'] = True
                nndist[tt].append(dist)

            
# Now that we have all the distances computed, plot them out in histogram

# Count number of obstypes with pairs of obs
cnt_obs = 0
for tt in typestr_list:
     if len(nndist[tt]) > 0:
         cnt_obs +=1

# Decide on how many rows and columns
ncol = 3
nrow = int(cnt_obs/ncol)
if nrow < cnt_obs/(1.*ncol):
    nrow +=1

# Setup figure.
fig = plt.figure( figsize=(ncol*3, nrow*3.5))
axes = [ [fig.add_subplot(nrow,ncol,i*ncol + j+1) for j in range(ncol)] for i in range(nrow)]
axes = [ val for sublist in axes for val in sublist]

# Now plot histograms
cc = 0
for t in np.arange(len(typestr_list)):
    tt = typestr_list[t]
    if len(nndist[tt]) == 0:
        continue
    edges = [ 0, np.percentile( nndist[tt], 95.)]
    hist, xe = np.histogram( nndist[tt], bins=np.linspace(edges[0],edges[1], 21))
    dx = abs(xe[1]-xe[0])
    axes[cc].bar( xe[:-1], hist, width=dx, align='edge', color='gray')
    mean = np.mean( nndist[tt])
    pctl = [ np.percentile( nndist[tt], pp) for pp in [25., 50., 75.]]
#    axes[cc].axvline( x=mean, color='k', linestyle = '-', linewidth=2.0)
#    for p in [0,2]:
#        axes[cc].axvline( x=pctl[p], color='r', linestyle='--')
#    axes[cc].axvline( x=pctl[1], color='k', linestyle='--', linewidth=2.0)
    axes[cc].set_title( "%s\n(avg: %d, med: %d)" % (tt, int(mean), int(pctl[1])) )
    axes[cc].set_xlabel( 'nn dist (km)')
    axes[cc].set_ylabel( 'frequency')
    cc += 1

fig.subplots_adjust( hspace=0.6, wspace=0.4)

plt.savefig( 'nn_distances.png', dpi=100)

   







            





