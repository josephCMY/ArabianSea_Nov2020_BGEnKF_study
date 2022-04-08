''' 
Functions used to compare obs against ERA5 (mean product)
'''

import pickle
import numpy as np
from netCDF4 import Dataset as ncopen
import datetime
from mpl_toolkits.basemap import Basemap
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import wrf
import time
from scipy.interpolate import interp1d, griddata
import sys
import math


'''
# User controls
ctrl_params = {}
ctrl_params['datestr'] = sys.argv[2]
ctrl_params['date'] = datetime.datetime.strptime( datestr , "%Y%m%d%H%M")
ctrl_params['obs_vlist'] = ['pres','wspd','wdir','z','t','dwpt', 'rh']
ctrl_params['era5_vlist'] = ['u','v','q','t','z']
ctrl_params['nSamp'] = 1000                 # number of monte carlo samples to interpolate observations

# Prepare observation fname
obs_gts = {}
obs_gts['obs_raw_fname'] = "SPECIFY PATH HERE"
obs_gts['era5_fname'] = "SPECIFY PATH HERE"
obs_gts['obs_out_fname'] = "SPECIFY PATH HERE"
'''



'''
Function to read in control parameters
'''
def read_ctrl_params( ctrl_params ):
  # Read in the control parameters
  datestr       = ctrl_params['datestr']
  date          = ctrl_params['date']
  obs_vlist     = ctrl_params['obs_vlist']
  era5_vlist    = ctrl_params['era5_vlist']
  nSamp         = ctrl_params['nSamp']
  return datestr, date, obs_vlist, era5_vlist, nSamp



'''
Function to draw true white noise.
Used for propagating observation errors from one variable to another.
Also used to propagate observation errors due to interpolation.
'''
def draw_true_white_noise( nSamp, nVar ):
  white_noise = np.random.normal( loc=0., scale=1., size=( nSamp, nVar ) )
  white_noise -= np.mean( white_noise , axis=0)
  cov = np.cov( white_noise, rowvar = False, ddof=1)
  decorr_matrix= np.linalg.inv(np.linalg.cholesky( cov ))
  white_noise = np.matrix( decorr_matrix ) * np.matrix(white_noise.T)
  white_noise = np.array(white_noise.T)
  return white_noise



''' 
Function to read in GTS file containing validation
Will read in the GTS file and then convert to model variables using a Monte-Carlo
method. 
'''
def read_obs_gts( ctrl_params, obs_gts ):

  # Interprete ctrl_params
  datestr, date, obs_vlist, era5_vlist, nSamp \
  = read_ctrl_params( ctrl_params )


  # Read in all the lines in gts file
  f = open( obs_gts['obs_raw_fname'] , 'r')
  raw_lines = [ line for line in f ]
  nL = len(raw_lines)
  f.close()
  
  # Seek out line numbers where individual obs entries start
  obsflags = np.zeros( nL, dtype=bool )
  for ll in np.arange(  nL ):
    if raw_lines[ll].find( 'FM-') != -1:
      obsflags[ll] = True
  
  # Compute number of entries in obs file
  nObs = len( np.arange(nL)[obsflags] )
  obs_gts['nObs'] = nObs
  print( "%05d entries found for %s" % (nObs, date.strftime("%b %d, %H UTC")) )
    
  oolistflag = np.zeros( nObs, dtype=np.int)
  
  # For each obs entry, read in all data
  for oo in np.arange(nObs):
         
    ll =  np.arange(nL)[obsflags][oo]
  
    # Interpret obs header
    info_line = raw_lines[ll]
    info_line.strip()
    platform = info_line[:12].strip()
    datestr = info_line[13:32].strip()
    name = info_line[33:73].strip()
    lvl = int(info_line[74:80])
    lat = float(info_line[80:92])
    lon = float(info_line[103:115])
    ele = float(info_line[126:138])
    ID = info_line[155:195].strip()
    
    # Check for duplicate site observations
    oolistflag[oo] += 1
    if oolistflag[oo] > 1:
      print( "\nException case")
      print( "duplicate observation at site: %5.1fN, %5.1fE" % (lat, lon) )
      print( "Quitting" )
      quit()
  
    # Load into obs header into container
    obs_gts[oo] = {}
    obs_gts[oo]['obs'] = {}
    obs_gts[oo]['obs']['fm'] = platform
    obs_gts[oo]['obs']['keep_flag'] = True
    obs_gts[oo]['obs']['date'] = datestr
    obs_gts[oo]['obs']['name'] = name
    obs_gts[oo]['obs']['zlvl'] = lvl
    obs_gts[oo]['obs']['lat'] = lat
    obs_gts[oo]['obs']['lon'] = lon
    obs_gts[oo]['obs']['ele'] = ele
    obs_gts[oo]['obs']['id'] = ID
    obs_gts[oo]['ignore'] = False

    # Interpret sfc values (all gts obs have a line for sfc)
    sfc_line = raw_lines[ll+1]
    sfc_line.strip()
    sfc_line = sfc_line.split()
    obs_gts[oo]['obs']['sfc'] = {}
    obs_gts[oo]['obs']['sfc']['slp'] = { 'data': float( sfc_line[0]), 'qc': int( sfc_line[1]), \
                                         'err': float( sfc_line[2]) }
    obs_gts[oo]['obs']['sfc']['pw'] = {'data': float( sfc_line[3]), 'qc': int( sfc_line[4]), \
                                             'err': float( sfc_line[5]) }
  
    # If obs has vertical levels, we will interpret vertical levels too
    obs_gts[oo]['obs']['vert'] = {}
    # Allocate memory
    for vname in obs_vlist:
      obs_gts[oo]['obs']['vert'][vname] = { 'data': np.zeros( lvl, dtype=np.float ),\
                                               'qc': np.zeros( lvl, dtype=np.int ),  \
                                                'err': np.zeros( lvl, dtype=np.float ) }
    for vname in ['q','u','v']:
      obs_gts[oo]['obs']['vert'][vname] = { 'data': np.zeros( lvl, dtype=np.float ),\
                                               'qc': np.zeros( lvl, dtype=np.int ),  \
                                               'err': np.zeros( lvl, dtype=np.float ) }
  
    # Read in the values
    vert_lines = raw_lines[(ll+2):(ll+2+lvl)]
    for kk in np.arange(lvl):
      subline = vert_lines[kk]
      subline.strip()
      subline=subline.split()
      for ii in np.arange(len(obs_vlist)):
        vname = obs_vlist[ii]
        obs_gts[oo]['obs']['vert'][vname]['data'][kk] = float(subline[ii*3])
        obs_gts[oo]['obs']['vert'][vname]['qc'][kk] = int(subline[ii*3+1])
        obs_gts[oo]['obs']['vert'][vname]['err'][kk] = float(subline[ii*3+2])
        if float(subline[ii*3]) == -888888:
          obs_gts[oo]['obs']['vert'][vname]['data'][kk] = np.nan
  
    ''' Use Monte Carlo approach to construct q, u and v '''
    white_noise = draw_true_white_noise( nSamp, 5*lvl )
    # Distribute white noise to form particles
    samples = { 't': white_noise[:, :lvl], 'wspd': white_noise[:,lvl:(2*lvl)],
                'wdir': white_noise[:, (2*lvl):(3*lvl)], 
                'rh': white_noise[:, (3*lvl):(4*lvl)],
                'pres': white_noise[:, (4*lvl):(5*lvl)]}
    # Convert white noise to relevant samples of various quantities
    for vname in ['t','wspd','wdir','rh', 'pres']:
        samples[vname] *= obs_gts[oo]['obs']['vert'][vname]['err']
        samples[vname] += obs_gts[oo]['obs']['vert'][vname]['data']
    # Construct U and V
    samples['u'] = samples['wspd']*np.sin( samples['wdir']*math.pi/180.) * -1
    samples['v'] = samples['wspd']*np.cos( samples['wdir']*math.pi/180.) * -1
    # Construct QVAPOR
    es = 611. * np.exp( 6808*( 1/273. -1/samples['t'])
                        - 5.09*np.log( samples['t']/273. ) )
    e = es * samples['rh'] / 100.
    pd = samples['pres'] - e
    rho_d = pd / ( 287. * samples['t'] )
    rho_v = e / ( 461. * samples['t'] )
    samples['q'] = rho_v / rho_d
    # Save statistics of converted obs quantities
    for vname in ['u','v','q']:
      obs_gts[oo]['obs']['vert'][vname]['data'][:] = np.mean( samples[vname], axis=0)
      obs_gts[oo]['obs']['vert'][vname]['err'][:] = np.std( samples[vname], axis=0, ddof=1 )
  

  # End of read GTS function
  return


'''
Function to construct simulated soundings from ERA5 reanalysis
'''
def era5_soundings( ctrl_params, obs_gts ):

  # Interprete ctrl_params
  datestr, date, obs_vlist, era5_vlist, nSamp \
  = read_ctrl_params( ctrl_params )
  nObs = obs_gts['nObs']

  # Allocate memory
  for oo in np.arange(nObs):
    obs_gts[oo]['era5'] = {}
    for vname in era5_vlist: 
      obs_gts[oo]['era5'][vname] = np.zeros( obs_gts[oo]['obs']['zlvl'] )

  # Going into each ensemble member and extracting relevant variables at relevant locations
  era5f = ncopen( obs_gts['era5_fname'], 'r')

  # Import coordinates and pin down locations closest to the obs site
  lat = np.array( era5f.variables['latitude'] )
  lon = np.array( era5f.variables['longitude'] )
  era5_pres= np.array( era5f.variables['level'] ) * 100
#  print( era5_pres)
  lon, lat = np.meshgrid( lon, lat )
  dlat = np.fabs(np.mean(lat[1:] - lat[:-1]) )
  dlon = np.fabs(np.mean(lon[:,1:] - lon[:,:-1] ))
#  print([dlat, dlon])
  # Pinning down locations
  obs_neighbourhood = {}
  for oo in np.arange(nObs):
    obs_neighbourhood[oo] = {}
    # Generating mask to pin down location
    mask = (lat <= obs_gts[oo]['obs']['lat'] + dlat)
    mask = mask * (lat > obs_gts[oo]['obs']['lat']-dlat)
    mask = mask * (lon <= obs_gts[oo]['obs']['lon']+dlon)
    mask = mask * (lon > obs_gts[oo]['obs']['lon']-dlon)
    # Generating location indices
    mask1d = np.sum( mask , axis=0) > 0
    indlon = np.arange( len(mask1d) )[mask1d]
    mask1d = np.sum( mask , axis=1) > 0
    indlat = np.arange( len(mask1d) )[mask1d]
    obs_neighbourhood[oo]['indlat'] = indlat
    obs_neighbourhood[oo]['indlon'] = indlon
    # Special flag to ignore obs that are outside checking area
    if np.sum(mask ) == 0:
      obs_gts[oo]['ignore'] = True
      print( 'Ignoring obs @ %6.2fN %6.2fE' % 
             (obs_gts[oo]['obs']['lat'], obs_gts[oo]['obs']['lon'] ) )

  # Import data in the neighbourhood of the observations
  for vname in era5_vlist:
    arr = np.array( era5f.variables[vname] )[0,:]
#    print( arr.shape)
    for oo in np.arange(nObs):
      # Skip over obs outside domain
      if obs_gts[oo]['ignore']:
        continue
      indlat = obs_neighbourhood[oo]['indlat']
      indlon = obs_neighbourhood[oo]['indlon']
      obs_neighbourhood[oo][vname] = arr[:,indlat][:,:,indlon]
#        print( obs_neighbourhood[oo][ee][vname].shape)
  era5f.close()

  # Compute horizontal interpolation weights using griddata
  xy_coord = []
  for oo in np.arange(nObs):
    # Skip over obs outside domain
    if obs_gts[oo]['ignore']:
      continue
    # Rearranging things to use griddata
    indlat = obs_neighbourhood[oo]['indlat']
    indlon = obs_neighbourhood[oo]['indlon']
    obsloc = np.array([obs_gts[oo]['obs']['lat'], obs_gts[oo]['obs']['lon']])
    latloc = lat[ indlat][:, indlon ]
    lonloc = lon[ indlat][:, indlon ]
    latloc1d = np.reshape( latloc, np.prod( latloc.shape) )
    lonloc1d = np.reshape( lonloc, np.prod( lonloc.shape) )
    comb_coord = np.array( [latloc1d, lonloc1d] ).T
    # Array to hold griddata weights
    obs_neighbourhood[oo]['weights'] = np.zeros( [len(indlat), len(indlon)] )
    for i in range( len(indlat) ):
      for j in range( len(indlon) ):
        tmp = np.zeros( [len(indlat), len(indlon)] )
        tmp[i,j] = 1.
        tmp = np.reshape( tmp, np.prod(tmp.shape) )
        ww = griddata( comb_coord, tmp, obsloc )
        obs_neighbourhood[oo]['weights'][i,j] = ww
    # Check if applying the weights back out the correct lat and lon
#    print( obs_neighbourhood[oo]['weights'] )
    memlat = np.sum(latloc * obs_neighbourhood[oo]['weights'])
    memlon = np.sum(lonloc * obs_neighbourhood[oo]['weights'])
    if np.fabs(memlat - obsloc[0] ) > dlat:
        print( 'ERROR: Obs and interpolated lats dont match.')
        print( memlat)
        print( obsloc[0] )
        quit()
    if np.fabs(memlon - obsloc[1] ) > dlon:
        print( 'ERROR: Obs and interpolated lons dont match.')
        quit()

  #   Now perform the interpolation accordingly
  for oo in np.arange(nObs):
    # Skip over obs outside domain
    if obs_gts[oo]['ignore']:
      continue

    indlat = obs_neighbourhood[oo]['indlat']
    indlon = obs_neighbourhood[oo]['indlon']

    # Interpolate using whichever vertical coordinates that the observation uses
    for kk in np.arange( obs_gts[oo]['obs']['zlvl'] ):

      # Check if hgt or pressure is available
      flag_z = ( obs_gts[oo]['obs']['vert']['z']['data'][kk]  != np.nan)
      flag_p = ( obs_gts[oo]['obs']['vert']['pres']['data'][kk] != np.nan)

      # Select appropriate vertical coordinate system and interpolate
      if flag_p:
        era5_vert_coord = np.zeros( [2,2,len(era5_pres)] )
        era5_vert_coord[:,:,:] += era5_pres[:]
        obs_lvl = obs_gts[oo]['obs']['vert']['pres']['data'][kk]
      elif flag_z:
        # hgt coordinate needs to be inverted as the vert coordinate in ERA5 is 
        # in terms of ascending pressure.
        era5_vert_coord = -obs_neighbourhood[oo]['z']/9.81
        obs_lvl = -obs_gts[oo]['obs']['vert']['z']['data'][kk]

      # Now interpolate to the obs site
      # Go into each variable
      for vname in era5_vlist:
        # Interpolate neighbourhood in the vertical and horizontal
        for i in range(len(indlat)):
          for j in range(len(indlon)):
            raw_data = obs_neighbourhood[oo][vname][:,i,j]
            # Vertical interpolation will not output NaN unless the raw_data
            # itself contains NaN.
            # If obs_lvl is beyond era_vert_coord, we will output the leftmost or
            # rightmost value in raw_data.
            tmp = np.interp( obs_lvl, era5_vert_coord[i,j], raw_data )
            if obs_gts[oo]['obs']['zlvl'] > 20 and ctrl_params['check_interp'] \
               and obs_lvl < 100000 and vname == 't' and i ==0 and j==0:
              if flag_p:
                check_p = np.interp( obs_lvl, era5_vert_coord[i,j], era5_pres)
                print( 'obs pres, interp pres  =  %6.1f hPa, %6.1f hPa' % (check_p/100, obs_lvl/100) )
              elif flag_z:
                check_z = np.interp( obs_lvl, era5_vert_coord[i,j], era5_vert_coord[i,j])
                print( 'obs hgt, interp hgt  =  %6.1f m, %6.1f m' % (check_z, obs_lvl) )
              if kk == 10:
                print('Exiting vertical interpolation sanity check')
                quit()
            
            # Horizontal interpolation
            obs_gts[oo]['era5'][vname][kk] = tmp *obs_neighbourhood[oo]['weights'][i,j]
 
  
  # Return after all sites have been contructed from ERA5 via interpolation
  return





'''
Function to compare the GTS obs against ERA5
'''
def compare_era5_obs( ctrl_params, obs_gts ):
  
  # Interpret ctrl_params
  datestr, date, obs_vlist, era5_vlist, nSamp \
  = read_ctrl_params( ctrl_params )
  nObs = obs_gts['nObs']

  # Iterate comparison across all obs
  for oo in np.arange(nObs):
    # Skip over obs outside domain
    if obs_gts[oo]['ignore']:
      continue

    for kk in np.arange( obs_gts[oo]['obs']['zlvl'] ):
      for vname in era5_vlist:
        # Skip over if entry's variable is NaN
        if obs_gts[oo]['obs']['vert'][vname]['data'][kk] == np.nan:
          continue
        # Compare GTS against ERA5. 
        diff = np.fabs( obs_gts[oo]['obs']['vert'][vname]['data'][kk] 
                        - obs_gts[oo]['era5'][vname][kk] )
        # Reject if distance is more than 2 times the obs sigma
        reject = diff > 3.0 * obs_gts[oo]['obs']['vert'][vname]['err'][kk]
        # Now reject accordingly
        if reject:
          if vname == "ua" or "va":
            obs_gts[oo]['obs']['vert']['wspd']['data'][kk] = np.nan
            obs_gts[oo]['obs']['vert']['wdir']['data'][kk] = np.nan
          elif vname == "q":
            obs_gts[oo]['obs']['vert']['rh']['data'][kk] = np.nan
            obs_gts[oo]['obs']['vert']['dwpt']['data'][kk] = np.nan
          elif vname == "t":
            obs_gts[oo]['obs']['vert']['t']['data'][kk] = np.nan
        # Move onto next variable in entry.
      # Move onto next level in entry

    # Identifying all useless vertical obs lines
    flag_nan= np.ones( obs_gts[oo]['obs']['zlvl'], dtype=np.bool )
    for vname in ['wspd','wdir','rh','t']:
      flag_nan *= np.isnan( obs_gts[oo]['obs']['vert'][vname]['data'] )

    # Removing useless vertical obs lvls
    for vname in obs_vlist:
      obs_gts[oo]['obs']['vert'][vname]['data'] \
        = obs_gts[oo]['obs']['vert'][vname]['data'][np.invert(flag_nan)]
      obs_gts[oo]['obs']['vert'][vname]['err'] \
        = obs_gts[oo]['obs']['vert'][vname]['err'][np.invert(flag_nan)]
      obs_gts[oo]['obs']['vert'][vname]['qc'] \
        = obs_gts[oo]['obs']['vert'][vname]['qc'][np.invert(flag_nan)]

    # Adjust total number of zlvls
    obs_gts[oo]['obs']['zlvl'] = np.sum( np.invert(flag_nan) )

    # Move onto next observation.

  # Finished comparing all obs. Returning.
  return       






'''
Function to write obs back into an obs gts file
'''
def write_obs_gts( ctrl_params, obs_gts ):

  # Interprete ctrl_params
  datestr, date, obs_vlist, era5_vlist, nSamp \
  = read_ctrl_params( ctrl_params )
  nObs = obs_gts['nObs']

  # Turn all nans in the obs_gts into -888888
  for oo in np.arange(nObs):
    # Processing surface obs
    for vname in ['slp', 'pw']:
      if obs_gts[oo]['obs']['sfc']['slp']['data'] == np.nan:
        obs_gts[oo]['obs']['sfc']['slp']['data'] = -888888.
#        obs_gts[oo]['obs']['sfc']['slp']['qc'] = -88
    # Processing vertical obs
    for vname in obs_vlist:
      nanflag = np.isnan( obs_gts[oo]['obs']['vert'][vname]['data'] )
      obs_gts[oo]['obs']['vert'][vname]['data'][nanflag] = -888888.
 #     obs_gts[oo]['obs']['vert'][vname]['qc'][nanflag] = -88



  # Preparing body string
  body_str = ''
  for oo in np.arange(nObs):
    # Info string
    info_str = ("%-12s %-19s %-40s %6d%12.3f" + " "*11 +"%12.3f" + " "*11 +"%12.3f"+" "*17+"%-40s\n")  % ( obs_gts[oo]['obs']['fm'], obs_gts[oo]['obs']['date'], obs_gts[oo]['obs']['name'], obs_gts[oo]['obs']['zlvl'], obs_gts[oo]['obs']['lat'], obs_gts[oo]['obs']['lon'], obs_gts[oo]['obs']['ele'], obs_gts[oo]['obs']['id'] )

    # srfc string
    srfc_str = "%12.3f%4d%7.2f%12.3f%4d%7.3f\n" % ( obs_gts[oo]['obs']['sfc']['slp']['data'], 
                                                          obs_gts[oo]['obs']['sfc']['slp']['qc'],
                                                          obs_gts[oo]['obs']['sfc']['slp']['err'],
                                                          obs_gts[oo]['obs']['sfc']['pw']['data'], 
                                                          obs_gts[oo]['obs']['sfc']['pw']['qc'],
                                                          obs_gts[oo]['obs']['sfc']['pw']['err'] )

    # zlvl string
    vert_str = ""
    if obs_gts[oo]['obs']['zlvl'] > 0:
        for kk in range( obs_gts[oo]['obs']['zlvl'] ):
            loc_str = ""
            for vv in ["pres","wspd","wdir"]:
                locloc_str = "%12.3f%4d%7.2f" % (obs_gts[oo]['obs']['vert'][vv]['data'][kk],
                                                    obs_gts[oo]['obs']['vert'][vv]['qc'][kk],
                                                    obs_gts[oo]['obs']['vert'][vv]['err'][kk])
                loc_str += locloc_str
            loc_str += " "*11
            for vv in ['z','t','dwpt']:
                locloc_str = "%12.3f%4d%7.2f" % (obs_gts[oo]['obs']['vert'][vv]['data'][kk],
                                                    obs_gts[oo]['obs']['vert'][vv]['qc'][kk],
                                                    obs_gts[oo]['obs']['vert'][vv]['err'][kk])
                loc_str += locloc_str
            loc_str += " "*11
            for vv in ['rh']:
                locloc_str = "%12.3f%4d%7.2f" % (obs_gts[oo]['obs']['vert'][vv]['data'][kk],
                                                    obs_gts[oo]['obs']['vert'][vv]['qc'][kk],
                                                    obs_gts[oo]['obs']['vert'][vv]['err'][kk])
                loc_str += locloc_str
            vert_str += loc_str
            vert_str += "\n"

    # Append strings tgt
    body_str += info_str + srfc_str + vert_str 

  # Prepare header string
  header_dict = {}
  header_dict['TOTAL'] = 0
  header_dict['MISS.'] = -888888.
  tstr_list = ['SYNOP','METAR','SHIP','BUOY','BOGUS','TEMP','AMDAR','AIREP','TAMDAR',\
               'PILOT','SATEM','SATOB','GPSPW','GPSZD','GPSRF','GPSEP','SSMT1','SSMT2',\
               'TOVS', 'QSCAT','PROFL','AIRSR','OTHER']

  for tstr in tstr_list:
      header_dict[tstr] = 0

  # Count number of obs
  for oo in np.arange(nObs):
    tstr = (obs_gts[oo]['obs']['fm'].split())[1]
    if tstr in header_dict:
        header_dict[tstr] += 1
    else:
        header_dict['OTHER'] +=1
  for tstr in tstr_list:
      header_dict['TOTAL'] += header_dict[tstr]

  # Print out headers
  head_str = ""
  head_str += "TOTAL =%7d, MISS. =-888888.,\n" % (header_dict['TOTAL'] )
  cnt = 0
  for tstr in tstr_list:
      head_str += "%-6s=%7d, " % (tstr, header_dict[tstr])
      cnt += 1
      if cnt == 6:
          head_str += "\n"
          cnt=0
  head_str += "\n"
  header_invariants = "PHIC  =  -2.70, XLONC = 107.50, TRUE1 =   0.00, TRUE2 = -20.00, XIM11 =   1.00, XJM11 =   1.00,\nbase_temp= 290.00, base_lapse=  50.00, PTOP  =  2000., base_pres=100000., base_tropo_pres= 20000., base_strat_temp=   215.,\nIXC   =    600, JXC   =   1300, IPROJ =      3, IDD   =      1, MAXNES=      1,\nNESTIX=    450, \nNESTJX=    560, \nNUMC  =      1, \nDIS   =   9.00, \nNESTI =      1, \nNESTJ =      1, \nINFO  = PLATFORM, DATE, NAME, LEVELS, LATITUDE, LONGITUDE, ELEVATION, ID.\nSRFC  = SLP, PW (DATA,QC,ERROR).\nEACH  = PRES, SPEED, DIR, HEIGHT, TEMP, DEW PT, HUMID (DATA,QC,ERROR)*LEVELS.\nINFO_FMT = (A12,1X,A19,1X,A40,1X,I6,3(F12.3,11X),6X,A40)\nSRFC_FMT = (F12.3,I4,F7.2,F12.3,I4,F7.3)\nEACH_FMT = (3(F12.3,I4,F7.2),11X,3(F12.3,I4,F7.2),11X,3(F12.3,I4,F7.2))\n#------------------------------------------------------------------------------#\n"

  head_str += header_invariants
 
  # Combine all strings
  output = head_str + body_str
  
  # Output QC'ed file
  outf = open( obs_gts['obs_out_fname'] ,'w')
  outf.write( output)
  outf.close()
 
  # Returning from writing out stuff
  return

