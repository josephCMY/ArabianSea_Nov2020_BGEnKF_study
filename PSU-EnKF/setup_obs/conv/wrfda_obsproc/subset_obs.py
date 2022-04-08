''' 
Script to preprocess NCAR RDS GTS files for accelerated obsproc.exe

Procedure
    1) For each 6 hour interval, read in obs that are within
       the reanalysis zone. Aggregate said obs.       
    2) Separate obs based on time windows and output.
'''

''' Basic libaries '''
import numpy as np
from datetime import datetime, timedelta
import sys
from netCDF4 import Dataset as ncopen

# User inputs
date_st = datetime.strptime(sys.argv[1], "%Y%m%d%H%M")
date_ed = datetime.strptime(sys.argv[2], "%Y%m%d%H%M")

# NCAR obs interval
ncar_t_int = timedelta( minutes=360 )

# EnKF cycle period
enkf_t_int = timedelta( minutes=60 )

# Set up list of dates
date=date_st + timedelta( minutes=0)
datelist = []
while date <= date_ed:
  datelist.append(date)
  date += ncar_t_int

# Deterine lat lon coordinate limits from geo_em.d01.nc
geo_f = ncopen('geo_em.d01.nc','r')
lat_lims = np.array( geo_f.variables['XLAT_V'][:])
lat_lims = [ lat_lims.min()+0.25, lat_lims.max()-0.25 ]
lon_lims = np.array( geo_f.variables['XLONG_U'][:])
lon_lims = [ lon_lims.min()+0.25, lon_lims.max()-0.25 ]
geo_f.close()

# Setting up object to hold observations
obs_data = {}
obs_data['txt'] = []
obs_data['date'] = []

# Searching for relevant obs across dates
print("Reading in global obs")
for date in datelist:

  print( "Handling obs on %s" % date.strftime( "%Y%m%d%H"))

  # Iterate over OBS and SURFACE_OBS files
  for prefix in ['raw_conv/OBS:','raw_conv/SURFACE_OBS:']:

    # Load all lines
    fname = "%s%s" % (prefix, date.strftime( "%Y%m%d%H"))
    print(fname)
    f = open( fname, 'r', encoding="ISO-8859-1")
    raw_lines = [line for line in f ]
    nL = len(raw_lines)

    # Seek out line numbers where individual obs entries start
    obsflags = np.zeros( nL, dtype=bool )
    for ll in np.arange(  nL ):
      if raw_lines[ll].find( 'FM-') != -1:
        obsflags[ll] = True

    # Determine line numbers
    fm_ind = np.arange(nL)[obsflags]
    nObs = len(fm_ind)

    # Go to each obs
    for oo in np.arange(nObs):

      # Starting line number
      l0 = fm_ind[oo]
      # Determine ending line number
      if oo < nObs -1:
        l1 = fm_ind[oo+1]
      else:
        l1 = nL

      # Check if position is relevant
      header = raw_lines[l0]
      header = header.split()
      lat = float(header[0])
      lon = float(header[1])

      # Skip obs if outside of domain
      if lon < lon_lims[0] or lon > lon_lims[1] \
         or lat < lat_lims[0] or lat > lat_lims[1]:
        continue        

      # Store observation into the list
      # Tag obs time
      header = raw_lines[l0]
      obs_time = header[326:(326+14)]
      obs_time = datetime.strptime( obs_time, "%Y%m%d%H%M%S")
      obs_lines = raw_lines[l0:l1]
      obs_data['txt'].append( obs_lines )
      obs_data['date'].append( obs_time)
      
obs_data['date'] = np.array( obs_data['date'] )

print( "Proceeding to generate subsetted obs files")

''' Now outputting observations based on EnKF cycling times '''

# Identify EnKF times
out_datelist = []
date = date_st + timedelta(minutes=0)
while date <= date_ed:
  out_datelist.append(date)
  date += enkf_t_int

# Selecting relevant obs and outputting
for date in out_datelist:

  # Determining bounding dates
  date0 = date - enkf_t_int/2
  date1 = date + enkf_t_int/2

  # Determine indices of relevant obs
  dateflags = (date0 <= obs_data['date'] )* (obs_data['date'] < date1)
  dateflags = np.array( dateflags, dtype=np.bool)
  obs_inds = np.arange( len( dateflags) )[dateflags]

  # Generate file
  outfname = ( "raw_conv/subsetted_obs:%s" \
               % date.strftime("%Y%m%d%H%M") )
  outf = open( outfname,'w')

  # Write to file
  for oi in obs_inds:
    write_lines = obs_data['txt'][oi]
    for line in write_lines:
      outf.write(line)

  # Close file
  outf.close()
  print( "Generated subsetted file for %s"  % date.strftime("%Y%m%d%H%M") )



  






