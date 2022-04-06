''' Script to convert SEVIRI native file into netCDF '''


# import modules
from satpy import Scene
import sys
import numpy as np
import datetime
from netCDF4 import Dataset as ncopen
import pickle


# ============================================================================
# USER INPUTS
# ----------------------------------------------------------------------------

# File name of SEVIRI native
fname = sys.argv[1]

# Geographical limits bounding the desired area
latmin = -20.
latmax =  20.
lonmin =  50.
lonmax = 100.

# List of desired channel names (names must follow native file's convention)
chlist = ['WV_062','IR_108']



# ===========================================================================
# Loading native file and subsetting to specified region
# ---------------------------------------------------------------------------

# Open SEVIRI file
scn = Scene( reader='seviri_l1b_native', filenames=[fname])

# Load channel data from file
scn.load(chlist)

# Load longitude and latitude
ch = chlist[0]
lon2d, lat2d = scn[ch].attrs['area'].get_lonlats()
shp2d = lat2d.shape
halfway_inds = np.array( np.array(shp2d, dtype=int)/2, dtype=int )

# Set up domain masks
latvec = lat2d[ :, halfway_inds[1] ]
latmask = ( latvec >= latmin ) * (latvec <= latmax )
lonvec = lon2d[ halfway_inds[0], : ]
lonmask = ( lonvec >= lonmin ) * ( lonvec <= lonmax )

# Subset data to desired area
sub_data = {}
sub_data['lon'] = np.array( lon2d[latmask][:,lonmask] )
sub_data['lat'] = np.array( lat2d[latmask][:,lonmask] )
for ch in chlist:
  sub_data[ch] = np.array( (scn[ch].values)[latmask][:,lonmask] )

# Load time period in which the sensing was done
date_st = scn[ch].attrs['start_time']
date_ed = scn[ch].attrs['end_time']



# =========================================================================
# Output subsetted data as a netcdf file
# -------------------------------------------------------------------------

# Devise output ncfile name and open ncfile
out_fname = ( 'extracted_seviri_ncfiles/seviri_%sUTC.nc' 
              % date_ed.strftime("%Y-%m-%d_%H%M") )
f = ncopen( out_fname, 'w')

# Construct dimensions
shp2d = sub_data[ch].shape
lat_dim = f.createDimension('lat', shp2d[0] )
lon_dim = f.createDimension('lon', shp2d[1] )

# Construct attributes
f.description = "Extracted SEVIRI BT data"
f.sensor_start_time = date_st.strftime("%Y-%m-%d_%H:%M:%S")
f.sensor_end_time   = date_ed.strftime("%Y-%m-%d_%H:%M:%S")

# Construct and store variables into the ncfile
lat = f.createVariable( "latitude" , np.float32, ("lat","lon",))
lat.units = 'deg_north'
lat[:] = sub_data['lat']*1.

lon = f.createVariable( "longitude", np.float32, ("lat","lon",))
lon.units = 'deg_south'
lon[:] = sub_data['lon']

for ch in chlist:
  bt = f.createVariable( ch, np.float32, ("lat","lon",))
  bt.units = 'Kelvins'
  bt[:] = sub_data[ch]

# Close ncfile and flush to harddrive
f.close()





