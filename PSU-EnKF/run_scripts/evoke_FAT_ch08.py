'''
Script to evoke FAT on an ensemble member using OLR observations.
'''

import numpy as np
from mpl_toolkits.basemap import Basemap
import datetime 
import sys
from netCDF4 import Dataset as ncopen
import glob
import wrf
import funclib_olr_validation as flib
from funclib_FAT import multismooth_FAT as msfat
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import pickle
from scipy.interpolate import RectBivariateSpline
import matplotlib.dates as mdates
myFmt = mdates.DateFormatter('%b-%d\n%HUTC')


''' User controls '''
# Date to make plots (user-specified)
date = datetime.datetime.strptime( sys.argv[1], '%Y%m%d%H%M')
varlist = [ 'OLR','U','V','T','Q','U10','V10','T2','Q2','QVAPOR', \
            'QRAIN','QCLOUD','QICE','QSNOW','QGRAUP','W' ]

wrf_rmsd = []
era5_rmsd = []

''' Loading WRF ensemble OLR-based BT '''
print( 'Processing RMSD and bias for %s' % date.strftime("%d %b, %H UTC"))
wrf_data = {}


##### Loading WRF OLR-BT #####
fname = 'wrfinput_d01_prior'
wrffile_prior = ncopen( fname, 'r')
wrf_data['xf olr'] = wrffile_prior.variables['OLR'][:]


# Load coordinates
wrf_data['lat'] = wrf.getvar(wrffile_prior,'XLAT' ,meta=False) 
wrf_data['lon'] = wrf.getvar(wrffile_prior,'XLONG',meta=False) 



##### Performing 1deg x 1deg averaing #####
avg_data0 = flib.avg_1deg( wrf_data['lon'], wrf_data['lat'], wrf_data['xf olr'] )
smth_data = {}
smth_data['xf olr'] = avg_data0['out'] * 1.
smth_data['lon'] = avg_data0['lon'] * 1.
smth_data['lat'] = avg_data0['lat'] * 1.
smth_data['xf olr-bt'] = np.power( smth_data['xf olr']/(5.67e-8), 0.25 )


''' Load the CERES data '''
ceres_data =  flib.load_ceres_OLR([date])
ceres_data['lon'], ceres_data['lat'] = np.meshgrid( ceres_data['lon'], ceres_data['lat'] )
ceres_data['bt'] = np.power( ceres_data['olr'] / (5.67e-8), 0.25 )[0]
# Subset CERES to WRF domain
lon_flag = (ceres_data['lon'][0,:] >= smth_data['lon'].min() ) * (ceres_data['lon'][0,:] <= smth_data['lon'].max() )
lat_flag = (ceres_data['lat'][:,0] >= smth_data['lat'].min() ) * (ceres_data['lat'][:,0] <= smth_data['lat'].max() )
for key in ['lon','lat','bt'] :
  ceres_data[key] = ceres_data[key][lat_flag,:]
  ceres_data[key] = ceres_data[key][:,lon_flag]
# Compute CERES obs error using Monte Carlo
olr_mean = (ceres_data['bt']**4) * (5.67e-8)
mc_shp = [100, olr_mean.shape[0], olr_mean.shape[1] ]
olr_samples = np.random.normal( loc=0, scale=1., size= mc_shp)
olr_samples -= np.mean( olr_samples, axis=0)
olr_samples /= np.std( olr_samples, ddof=1, axis=0)
olr_samples = olr_samples * 2.81 + olr_mean
olr_samples = np.power( olr_samples / (5.67e-8), 0.25 )
ceres_data['bt_err'] = np.std( olr_samples, ddof=1, axis=0)




##### Perform multsmooth FAT ######

# Perturb the observation field
yo2d_pert = np.random.normal( size = ceres_data['bt'].shape )
yo2d_pert -= np.mean( yo2d_pert)
yo2d_pert *= ceres_data['bt_err']
yo2d = ceres_data['bt'] + yo2d_pert

# Normalize the observation and background fields
yo2d /= ceres_data['bt_err']
yf2d = smth_data['xf olr-bt'][0] / ceres_data['bt_err'] 
shp = yf2d.shape
print( shp )
quit()

# Evoke msfat code (shear-scale set to 0.167)
ua2d, va2d = msfat( yo2d, yf2d, yo2d*0., yo2d*0., 0.167, maxiter=100,
                    maxsmth=4 )

# Map the msfat solution from geo degree space to physical space
# The grid spacing is 1 deg
# Rescaling meridional displacement
va2d *= (1. * 3.1415 / 180.) * 111.
# Rescaling the zonal displacement
ua2d *=  ( (1. * 3.1415 / 180.) * 111. 
           * np.cos( ceres_data['lat'] * 3.1415/180.) )

# Run interpolation from OLR location to the wrf locations
ufunc = RectBivariateSpline( ceres_data['lat'][:,0],
                             ceres_data['lon'][0,:], ua2d )
vfunc = RectBivariateSpline( ceres_data['lat'][:,0],
                             ceres_data['lon'][0,:], va2d )
ua2d_fine = ufunc( wrf_data['lat'], wrf_data['lon'], grid=False )
va2d_fine = vfunc( wrf_data['lat'], wrf_data['lon'], grid=False )

# Ensure that the out-of-site stuff is set to zero
mask = ( (wrf_data['lat'] >= smth_data['lat'].min())
        *(wrf_data['lat'] <= smth_data['lat'].max()) )
mask*= ( (wrf_data['lon'] >= smth_data['lon'].min())
        *(wrf_data['lon'] <= smth_data['lon'].max()) )
mask = np.invert( mask )
ua2d_fine[mask] = 0.
va2d_fine[mask] = 0.



######## Apply FAT field solution ###### 
fname = 'wrfinput_d01_fatted'
wrffile_fat = ncopen( fname,'a')
xmesh0, ymesh0 = np.meshgrid( np.arange(shp[1]), np.arange(shp[0]) )
xmesh1 = xmesh0 + ua2d_fine
ymesh1 = ymesh0 + va2d_fine
for vname in ['OLR']: #varlist:
  fld_prior = wrffile_prior.variables[vname][0]
  fld_fat   = wrffile_fat
  fshp = fld.shape
  if len( fshp ) == 2:
    fld_ifunc = RectBivariateSpline( np.arange(shp[0]), 
                                     np.arange(shp[1]), fld )
    fld_fatted = fld_ifunc( ymesh1, xmesh1, grid=False)
    fld




