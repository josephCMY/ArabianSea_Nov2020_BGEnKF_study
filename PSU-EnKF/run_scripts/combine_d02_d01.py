''' 
Python script to combined wrfout_d02 with wrfout_d01, so that CRTM can be performed.

Variables to combine: P, PB, PH, PHB, T, QVAPOR, QCLOUD, QRAIN, QSNOW, QICE, QGRAUP

Note that the script is only designed to combined unstaggered variables!
Note also that wrfout_d01 will be overwritten by this script


The code will figure out where the 
'''

import numpy as np
from netCDF4 import Dataset as ncopen
import sys

''' Control parameters '''
d01_dx = 9.0
d01_dy = 9.0
d02_dx = 3.0
d02_dy = 3.0
#varname_list =['P','XLAT','XLONG']
varname_list = ['P','PB','PH','PHB','T','QVAPOR','QCLOUD','QRAIN','QICE','QSNOW', \
                'QGRAUP', 'PSFC','HGT','XLAT','XLONG' ]
grid_ratio = { 'dx': int( d01_dx / d02_dx), 'dy': int( d01_dy / d02_dy) }

''' User inputs '''
fname_d01 = sys.argv[1]         # wrfout_d01 filename
fname_d02 = sys.argv[2]         # wrfout_d02 filename


''' Open files '''
flist  = { 'd01': ncopen( fname_d01, 'r+'),
           'd02': ncopen( fname_d02, 'r') }


''' Read in all of the desired variables '''
data = {}
for dom in ['d01','d02']:
  data[dom] = {}
  for vname in varname_list:
    data[dom][vname] = np.array( flist[dom].variables[vname] ) * 1.0


''' Now take 2D spatial average for the grid points in domain 2 '''
dom = 'd02'
for vname in varname_list:
    
  # Setting up array to help with broadcasting the averaging process
  # Going for a [ dy_ratio, dx_ratio, ...., ydim/ dy_ratio, xdim/dx_ratio] array
  shp = data[dom][vname].shape
  shp_tmp = np.array(shp )
  shp_tmp[-1] /= grid_ratio['dx']
  shp_tmp[-2] /= grid_ratio['dy']
  shp = [grid_ratio['dy'], grid_ratio['dx']]
  for ii in shp_tmp:
      shp.append( int(ii) )
  tmp = np.zeros( shp, dtype='float' )
  
  # Filling in temporary array
  for ii in np.arange(shp[0]):
    for jj in np.arange(shp[1]):
      # Dealing with 2D and 3D variables
      if len(shp) == 5:
        tmp[ii,jj,:,:,:] = data[dom][vname][:,ii::grid_ratio['dy'], jj::grid_ratio['dx'] ]
#        print( vname + ' '+str(len(shp)) )
      elif len(shp) == 6:
        tmp[ii,jj,:,:,:,:] = data[dom][vname][:,:,ii::grid_ratio['dy'], jj::grid_ratio['dx'] ]
#        print( vname + ' '+str(len(shp)) )

  
  # Take the average
  data[dom][vname] = np.mean( np.mean( tmp, axis=0 ), axis=0)
#  print(data[dom][vname].shape)


''' Now determine where the d02 variables fit into d01 '''
dist = (  np.power( data['d01']['XLAT'][0] - data['d02']['XLAT'][0,0,0],2) 
        + np.power( data['d01']['XLONG'][0] - data['d02']['XLONG'][0,0,0],2) )
mindist = dist.min()
if mindist > 1e-3:
  print( 'problematic mindist: %e' %mindist)
  print('exiting')
  quit()
ind_st = np.where( dist == mindist )
ind_st = [ind_st[0][0], ind_st[1][0]]

''' Now perform the overwriting '''
for vname in varname_list:
  arr = data['d02'][vname]
  shp = arr.shape
  nx = int(shp[-1])
  ny = int(shp[-2])
  print(vname, shp, nx, ny, ind_st)
  if len(shp) == 3:
    flist['d01'].variables[vname][:,ind_st[0]:(ind_st[0]+ny), ind_st[1]:(ind_st[1]+nx)]\
    = data['d02'][vname][:,:,:]

  elif len(shp) == 4:
    flist['d01'].variables[vname][:,:,ind_st[0]:(ind_st[0]+ny), ind_st[1]:(ind_st[1]+nx)]  \
    = data['d02'][vname][:,:,:,:]
    tmp2 = flist['d01'].variables[vname][0,:,ind_st[0] +50,ind_st[1]+50 ] 
#    print( vname, tmp1, tmp2)

  
# Close
for dom in ['d01','d02']:
    flist[dom].close()
