'''
Script to run multiscale validation on ensemble
'''
import numpy as np
import math
import datetime
import sys
import pickle
from netCDF4 import Dataset as ncopen
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

''' Script control parameters '''
chlist = ['ch08','ch10','ch14'] #['ch08','ch09','ch10','ch13','ch14','ch15']
new_chlist = { 'ch08': 'ahi_himawari8_ch008', 
               'ch10': 'ahi_himawari8_ch010',
               'ch14': 'ahi_himawari8_ch014' }
ref_fname = 'him8_xa_001.bin'
new_fname = 'ap_crtm.nc'
nY = 449
nX = 559

data = {}
''' Load the reference file '''
data['ref'] = {}
arr = np.fromfile( ref_fname, dtype = '>f4' )
arr = arr.reshape( len(chlist) +2, nY, nX )
data['ref']['lon'] = arr[0,:,:]
data['ref']['lat'] = arr[1,:,:]
for i in range(len(chlist)):
  ch = chlist[i]
  data['ref'][ch] = arr[i+2,:,:]


''' Load the new file '''
data['new'] = {}
f = ncopen( new_fname,'r') 
for ch in chlist:
  data['new'][ch] = f.variables[ new_chlist[ch] ][:]

data['new']['lon'] = f.variables[ 'longitude' ][:]
data['new']['lat'] = f.variables[ 'latitude' ][:]

f.close()


''' Compare fields '''
vlist = ['lon','lat']
for ch in chlist:
  vlist.append( ch )
for vv in vlist:
  print( vv,'max abs diff', (np.fabs(data['new'][vv] - data['ref'][vv])).max() )
  print( vv,'avg diff', np.mean(data['new'][vv] - data['ref'][vv]))

print( data['new']['ch14'][100,100], data['ref']['ch14'][100,100] )
quit()

''' Plot out the ch14 fields '''
fig = plt.figure( figsize=(6,6))
axs = []
for i in range(3):
    axs.append( fig.add_subplot(2,2,i+1) )
ch = 'ch14'
crange=np.linspace(200,280,11)

# Ref field
axs[0].contourf( data['new']['lon'], data['new']['lat'], data['ref'][ch],
                 crange, cmap = 'RdBu_r', extend='min')
axs[0].set_title('ref')
axs[1].contourf( data['new']['lon'], data['new']['lat'], data['new'][ch],
                 crange, cmap = 'RdBu_r', extend='min')
axs[1].set_title('new')
axs[2].contourf( data['new']['lon'], data['new']['lat'], 
                 data['new'][ch] - data['ref'][ch],
                 np.linspace(-5,5,11), cmap = 'RdBu_r', extend='min')
axs[2].set_title('new - old')
plt.savefig('sanity_check.png' )
