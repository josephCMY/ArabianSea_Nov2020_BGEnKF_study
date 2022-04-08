''' 
Script to download Himawari-8 observations from P-Tree 
Got original script from Yinghui Lu and modified accordingly

'''

import datetime
import os
import urllib.request
import sys


# Function to retrieve radiance files from the P-Tree server
# timeIn is the observation time as a datetime obj. ch is the channel number as an integer
def retrieve_him8_radiance(time):

  localPath = 'PATH_TO_DIRECTORY_TO_STORE_HIM8_DATA'
  
  # USER_NAME should be email address used in registration, and substuting '@' with '_'
  # Path to retrieve radiance data from
  remotePath = 'ftp://USER_NAME:SP+wari8@ftp.ptree.jaxa.jp/' 

  # I guess the following three lines determines which type of file (full-disk vs. Japan, for example)
  npixel= 2401
  nline = 2401
  # Special note: R21 retrieves all full-disk channels. Cant find files for individual channels on the server
  bandstr = 'R21'

  # Create radiance file name
  filename = 'NC_H08_' + time.strftime('%Y%m%d_%H%M_') + '%s_FLDK.%05d_%05d.nc' % (bandstr, npixel, nline) 

  # Subpath of the JAXA server
  subPath = 'jma/netcdf/'+time.strftime('%Y%m/%d/')

  # Download data into local path
  urllib.request.urlretrieve( (remotePath + subPath + filename),
                              (localPath + filename))
  return


if __name__ == '__main__':

  retrieve_him8_radiance(datetime.datetime.strptime(sys.argv[1],'%Y%m%d%H%M%S'))
