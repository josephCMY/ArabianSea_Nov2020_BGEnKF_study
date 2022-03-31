#!/bin/bash

# Master script to prepare conventional observations for DA
# System is controlled through config file.


## Part 1: request raw conventional files (NCAR-RDA GTS LITTLE_R and CIMSS AMV)
#./run_download_raw_conv.sh >& log.request

# Part 2: run NCAR-RDA GTS obs through WRFDA obsproc.exe for obs errors
cd wrfda_obsproc
./run_obsproc.sh >& log.obsproc
cd ..

# Part 3: check NCAR-RDA GTS against ERA5.
cd quality_control
./run_quality_control.sh >& log.quality_control
cd ..




