#!/bin/bash
# Script to retrieve conventional data

. config

# Retrieve conventioanl data from NCAR Research Data Archive
date_st=201708010000
date_ed=201709010000
date=$date_st
while [[ $date -le $date_ed ]]; do
  ccyy=`echo $date | cut -c1-4`
  mmddhh=`echo $date | cut -c5-10`
  ./download_ds351.0.csh $ccyy $mmddhh raw_conv
  ./download_ds461.0.csh $ccyy $mmddhh raw_conv
  echo Retrieved data for $date
  date=`advance_time $date 360`
done


# CIMSS TC AMV not used here since NCAR RDA has AMV data.
## Retrieve AMV from CIMSS TC group archive
#AMV_zones="Australia EuropeAfrica Indian NWPacific"
#date=$date_st
#while [[ $date -le $date_ed ]]; do
#
#  # Interpreting date
#  ccyymmdd=`echo $date | cut -c1-8`
#  hh=`echo $date | cut -c9-10`
#
#  # Retrieve AMV for each CIMSS zone
#  for zone in $AMV_zones; do
#    src_html=http://tropic.ssec.wisc.edu/archive/data/$zone/$ccyymmdd/AllWindsQIText/$ccyymmdd.$hh.$zone.AllWindsQIText
#    wget -O raw_conv/CIMSS_AMV_"$zone"_$ccyymmdd$hh $src_html
#  done
#
#  # Increment time
#  echo Retrieved CIMSS AMV for $date
#  date=`advance_time $date 180`
#done
