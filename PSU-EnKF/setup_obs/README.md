# Downloading and processing of observations

Observation sources:
1) In-situ observations from NCAR RDA
2) Infrared observations from Meteosat-8 (IODC)
3) Infrared observations from Himawari-8


## Procedure for in-situ observations
1) Download observations from listed sources
2) Process in-situ and AMV obs via obsproc into relevant LITTLE\_R .3DVAR files.
3) Compare .3DVAR files against ERA5. Reject all obs that are far from ERA5 (3x obs sigma)

Note: We are no longer using AMV from CIMSS. Just sticking to NCAR GTS AMV.


## Procedure for IR observations

