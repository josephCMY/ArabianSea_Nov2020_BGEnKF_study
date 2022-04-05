# Procedure to run the WPS and real.exe process for nature run 

1) Adjusted the `setup_namelist_wps.sh` and `setup_namelist_real.sh` to be consistent with new experiment domain
2) Manually ran `setup_namelist_wps.sh` to construct appropriate `namelist.wps`.
3) Created a directory `make_geo_em` containing symlinks to `namelist.wps`, `geogrid.exe` and WPS' `geogrid` directory.
4) In the `make_geo_em` directory, manually ran geogrid with `srun -n 16 geogrid.exe`.
5) Called `./link_grib.csh raw_era5/*grib` to symlink ERA5-related GRIB files
6) Called `srun -n 16 ungrib.exe` to ungrib the ERA5 GRIB files. 

