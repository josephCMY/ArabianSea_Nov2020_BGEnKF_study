# Section 2: Setting up WRF, WRF Preprocessing System (WPS) and WRF Data Assimilation system (WRFDA) on Expanse
Authors: Man-Yau (Joseph) Chan and Yunji (Jerry) Zhang

[Click here to return to the main navigation page](../../../README.md)

[Click here to see all Expanse-specific PSU-EnKF setup steps](Setting_Up_On_Expanse.md)

&nbsp;
## Overview

The PSU-EnKF uses Advanced Research Weather Research and Forecasting model (WRF-ARW, henceforth WRF) as the forecast model. 
To prepare the initial and boundary conditions to run WRF, the WPS is needed.
To process the observation files needed by the PSU-EnKF, the WRFDA's observation processing program (`obsproc.exe`) and related files are needed.

In this section, we will cover compiling each of these three components (WRF, WPS and WRFDA). Here are the hyperlinks to navigate to those steps:

2.1) [Instructions to compile WRF on Expanse](#21-compiling-wrf-on-expanse)

2.2) [Instructions to compile WPS on Expanse](#22-compiling-wps-on-expanse)

2.3) [Instructions to compile WRFDA on Expanse](#23-compiling-wrfda-on-expanse)

&nbsp;
&nbsp;
&nbsp;
&nbsp;
&nbsp;
&nbsp;
&nbsp;
&nbsp;
&nbsp;


## 2.1 Compiling WRF on Expanse

### Step 2.1.1 Obtain the WRF model code
Joseph's bare bones system comes with the source code for WRF version 3.8.1. Skip this step if you want to use WRF v3.8.1.

If you want to use some other WRF version, do the following:
1) Go to the [UCAR website](https://www2.mmm.ucar.edu/wrf/users/download/get_sources.html) using your web browser
2) Find the URL containing the source code of the desired WRF version.
3) On your terminal, issue `wget WRF_SOURCE_CODE_URL` to download the compressed source code.
4) De-compress the downloaded source code.


### Step 2.1.2 Set up the environment to compile WRF
Enter a compute node by issuing:
```
srun --partition=shared --pty --account=pen116 --nodes=1 --ntasks-per-node=8 --wait=0 --export=ALL -t 02:00:00 /bin/bash
```

Then, enter the directory containing the WRF source code using `cd`.


Activate the `intel` module set (which was created using instructions [here](Expanse_Mod_Setup.md)). This activation can be done by issuing:
```
module restore intel
```

### Step 2.1.3 Set up a local NetCDF4 library
Then, set up a local NetCDF4 library on your $HOME directory by issuing the following commands:
```
mkdir $HOME/lib
mkdir $HOME/include
cp $NETCDF_CHOME/include/* $HOME/include
cp $NETCDF_CHOME/lib/* $HOME/lib
cp $NETCDF_FORTRANHOME/include/* $HOME/include
cp $NETCDF_FORTRANHOME/lib/* $HOME/lib
export NETCDF=$HOME
```
The local NetCDF4 library is needed because Expanse has separate libraries for NetCDF compiled with C and FORTRAN, but WRF need both to compile (libnetcdf.a and libnetcdff.a). 


Some Intel compiler specific environmental variables also need to be defined to compile WRF. Just issue:
```
export I_MPI_CC=icc
export I_MPI_CXX=icpc
export I_MPI_F90=ifort
export I_MPI_F77=ifort
export I_MPI_FC=ifort
```

You will also need to set up some NetCDF-specific environmental variables to compile WRF. Just issue:
```
export WRFIO_NCD_LARGE_FILE_SUPPORT=1
export HDF5=$HDF5HOME
export PNETCDF=$PARALLEL_NETCDFHOME
export NETCDF_classic=1
```


### Step 2.1.4 Configure WRF and make Expanse-specific edits to the configuration
Assuming that this test case (`QuickStart_TropicalSquallCase`) lives inside the parent directory (absolute path `PARENT_DIR`), enter the WRF source code directory by issuing
```
cd $PARENT_DIR/QuickStart_TropicalSquallCase/fortran_src/WRFV3
```

Now prepare the WRF configuration by issuing:
```
./clean -a	# Wipe out any pre-existing WRF compilation in the directory
./configure
```
The terminal will then prompt you about which compilers to use and what type of nesting to use.
For Expanse, compiler option 66 works (67 might also work). 
To use Joseph's bare bones system, nesting option 1 (basic) is recommended. Joseph has never tested using options 2 (preset moves) and 3 (vortex following) for his bare bones system.


Some Expanse-specific adjustments need to be made to the configuration file. Just issue:
```
sed -i "s/-openmp/-qopenmp/g" configure.wrf
sed -i "s/xHost/march=core-avx2/g" configure.wrf
sed -i "s/-xCORE-AVX2/-I\/cm\/shared\/apps\/spack\/cpu\/opt\/spack\/linux-centos8-zen2\/intel-19.1.1.217\/libtirpc-1.2.6-vunl6q5hylerbcw7l6veqliruhyfjimz\/include\/tirpc/g" configure.wrf
sed -i "s/-lz/-lz -L\/cm\/shared\/apps\/spack\/cpu\/opt\/spack\/linux-centos8-zen2\/intel-19.1.1.217\/libtirpc-1.2.6-vunl6q5hylerbcw7l6veqliruhyfjimz\/lib -ltirpc/g" configure.wrf
```

### Step 2.1.5 Compile WRF
To compile WRF, just issue:
```
./compile em_real >& log.compile & tail -f log.compile &
```

The compilation will take roughly an hour to complete. When the compilation has successfully completed, a message like this will appear on your terminal
```
==========================================================================
build started:   Wed Mar 16 09:49:15 PDT 2022
build completed: Wed Mar 16 10:40:30 PDT 2022
 
--->                  Executables successfully built                  <---
 
-rwxr-xr-x 1 chanmy pen116 46697024 Mar 16 10:40 main/ndown.exe
-rwxr-xr-x 1 chanmy pen116 46688608 Mar 16 10:40 main/real.exe
-rwxr-xr-x 1 chanmy pen116 45974240 Mar 16 10:40 main/tc.exe
-rwxr-xr-x 1 chanmy pen116 53066920 Mar 16 10:39 main/wrf.exe
 
==========================================================================
```
Just hit CTRL-C on your keyboard to regain access of your terminal window.

Then, hit CTRL-D to exit the compute node.


&nbsp;
&nbsp;
&nbsp;
&nbsp;
&nbsp;
&nbsp;
&nbsp;
&nbsp;
&nbsp;







## 2.2 Compiling WPS on Expanse

### Step 2.2.1 Obtain the WPS code
Joseph's bare bones system comes with the source code for WPS version 3.8.1. Skip this step if you want to use WPS v3.8.1.

If you want to use some other WPS version, do the following:
1) Go to the [UCAR website](https://www2.mmm.ucar.edu/wrf/users/download/get_sources.html) using your web browser
2) Find the URL containing the source code of the desired WPS version.
3) On your terminal, issue `wget WPS_SOURCE_CODE_URL` to download the compressed source code.
4) De-compress the downloaded source code.

**IMPORTANT NOTE**: The WPS source code directory needs to be in the same parent directory as the WRF source code directory.


### Step 2.2.2 Enter compute node and activate intel module environment
Enter a compute node by issuing:
```
srun --partition=shared --pty --account=pen116 --nodes=1 --ntasks-per-node=8 --wait=0 --export=ALL -t 02:00:00 /bin/bash
```

Activate the `intel` module set (which was created using instructions [here](Expanse_Mod_Setup.md)). This activation can be done by issuing:
```
module restore intel
```


### Step 2.2.3 Download and compile ZLIB, PNG and JASPER packages

At this point, you should have already finished compiling WRF. If you followed the [instructions for compiling WRF](#21-compiling-wrf-on-expanse), you should have already set up a local NetCDF4 library on your $HOME directory. If the local NetCDF4 library is not set up, please checkout out the [instructions to do so](#21.3-set-up-a-local-netcdf4-library).


Go to your $HOME directory and download the ZLIB, PNG and JASPER packages by issuing:
```
cd ~
wget https://www2.mmm.ucar.edu/wrf/OnLineTutorial/compile_tutorial/tar_files/libpng-1.2.50.tar.gz
wget https://www2.mmm.ucar.edu/wrf/OnLineTutorial/compile_tutorial/tar_files/zlib-1.2.7.tar.gz
wget https://www.ece.uvic.ca/~frodo/jasper/software/jasper-1.900.0.zip
```

Then, untar the resulting `tar.gz` files by issuing:
```
unzip jasper-1.900.0.zip
tar -xzvf libpng-1.2.50.tar.gz
tar -xzvf zlib-1.2.7.tar.gz
```

Compile the ZLIB package by issuing:
```
cd ~/zlib-1.2.7
make clean
./configure --prefix=~/
make 
make install
```

Compile the LIBPNG package by issuing:
```
cd ~/libpng-1.2.50
make clean
./configure --prefix=$HOME
make
make install
```

Compile the JASPER package by issuing:
``` 
cd ~/jasper-1.900.0
export CC=icc
export FC=ifort
make clean
./configure --prefix=$HOME
make
make install
```



### Step 2.2.4 Set up the environment to compile WPS

At this point, you should have already finished compiling WRF. If you followed the [instructions for compiling WRF](#21-compiling-wrf-on-expanse), you should have already set up a local NetCDF4 library on your $HOME directory. If the local NetCDF4 library is not set up, please checkout out the [instructions to do so](#21.3-set-up-a-local-netcdf4-library).

Set up an environmental variable to indicate where the local NetCDF file is:
```
export NETCDF=$HOME
```

Some Intel compiler specific environmental variables also need to be defined to compile WPS. Just issue:
```
export I_MPI_CC=icc
export I_MPI_CXX=icpc
export I_MPI_F90=ifort
export I_MPI_F77=ifort
export I_MPI_FC=ifort
```

You will also need to set up some NetCDF-specific environmental variables to compile WPS. Just issue:
```
export WRFIO_NCD_LARGE_FILE_SUPPORT=1
export HDF5=$HDF5HOME
export PNETCDF=$PARALLEL_NETCDFHOME
export NETCDF_classic=1
```

Now specify some environmental variables for the Jasper library
```
export JASPERLIB=$HOME/lib
export JASPERINC=$HOME/include
```

Now add the `$HOME/lib` and `$HOME/bin` to your environment list of paths. To do so, add the following lines to your `$HOME/.bashrc` file:
```
export PATH=$PATH:$HOME/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/lib
```



### Step 2.2.3: Configure WPS and make Expanse-specific edits to the configuration
Enter the directory containing the WPS source code via `cd PATH_TO_WPS_SOURCE_CODE`.

Now prepare the WPS configuration by issuing:
```
./clean -a	# Wiping out any pre-existing WPS compilation
./configure
```
The terminal will then prompt you about which compilers to use.
For Expanse, compiler option 19 works.


### Step 2.2.4: Compile WPS
To compile WPS, just issue:
```
./compile >& log.compile & tail -f log.compile 
```

The compilation will take < 10 minutes to complete. When the compilation is complete, you will see this message on your terminal:
```
if [ -h int2nc.exe ] ; then \
	/bin/rm -f int2nc.exe ; \
fi ; \
if [ -h ../int2nc.exe ] ; then \
	/bin/rm -f ../int2nc.exe ; \
fi ; \
if [ -e src/int2nc.exe ] ; then \
	ln -sf src/int2nc.exe . ; \
fi
```
Just hit CTRL-C on your keyboard to regain access of your terminal window.

If the compilation is successful, the following executables should exist within the WPS directory:
```
geogrid/src/geogrid.exe
metgrid/src/metgrid.exe
ungrib/src/ungrib.exe
```

If these three executables exist, you have succesfully compiled WPS! =)

Then, hit CTRL-D to exit the compute node.



&nbsp;
&nbsp;
&nbsp;
&nbsp;
&nbsp;
&nbsp;
&nbsp;
&nbsp;
&nbsp;



## Step 2.3 Compiling WRFDA on Expanse (UNDER DEVELOPMENT)
### Step 2.3.1 Obtain the WRFDA code
Joseph's bare bones system comes with the source code for WRFDA version 3.8.1. Skip this step if you want to use WRFDA v3.8.1.

If you want to use some other WPS version, do the following:
1) Go to the [UCAR website](https://www2.mmm.ucar.edu/wrf/users/download/get_sources.html) using your web browser
2) Find the URL containing the source code of the desired WRFDA version.
3) On your terminal, issue `wget WRFDA_SOURCE_CODE_URL` to download the compressed source code.
4) De-compress the downloaded source code.

**IMPORTANT NOTE**: The WRFDA source code directory needs to be in the same parent directory as the WRF source code directory.


### Step 2.3.2 Enter compute node and activate intel module environment
Enter a compute node by issuing:
```
srun --partition=shared --pty --account=pen116 --nodes=1 --ntasks-per-node=8 --wait=0 --export=ALL -t 02:00:00 /bin/bash
```

Activate the `intel` module set (which was created using instructions [here](Expanse_Mod_Setup.md)). This activation can be done by issuing:
```
module restore intel
```


### Step 2.3.4 Set up the environment to compile WPS
At this point, you should have already finished compiling WRF. If you followed the [instructions for compiling WRF](#21-compiling-wrf-on-expanse), you should have already set up a local NetCDF4 library on your $HOME directory. If the local NetCDF4 library is not set up, please checkout out the [instructions to do so](#21.3-set-up-a-local-netcdf4-library).

Set up an environmental variable to indicate where the local NetCDF file is:
```
export NETCDF=$HOME
```


Some Intel compiler specific environmental variables also need to be defined to compile WRFDA. Just issue:
```
export I_MPI_CC=icc
export I_MPI_CXX=icpc
export I_MPI_F90=ifort
export I_MPI_F77=ifort
export I_MPI_FC=ifort
```

You will also need to set up some NetCDF-specific environmental variables to compile WPS. Just issue:
```
export WRFIO_NCD_LARGE_FILE_SUPPORT=1
export HDF5=$HDF5HOME
export PNETCDF=$PARALLEL_NETCDFHOME
export NETCDF_classic=1
```


### Step 2.3.5 Configure and compile WRFDA 
Now prepare the WRF configuration by issuing:
```
./configure wrfda
```
Use compiler option 66 (_INTEL (ifort/icc): HSW/BDW_)

Once the configuration is complete, compile WRFDA by running:
```
/compile all_wrfvar >& log.compile &
tail -f log.compile
```

At this point, the guide for compiling WRFDA is incomplete. Only the most important program for the PSU-EnKF (`obsproc.exe`) is compiled. Of all the WRFDA programs, only `obsproc.exe` alone is needed for running the PSU-EnKF.

Frankly, `obsproc.exe` can be replaced with a Python code to process observations into a special file that the PSU-EnKF `enkf.mpi` can read. 


&nbsp; 

[Click here to return to the main navigation page](../../../README.md)

[Click here to see all Expanse-specific PSU-EnKF setup steps](Setting_Up_On_Expanse.md)
