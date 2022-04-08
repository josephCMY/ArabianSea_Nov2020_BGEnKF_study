# Section 4: Compiling the PSU-EnKF EnSRF algorithm code on Expanse
Author: Man-Yau (Joseph) Chan

[Click here to return to the main navigation page](../../../README.md)

[Click here to see all Expanse-specific PSU-EnKF setup steps](Setting_Up_On_Expanse.md)

&nbsp;
## Overview

In this section, we will cover compiling the PSU-EnKF Fortran EnSRF code on Expanse. 

### Step 4.1 Navigate to the EnSRF code directory
The EnSRF Fortran code lives in the `PSU-EnKF_w_satellites/fortran_src/EnSRF/src` directory. Just `cd` into this directory.


### Step 4.2 Compilation preparations
Activate the `intel` module set (which was created using instructions [here](Expanse_Mod_Setup.md#11-setting-up-the-intel-module-set)). This activation can be done by issuing:
```
module restore intel
```

If you have compiled your WRF model following the [instructions in Section 2](Expanse_Compile_WRF_Stuff.md), you should have a local NetCDF4 library on your $HOME directory. 
If you don't have the local NetCDF4 library created on your $HOME directory, [follow these instructions](Expanse_Compile_WRF_Stuff.md#step-213-set-up-a-local-netcdf4-library).


Once you have ensured that you have local NetCDF library on your $HOME directory, go back to the `PSU-EnKF_w_satellites/fortran_src/EnSRF/src` directory. Then, symbolically link the appropriate Makefile by issuing:
```
rm Makefile
ln -s Makefile.Expanse Makefile
```

### Step 4.3 Compile the PSU-EnKF EnSRF algorithm code
To compile the EnSRF code, just issue:
```
make clean
make
```
The compilation process should take less than 5 minutes (usually <1 min on Expanse).

If the compilation is successful, an `enkf.mpi` executable is created.

##
&nbsp;

[Click here to return to the main navigation page](../../../README.md)

[Click here to see all Expanse-specific PSU-EnKF setup steps](Setting_Up_On_Expanse.md)


