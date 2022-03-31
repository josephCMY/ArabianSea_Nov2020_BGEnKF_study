# Section 3: Setting up CRTM and QN2RE on Expanse
Authors: Man-Yau (Joseph) Chan and Yunji (Jerry) Zhang

[Click here to return to the main navigation page](../../../README.md)

[Click here to see all Expanse-specific PSU-EnKF setup steps](Setting_Up_On_Expanse.md)

&nbsp;
## Overview

The PSU-EnKF uses the Community Radiative Transfer Model (CRTM) to assimilate 
radiance observations from satellite-borne infrared and microwave radiometers.
The QN2RE is a necessary add-on written by Jerry and Scott Sieron to improve 
the accuracy of CRTM-simulated cloudy-sky microwave brightness temperatures. 

In this section, we will cover compiling the CRTM and QN2RE. Here are the 
hyperlinks to navigate to those steps

3.1) [Instructions to compile the CRTM on Expanse](#31-compiling-the-crtm-on-expanse)

3.2) [Instructions to compile the QN2RE on Expanse](#32-compiling-the-qn2re-on-expanse)

&nbsp;
&nbsp;
&nbsp;
&nbsp;
&nbsp;
&nbsp;
&nbsp;
&nbsp;

## 3.1 Compiling the CRTM on Expanse

### Step 3.1.1 Create and activate the `intel` module set
If you haven't already followed the instructions [here](Expanse_Mod_Setup.md#11-setting-up-the-intel-module-set) to create an `intel` module set, start by following said instructions to do so.

Once you have the `intel` module set, activate the `intel` module set by issuing
```
module restore intel
```

### Step 3.1.2 Enter the CRTM directory 
```
cd PSU-EnKF_w_satellites/fortran_src/crtm_v2.3.0/
```

### Step 3.1.3 Download coefficients tarball file for CRTM using provided bash script.
```
./download_coefficients.sh
```

### Step 3.1.4 Untar the coefficients tarball file.
```
tar -xvf coefficients.*.tar
```

### Step 3.1.5 Setup relevant environment
Just issue
```
. config-setup/ifort.setup
```



### Step 3.1.6 Set up the compilation configuration
Issue this on the command line:
```
./configure --disable-big-endian --prefix=`pwd`
```

### Step 3.1.7 CRTM compilation and linking
Issue this on the command line to compile the CRTM
```
make; make install
```

Then check if file `crtm_v2.3.0/lib/libcrtm.a` exists. If this file exists, compilation is successful.

Finally, generate symbolic links needed for `EnSRF` to call CRTM by running
```
ln -s crtm_v2.3.0/* .
```

### Step 3.1.8 Closing shop
To revert from the `intel` module set to the default module set, issue:
```
module restore default
```

&nbsp;
&nbsp;
&nbsp;
&nbsp;
&nbsp;
&nbsp;
&nbsp;
&nbsp;



## 3.2 Compiling the QN2RE on Expanse

### Step 3.2.1 Create and activate the `intel` module set
If you haven't already followed the instructions [here](Expanse_Mod_Setup.md#11-setting-up-the-intel-module-set) to create an `intel` module set, start by following said instructions to do so.

Once you have the `intel` module set, activate the `intel` module set by issuing
```
module restore intel
```

### Step 3.2.2 Setup local NetCDF library on your $HOME directory
If you haven't done it already, follow the instructions [here](Expanse_Compile_WRF_Stuff.md#step-213-set-up-a-local-netcdf4-library) to create a local NetCDF4 library on your $HOME directory.


### Step 3.2.3 Enter the QN2RE directory
```
cd PSU-EnKF_w_satellites/fortran_src/qn2re/
```

### Step 3.2.4 Symbolic link the appropriate makefile for Expanse
Just issue:
```
rm makefile
ln -s makefile.expanse makefile
```

### Step 3.2.5 Compile QN2RE
Just issue:
```
make
```

If the file `libqn2re.a` is created by the `make` command, you have successfully compiled QN2RE.

### Step 3.2.6 Closing shop
Just issue the following command to switch from the `intel` module set back to the default module set.
```
module restore default
```

##

&nbsp;

[Click here to return to the main navigation page](../../../README.md)

[Click here to see all Expanse-specific PSU-EnKF setup steps](Setting_Up_On_Expanse.md)


