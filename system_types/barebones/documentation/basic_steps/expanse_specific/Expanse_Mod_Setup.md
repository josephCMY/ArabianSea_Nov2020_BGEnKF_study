# Section 1: Setting up Expanse modules for first-time users
Written by Man-Yau (Joseph) Chan, based on Yunji (Jerry) Zhang's emailed instructions.


[Click here to return to the main navigation page](../../../README.md)

[Click here to see all Expanse-specific PSU-EnKF setup steps](Setting_Up_On_Expanse.md)

&nbsp;


## IMPORTANT: On the need to have two sets of Expanse modules
To run WRF-related programs, the Community Radiative Transfer Model (CRTM), and the PSU-EnKF's Fortran programs, a set of modules relating to the Intel compilers is needed.

However, this set of Intel modules is incompatible with some commonly-used software (e.g., NCAR Command Line (`ncl`), NetCDF View (`ncview`) and the NetCDF Operators (`nco`)). 

It is thus convenient to have two sets of modules: 
1) The `intel` module set: for running WRF-related programs, the CRTM, and the PSU-EnKF's Fortran programs.
2) The `default` module set: for running day-to-day analysis codes.

In this section, we will describe how the two module sets are created.


## 1.1: Setting up the Intel module set
Issue the following commands in your bash terminal:
```
module reset
module load intel
module load intel-mpi
module load intel-mkl
module load libtirpc
module load hdf5
module load netcdf-c
module load netcdf-fortran
module load parallel-netcdf
module load sdsc
```

To save the current set of modules as a module set called `intel`, issue: 
```
module save intel
```

Whenever you need to compile or run WRF-related programs, the CRTM, or PSU-EnKF Fortran programs, you need to activate this `intel` module set. To activate the `intel` module set, issue:
```
module restore intel
```

## 1.2: Setting up the Default module set for day-to-day analysis codes
Issue the following commands in your bash terminal
```
module reset
module load gcc
module load openmpi
module load libtirpc
module load hdf5
module load netcdf-c
module load netcdf-fortran
module load parallel-netcdf
module load sdsc
```
To save this set of modules as the `default` set, issue:
```
module save default
```

Note that whenever you log into Expanse, this `default` set is automatically loaded.


Whenever you want to run NCL, NCO, or NCVIEW, make sure the `default` module set is active.


## 1.3 Setting up Anaconda to use Joseph's test case
Anaconda is a convenient package management tool for Python users. While Expanse has an Anaconda module, the module is out-of-date. As such, **do not use Expanse's Anaconda module**!
Instead, Python users should install Anaconda themselves.

To install Anaconda, go to your $HOME directory by issuing:
```
cd
```

Then, download the Anaconda installer for Linux by issuing:
```
wget https://repo.anaconda.com/archive/Anaconda3-2021.11-Linux-x86_64.sh
```

To install Anaconda, just execute
```
bash Anaconda3-2021.11-Linux-x86_64.sh
```
Then, follow the instructions and the prompts that appear on your Terminal screen. 
When Anaconda asks you if you want to initialize the base Anaconda environment, please say yes.
For more detailed instructions, visit the [Anaconda installation guide](https://docs.anaconda.com/anaconda/install/linux).

Then, exit Expanse and re-enter it.

## 1.4 Preparing Python package to run Joseph's test case
You will need the netCDF4 Python package and the associated Python packages.

First, update your Anaconda by issuing:
```
conda update -n base -c defaults conda
```

Then, install the packages by issuing:
```
conda install -c conda-forge netcdf4
```
and let Anaconda run.

At some point, you will be prompted by the terminal about whether to install a whole list of packages. Just key in `y` and hit Enter, then let Anaconda do its own thing.

##

&nbsp; 

[Click here to return to the main navigation page](../../../README.md)

[Click here to see all Expanse-specific PSU-EnKF setup steps](Setting_Up_On_Expanse.md)

