# Setting up the PSU-EnKF for XSEDE's Expanse cluster

## Outline

**Major steps:**
1) Download and compile WRF using Intel compilers
2) Compile the CRTM and QN2RE
3) Compile the PSU-EnKF


**Important notes:**
This set of instructions are meant for first-time users of Expanse.
If you have already set up your own environment on Expanse, you can either:
1) "Factory reset" your environment (will include instructions later)
2) Fiddle with your environment to get things to work.




##  







## Part A: Preparing the user environment



### 2) For Python users: download and install Anaconda in your $HOME directory
Anaconda is a package management system commonly used by Python users. 
Anaconda makes it pretty easy to install new Python packages without needing administrator privileges.

To install Anaconda, start by issuing the following commands:

```
cd $HOME
wget https://repo.anaconda.com/archive/Anaconda3-2021.11-Linux-x86_64.sh
bash Anaconda3-2021.11-Linux-x86_64.sh
```

Upon running the third command (`bash Anaconda3-2021.11-Linux-x86_64.sh`), a license agreement will appear. 
Keep pressing the Enter key to scroll through the agreement. 

Then something like this will appear:
```
Do you accept the license terms? [yes|no]
```
Key in `yes`, hit enter twice, and the Anaconda installation process will begin.

Note that the Anaconda installation may take up to an hour to complete.

Once the installation is complete, the following messages and command prompt will appear:
```
installation finished.
Do you wish the installer to initialize Anaconda3
by running conda init? [yes|no]
[no] >>>
```
Key in `yes` and hit the Enter key. 

Then, you should see the following message.

```==> For changes to take effect, close and re-open your current shell. <==```

Exit Expanse and re-enter it.

Upon re-entering Expanse via a terminal, your command prompt should look like this:
```
(base) [chanmy@login01 ~]$ 
```
If you see the `(base)` to the left of the command prompt, you have succeeded in installing and initializing Anaconda.


### 3) Generating a new Anaconda environmen
Because I am paranoid and like redundancy, I recommend creating a new environment in case anything goes south. 

First, issue:
```
conda create --name myenv
```
When the terminal asks you to whether it should proceed, key in `y` and hit enter.


Then, issue the following command to ensure that the Anaconda environment `myenv` activates upon logging into Expanse:
```
cat "conda activate myenv" >> ~/.bashrc
```


### 3) Installing NCL and NCO
NCL and NCO are commonly used in meteorological research. NCL stands for the NCAR Command Line and NCO stands for NetCDF Operators.

Now that Anaconda is installed, it is pretty easy to install NCL and NCO. To install NCL, issue:
```
conda install -c conda-forge ncl
```
The terminal will then display a list of new packages that will be installed. Then, the terminal will prompt you for an input:
```Proceed ([y]/n)?```
Hit the `y` key and the enter key, then wait.

Thncl






## Part B: Compiling WRF

```
module reset
module load intel
module load intel-mpi
module load libtirpc
module load hdf5
module load netcdf-c
module load netcdf-fortran
module load parallel-netcdf
```

