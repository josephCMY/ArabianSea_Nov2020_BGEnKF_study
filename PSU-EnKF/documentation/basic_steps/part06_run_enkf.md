# Part 6: Setting up the directory system on your cluster's SCRATCH space

[Click here to return to the main navigation page](../../README.md)

[Click here to return to the list of basic steps](../../README.md#basic-steps-to-use-the-psu-enkf-system)

&nbsp;


## Description and step-by-step instructions
I typically compile the 6 programs listed in [part 1](#part-1-setting-up-the-environment-and-compiling-relevant-programs) on my cluster's WORK or HOME space. I also often store my PSU-EnKF's observations, ensemble initial conditions, and ensemble boundary conditions on my cluster's WORK or HOME space.

However, **it is usually unadvisable to run PSU-EnKF experiments on your WORK or HOME space**. The reason is that the PSU-EnKF tends to generate terabytes of data (even my smallest experiments generate ~8 TB each) and the WORK and HOME spaces tend to have < 2 TB. As such, I almost always run my PSU-EnKF experiments on my cluster's SCRATCH space.


In the instructions to follow, I will assume that:
1) The path to the `barebones` directory (_i.e._, the directory of this README.md file you are reading) is stored in the bash variable `PATH_TO_BAREBONES`.
2) You have successfully compiled the programs in [part 1](#part-1-setting-up-the-environment-and-compiling-relevant-programs) and these programs live within `$PATH_TO_BAREBONES/fortran_src`.
3) The bash variable `SCRATCH` contains the path to your SCRATCH space.
4) You have followed the earlier instructions to preprocess the observations and the ensemble's initial and boundary conditions.
5) You have followed the earlier instructions to spin-up your ensemble.
5) The name of your experment is `new_expt`.
6) The desired path to your new experiment directory is `$SCRATCH/some_path/new_expt`
7) The start date of your ensemble is `DATE_START`


### Step 1: Make the parent directory for your new experiment on SCRATCH
Just issue
```
mkdir -p $SCRATCH/some_path/new_expt
```

### Step 2: Symbolic link the compiled programs, bash scripts and prepared data from `$PATH_TO_BAREBONES`
Enter the experiment directory via
```
cd $SCRATCH/some_path/new_expt
```

Then issue
```
ln -s $PATH_TO_BAREBONES/fortran_src code
ln -s $PATH_TO_BAREBONES/run_scripts DA
ln -s $PATH_TO_BAREBONES/data data
ln -s $PATH_TO_BAREBONES/bdy bdy
```


### Step 3: Copying the spun-up WRF ensemble to SCRATCH
Finally, construct the directory containing the initial and spun-up ensemble WRF files by issuing
```
mkdir -p $SCRATCH/some_path/new_expt/fc/$DATE_START
cd $SCRATCH/some_path/new_expt/fc/$DATE_START
cp $PATH_TO_BAREBONES/ens_spinup/wrfinput_d01* .
```



&nbsp;

##
[Click here to return to the main navigation page](../../README.md)

[Click here to return to the list of basic steps](../../README.md#basic-steps-to-use-the-psu-enkf-system)

