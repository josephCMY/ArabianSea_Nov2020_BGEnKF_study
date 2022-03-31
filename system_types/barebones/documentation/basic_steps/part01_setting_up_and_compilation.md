# Basic step 1: Setting up the environment and compiling relevant programs


[Click here to return to the main navigation page](../../README.md)

[Click here to return to the list of basic steps](../../README.md#basic-steps-to-use-the-psu-enkf-system)

&nbsp;

## Description
The following programs need to be compiled to run the PSU-EnKF:
1) **The Weather Research and Forecast model (WRF)**: Used as the PSU-EnKF's forecast model
2) **The WRF Preprocessing System (WPS)**: Used to process initial and boundary conditions to run WRF
3) **The WRF Data Assimilation (WRFDA) system**: Used mainly to process Global Telecommunication System (GTS) conventional observation files into something the PSU-EnKF can ingest.
4) **The Community Radiative Transfer Model (CRTM)**: Used to simulate infrared and microwave radiance observations from space-borne instruments. Necessary to assimilate such observations.
5) **QN2RE**: Enhances the accuracy of cloud-affected microwave radiance simulations from the CRTM. 
6) **The PSU-EnKF EnSRF code**: Actual program used to assimilate observations into an ensemble of WRF model outcomes.

Generally speaking, you will need to set up the compiler and library environments on the cluster/computer before compiling these 6 items. 
These 6 items can be compiled on your WORK or HOME space (as long as you have at least 1 GB of free space available on your WORK/HOME space).
&nbsp;

**For Expanse users:** Click [here](expanse_specific/Setting_Up_On_Expanse.md) for step-by-step instructions to set up the PSU-EnKF on Expanse.


**For non-Expanse users:** Refer to the instructions for setting up on Expanse, and adjust the procedure for your cluster/computer. Note that you might need to use different C and Fortran compilers.

&nbsp;

##
[Click here to return to the main navigation page](../../README.md)

[Click here to return to the list of basic steps](../../README.md#basic-steps-to-use-the-psu-enkf-system)


