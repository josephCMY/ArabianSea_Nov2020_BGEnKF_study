# Basic step 2: Preparing a WRF ensemble for the PSU-EnKF
Author: Man-Yau (Joseph) Chan

[Click here to return to the main navigation page](../../README.md)

[Click here to return to the list of basic steps](../../README.md#basic-steps-to-use-the-psu-enkf-system)

&nbsp;


## Description

Before running the PSU-EnKF system, a WRF ensemble needs to be created. 
For this bare-bones system, I have written a system to construct the ensemble using the NCEP Global Ensemble Prediction System (GEFS). 
This construction system downloads the GEFS data archived on the Amazon Web Service and requires `rclone` to be set up.

**Attribution:** The construction system here is based off GEFS processing scripts that Yunji (Jerry) Zhang provided. What I did was rewrite his scripts for simplicity, readability and compatibility with my bare-bones system.

&nbsp;

**IMPORTANT NOTE: The construction system here is designed to work on Expanse.**

To adapt the construction system to work on other clusters, the scripts in `PSU-EnKF_w_satellites/system_types/barebones/setup_ens/gefs_based_ens` needs to be adjusted.



&nbsp;

## Instructions to construct GEFS-based WRF ensemble

1) [For first-time users: setting up shop](expanse_specific/Expanse_set_up_GEFS_system.md)
2) Navigate to `~/PSU-EnKF_w_satellites/system_types/barebones/setup_ens/gefs_based_ens`
3) Construct three WPS- and WRF-related namelists: one to run the WPS (example [here](../../setup_ens/gefs_based_ens/namelist.wps)), one to run `real.exe` (example [here](../../setup_ens/gefs_based_ens/namelist.real)) and one to perform the ensemble spin-up (example [here](../../setup_ens/gefs_based_ens/namelist.spinup))
4) Adjust the `config_ens_setup` file. (Link [here](../../setup_ens/gefs_based_ens/config_ens_setup))
5) Submit the `generate_ensemble_from_GEFS.sh` script as a SLURM job (just run `sbatch generate_ensemble_from_GEFS.sh`)

&nbsp;
##

[Click here to return to the main navigation page](../../README.md)

[Click here to return to the list of basic steps](../../README.md#basic-steps-to-use-the-psu-enkf-system)


