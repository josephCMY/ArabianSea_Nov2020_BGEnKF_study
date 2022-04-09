# Quick script to set up new experiments

# Path to directory containing scripts and what not
workdir=~/nonlinear_IR-DA/indian_ocean_osse/PSU_EnKF

for expt in JointSpace_EnKF_IR_only  JointSpace_EnKF_WV_only  JointSpace_BGEnKF_IR_only  StateSpace_EnKF_IR_only  StateSpace_EnKF_WV_only NoDA; do

  # Set up experiment directory
  expt_dir=/global/cscratch1/sd/my_chan/nonlinear_IR-DA/indian_ocean_osse/PSU_EnKF/$expt
  mkdir $expt_dir

  # Enter experiment directory
  cd $expt_dir

  # Symbolic links to the various directories
  ln -s $workdir/bdy .
  ln -s $workdir/scripts DA
  ln -s $workdir/code .
  ln -s $workdir/data .

  # Set up ensemble files
  mkdir -p fc/201110150000
  cd fc/201110150000
  cp $workdir/../setup_ens/determine_member_roles/firstcycle_ens/* .

  echo finished setting up $expt 

done 
