#!/bin/bash
#####header for stampede######
#SBATCH -J L96_roi
#SBATCH -N 3 
#SBATCH -n 120
#SBATCH -p skx-normal
#SBATCH -t 02:00:00
#SBATCH -o test.log
#SBATCH -e test.err


# Script to setup and run many iterations of the experiment using TACC launcher 

intlist="4 6 8 10 12 14 16"
roilist=`seq 1 20`
ointlist="1 2 4"

# Load launcher
module load launcher

# Configure launcher
EXE=$TACC_LAUNCHER_DIR/init_launcher
PRUN=$TACC_LAUNCHER_DIR/paramrun
export LAUNCHER_RMI=SLURM
export LAUNCHER_PLUGIN_DIR=$LAUNCHER_DIR/plugins
export LAUNCHER_SCHED=interleaved
export LAUNCHER_JOB_FILE=launcher_script.sh
export LAUNCHER_WORKDIR=`pwd`

# Iterative loop to set up and run python scripts
for roi in $roilist; do
  for nint in $intlist; do
    rm launcher_script.sh
    touch launcher_script.sh
    for oint in $ointlist; do
      for ii in `seq -f "%04g" 1 40`; do
        echo "python3 roi_sensitivities.py $nint $oint 8 $roi $ii >& outputs/log."$nint"_"$oint"_"$roi"_$ii" 
      done 
    done > launcher_script.sh
    # Start launcher
    $PRUN 
  done
done


exit
