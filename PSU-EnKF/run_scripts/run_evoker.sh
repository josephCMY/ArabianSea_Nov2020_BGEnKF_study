#!/bin/bash
# Script to make the run script for each experiment
# Will also submit the experiment into a desired queue

# Request for user inputs
echo "State config file path: "
read config_file
echo "State queue: "
read queue
echo "Number of nodes: "
read nnodes
echo "Number of tasks: "
read ntasks
echo "Wall time (hh:mm:ss): "
read walltime


# Read in config file stuff
. $config_file
config_file=`pwd`/$config_file

# Making a copy of the template run script
cp run_template.sh run_$EXP_NAME.sh
runfile=run_$EXP_NAME.sh

# Editing the run script
sed -i "s/SBATCH_EXP_NAME/$EXP_NAME/g" $runfile
sed -i "s/SBATCH_NNODES/$nnodes/g" $runfile
sed -i "s/SBATCH_NTASKS/$ntasks/g" $runfile
sed -i "s/SBATCH_TIME/$walltime/g" $runfile
sed -i "s/SBATCH_LOG/$EXP_NAME/g" $runfile
sed -i "s/SBATCH_ERR/$EXP_NAME/g" $runfile
sed -i "s|CONFIG_FILE_PATH|$config_file|g" $runfile
sed -i "s/SBATCH_QUEUE/$queue/g" $runfile


