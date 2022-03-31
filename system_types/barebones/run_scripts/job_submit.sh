#!/bin/bash
#HOSTTYPE, HOSTPPN are defined in config file
#total_ntasks is defined in run.sh
. $CONFIG_FILE
n=$1  # num of tasks job uses
o=$2  # offset location in total_ntasks, useful for several jobs to run together
ppn=$3  # proc per node for the job
exe=$4  # executable

### TACC Stampede2
if [[ $HOSTTYPE == "stampede" ]]; then
  export SLURM_TASKS_PER_NODE="$ppn(x$SLURM_NNODES)"
  ibrun -n $n -o $o task_affinity $exe
  export SLURM_TASKS_PER_NODE="$((SLURM_NTASKS/$SLURM_NNODES))(x$SLURM_NNODES)"
fi


### XSEDE Expanse cluster
if [[ $HOSTTYPE == "expanse" ]]; then
  export SLURM_TASKS_PER_NODE="$ppn(x$SLURM_NNODES)"
  module restore intel
  ibrun -n $n -o $o $exe
  export SLURM_TASKS_PER_NODE="$((SLURM_NTASKS/$SLURM_NNODES))(x$SLURM_NNODES)"
  module restore default
fi


### NERSC Cori

# If we are runing the PSU EnKF system from within a Cori KNL node
if [[ $HOSTTYPE == 'cori-knl' ]]; then
  nnodes=$(($n/$ppn))
  if [[ $(($nnodes * $ppn)) -lt $n ]]; then
    nnodes=$(($nnodes+1))
  fi
  srun -N $nnodes -n $n --ntasks-per-node=$ppn --cpus-per-task=1 --cpu_bind=cores $exe 
  echo Executing $exe with $nnodes nodes $n procs and $ppn procs-per-node 
fi  


# If we are running the PSU EnKF system from within a Cori login node
if [[ $HOSTTYPE == 'cori-login' ]] ; then

  # Ensuring the correct number of nodes will be called
  nnodes=$(($n/$ppn))
  if [[ $(($nnodes * $ppn)) -lt $n ]]; then
    nnodes=$(($nnodes+1))
  fi
 
  max_run_time=60
  # Determining the maximum amount of runtime to allocate
  if [[ $exe == 'enkf.mpi' ]] ; then
    max_run_time=60
  fi
  if [[ $exe == 'wrf.exe' ]]; then
    max_run_time=60
  fi

  # Execute
  salloc -C knl -N $nnodes -n $n -q interactive -t $max_run_time --job-name="$exe" srun -N $nnodes --ntasks-per-node=$ppn --cpus-per-task=1 --cpu_bind=cores $exe

fi



# If we are running the PSU EnKF system from within a Cori login node, but using regular queue
if [[ $HOSTTYPE == 'cori-regular' ]] ; then

  # Ensuring the correct number of nodes will be called
  nnodes=$(($n/$ppn))
  if [[ $(($nnodes * $ppn)) -lt $n ]]; then
    nnodes=$(($nnodes+1))
  fi
 
  max_run_time=60
  # Determining the maximum amount of runtime to allocate
  if [[ $exe == 'enkf.mpi' ]] ; then
    max_run_time=60
  fi
  if [[ $exe == 'wrf.exe' ]]; then
    max_run_time=60
  fi

  # Execute
  salloc -C knl -N $nnodes -n $n -q regular -t $max_run_time --job-name="$exe" srun -N $nnodes --ntasks-per-node=$ppn --cpus-per-task=1 --cpu_bind=cores $exe

fi
