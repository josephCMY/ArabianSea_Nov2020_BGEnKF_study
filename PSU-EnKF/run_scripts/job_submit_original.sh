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
  #ibrun tacc_affinity -n $n -o $o $exe << this is the command I tried
  ibrun -n $n -o $o task_affinity $exe
  export SLURM_TASKS_PER_NODE="$((SLURM_NTASKS/$SLURM_NNODES))(x$SLURM_NNODES)"
fi

### NERSC Cori
# Specific to using knl
if [[ $HOSTTYPE == 'cori-knl' ]]; then
  nnodes=$(($n/$ppn))
  if [[ $(($nnodes * $ppn)) -lt $n ]]; then
    nnodes=$(($nnodes+1))
  fi
  srun -N $nnodes -n $n --ntasks-per-node=$ppn --cpus-per-task=1 --cpu_bind=cores $exe 
  echo Executing $exe with $nnodes nodes $n procs and $ppn procs-per-node 
fi  




  script="
#!/bin/bash

#####header for Cori######
#SBATCH -J $exe 
#SBATCH -q normal
#SBATCH -C knl
#SBATCH -N $nnodes
#SBATCH -n $n
#SBATCH --time=$runtime
#SBATCH -o log.$exe
#SBATCH -e err.$exe

# OpenMP settings
export OMP_NUM_THREADS=1
export OMP_PLACES=cores
export OMP_PROC_BIND=true


srun -N $nnodes -n $n -c 2 --cpu_bind=cores $exe
fi

"

