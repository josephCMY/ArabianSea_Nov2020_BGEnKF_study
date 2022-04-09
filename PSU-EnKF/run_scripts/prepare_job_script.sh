#!/bin/bash
id=$1

#=============TIME CONTROL PART=============
echo "#!/bin/bash"
cat << EOF
#####header for stampede######
#SBATCH -J WRF$id
#SBATCH -N 2
#SBATCH -n 64
#SBATCH -p normal
#SBATCH -t 36:00:00
#SBATCH -o out_wrf
#SBATCH -e error_wrf
#SBATCH --mail-user=mum373@psu.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

ibrun tacc_affinity ./wrf.exe

EOF

