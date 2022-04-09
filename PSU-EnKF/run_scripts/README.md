# Control scripts to perform EnKF experiments 

**Special notes**:

1) Cori doesnt seem to like running multiple instances of WRF at the same time in the same job. Ie, instead of running 3 instances of WRF at the same time, it seems to prefer doing 2 instances of WRF, then finish up the last instance. **Solution: submit each WRF integration as a job on SBATCH.**
