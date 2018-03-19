#!/bin/bash 
 
#SBATCH --job-name=3_18_1 
#SBATCH --output=reports/all_NEE_Sylvania_1_slurm.out 
#SBATCH --time=12:00:00 
#SBATCH -N 1 
#SBATCH --account=s1688 
#SBATCH --constraint=hasw 
##SBATCH --qos=long 
 
cd /discover/nobackup/projects/summa/misc_projects/plumber/benchmarking 
./run_benchmark_bootstrap.sh 1 3 18 1 
exit 0 
 
