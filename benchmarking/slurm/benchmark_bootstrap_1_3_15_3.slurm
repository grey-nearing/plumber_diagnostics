#!/bin/bash 
 
#SBATCH --job-name=3_15_3 
#SBATCH --output=reports/all_NEE_Merbleue_3_slurm.out 
#SBATCH --time=12:00:00 
#SBATCH -N 1 
#SBATCH --account=s1688 
#SBATCH --constraint=hasw 
##SBATCH --qos=long 
 
cd /discover/nobackup/projects/summa/misc_projects/plumber/benchmarking 
./run_benchmark_bootstrap.sh 1 3 15 3 
exit 0 
 
