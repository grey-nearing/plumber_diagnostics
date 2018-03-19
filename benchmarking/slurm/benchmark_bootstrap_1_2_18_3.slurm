#!/bin/bash 
 
#SBATCH --job-name=2_18_3 
#SBATCH --output=reports/all_Qh_Sylvania_3_slurm.out 
#SBATCH --time=12:00:00 
#SBATCH -N 1 
#SBATCH --account=s1688 
#SBATCH --constraint=hasw 
##SBATCH --qos=long 
 
cd /discover/nobackup/projects/summa/misc_projects/plumber/benchmarking 
./run_benchmark_bootstrap.sh 1 2 18 3 
exit 0 
 
