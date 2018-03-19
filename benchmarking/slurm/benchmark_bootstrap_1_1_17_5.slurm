#!/bin/bash 
 
#SBATCH --job-name=1_17_5 
#SBATCH --output=reports/all_Qe_Palang_5_slurm.out 
#SBATCH --time=12:00:00 
#SBATCH -N 1 
#SBATCH --account=s1688 
#SBATCH --constraint=hasw 
##SBATCH --qos=long 
 
cd /discover/nobackup/projects/summa/misc_projects/plumber/benchmarking 
./run_benchmark_bootstrap.sh 1 1 17 5 
exit 0 
 
