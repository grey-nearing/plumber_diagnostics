#!/bin/bash 
 
#SBATCH --job-name=2_20_2 
#SBATCH --output=reports/all_Qh_UniMich_2_slurm.out 
#SBATCH --time=12:00:00 
#SBATCH -N 1 
#SBATCH --account=s1688 
#SBATCH --constraint=hasw 
##SBATCH --qos=long 
 
cd /discover/nobackup/projects/summa/misc_projects/plumber/benchmarking 
./run_benchmark_bootstrap.sh 1 2 20 2 
exit 0 
 
