#!/bin/bash -x 
 
#SBATCH --job-name=conv_2 
#SBATCH --output=reports/convBoot_dry_slurm.out 
#SBATCH --time=24:00:00 
#SBATCH --ntasks=3 
#SBATCH --account=s1688 
#SBATCH --constraint=hasw 
#SBATCH --qos=long 
 
cd /discover/nobackup/projects/summa/plumber/convergence 
 
/usr/local/other/PoDS/PoDS/pods.py -x "/discover/nobackup/projects/summa/plumber/convergence/convergence_bootstrap_dry.pods" -n 3 
 
exit 0 
 
