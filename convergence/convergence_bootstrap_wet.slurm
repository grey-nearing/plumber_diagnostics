#!/bin/bash -x 
 
#SBATCH --job-name=conv_3 
#SBATCH --output=reports/convBoot_wet_slurm.out 
#SBATCH --time=24:00:00 
#SBATCH --ntasks=3 
#SBATCH --account=s1688 
#SBATCH --constraint=hasw 
#SBATCH --qos=long 
 
cd /discover/nobackup/projects/summa/plumber/convergence 
 
/usr/local/other/PoDS/PoDS/pods.py -x "/discover/nobackup/projects/summa/plumber/convergence/convergence_bootstrap_wet.pods" -n 3 
 
exit 0 
 
