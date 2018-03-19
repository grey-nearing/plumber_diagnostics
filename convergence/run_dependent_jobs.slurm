#!/bin/bash -x 
JOB1=$(sbatch convergence_bootstrap_all.slurm | cut -f 4 -d' ') 
echo $JOB1 
JOB2=$(sbatch -d afterany:$JOB1 convergence_bootstrap_dry.slurm | cut -f 4 -d' ') 
echo $JOB2 
JOB3=$(sbatch -d afterany:$JOB2 convergence_bootstrap_wet.slurm | cut -f 4 -d' ') 
echo $JOB3 
squeue -u gnearing 
