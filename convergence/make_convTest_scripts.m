clear all
close all
clc

%% **** Set Up Experiment **************************************

% type names
typeNames = [{'all'},{'dry'},{'wet'}];
Ntypes = length(typeNames);

% target names
targNames = [{'Qe'},{'Qh'},{'NEE'}];
Ntargs = length(targNames);

% number of simultaneous jobs
Nprocs = Ntargs;   

%% **** Make PODS Scripts **************************************

% loop through types 
for iType = 1:Ntypes

 % exec file name
 fpods = strcat('convergence_bootstrap_',typeNames{iType},'.pods');

 % open file for writing
 fid = fopen(fpods,'w');

 % write instructions
 for iTarg = 1:Ntargs 
  cmd = strcat({'sh run_convergence_bootstrap.sh '},num2str(iType),{' '},num2str(iTarg)); 
  fprintf(fid,'%s \n',cmd{1});
 end

 % close file
 fclose(fid);

end

%% **** Make SLURM Scripts ********************************

% loop through types 
for iType = 1:Ntypes

 % exec file name
 fpods = strcat('convergence_bootstrap_',typeNames{iType},'.pods');

 % slurm file name 
 fslurm = strcat('convergence_bootstrap_',typeNames{iType},'.slurm');

 % open file
 fid = fopen(fslurm,'w');

 % write lines
 fprintf(fid,'#!/bin/bash -x \n \n');

 fprintf(fid,'#SBATCH --job-name=conv_%d \n',iType);
 fprintf(fid,'#SBATCH --output=reports/convBoot_%s_slurm.out \n',typeNames{iType});
 fprintf(fid,'#SBATCH --time=24:00:00 \n');
 fprintf(fid,'#SBATCH --ntasks=%d \n',Nprocs);
 fprintf(fid,'#SBATCH --account=s1688 \n');
 fprintf(fid,'#SBATCH --constraint=hasw \n');
 fprintf(fid,'#SBATCH --qos=long \n');
 fprintf(fid,' \n');  

 fprintf(fid,'cd /discover/nobackup/projects/summa/plumber/convergence \n');
 fprintf(fid,' \n');  

 str = strcat('/usr/local/other/PoDS/PoDS/pods.py -x "/discover/nobackup/projects/summa/plumber/convergence/',fpods,'" -n ');
 fprintf(fid,'%s %d \n',str,Nprocs);
 fprintf(fid,' \n');

 fprintf(fid,'exit 0 \n');
 fprintf(fid,' \n');

 % close file
 fclose(fid);

end

%% **** Make Dependency File *************************

% file name
fname = strcat('run_dependent_jobs.slurm');

% open file
fid = fopen(fname,'w');

% header lines
fprintf(fid,'#!/bin/bash -x \n');

% job counter
jobID = 0;

% write to file
for iType = 1:Ntypes

 % augment job counter
 jobID = jobID + 1;

 % slurm file name 
 fslurm = strcat('convergence_bootstrap_',typeNames{iType},'.slurm');

 % dependency name
 if jobID > 1; dname = jname; end

 % job name
 jname = strcat('JOB',num2str(jobID));

 % write job line to file
 if jobID == 1;
  fprintf(fid,'%s=$(sbatch %s | cut -f 4 -d'' '') \n',jname,fslurm);
 else
  fprintf(fid,'%s=$(sbatch -d afterany:$%s %s | cut -f 4 -d'' '') \n',jname,dname,fslurm);
 end

 % run it
 fprintf(fid,'echo $%s \n',jname);

end

% close file
fprintf(fid,'squeue -u gnearing \n');
fclose(fid);



