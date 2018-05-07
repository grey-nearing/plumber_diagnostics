clear all; close all; clc
restoredefaultpath; addpath(genpath(pwd));
fignum = 0;

%% **** Set Up Experiment *************************************************

% number of bootstraps
Nboots = 1;

% type names
typeNames = [{'all'},{'dry'},{'wet'}];%
Ntypes = length(typeNames);

% target names
targNames = [{'Qe'},{'Qh'},{'NEE'}];
Ntargs = length(targNames);

% site names
siteNames = [{'Amplero'},{'Blodgett'},{'Bugac'   },{'ElSaler2'},{'ElSaler' },...
             {'Espirra'},{'FortPeck'},{'Harvard' },{'Hesse'   },{'Howlandm'},...
             {'Howard' },{'Hyytiala'},{'Kruger'  },{'Loobos'  },{'Merbleue'},...
             {'Mopane' },{'Palang'  },{'Sylvania'},{'Tumba'   },{'UniMich' }];
Nsites = length(siteNames);

% number of simultaneous jobs
Nprocs = 1;   

%% **** Make SLURM Scripts ************************************************

% loop through types 
for iType = 2:Ntypes
 for iBoot = 1:Nboots
  for iTarg = 1:Ntargs 
   for iSite = 1:Nsites

    % slurm file name
    fslurm = strcat('slurm/benchmark_bootstrap_',num2str(iType),'_',...
        num2str(iTarg),'_',num2str(iSite),'_',num2str(iBoot),'.slurm');

    % open file
    fid = fopen(fslurm,'w');

    % write lines
    fprintf(fid,'#!/bin/bash \n \n');

    fprintf(fid,'#SBATCH --job-name=%d_%d_%d \n',iTarg,iSite,iBoot);
    fprintf(fid,'#SBATCH --output=reports/%s_%s_%s_%d_slurm.out \n',...
        typeNames{iType},targNames{iTarg},siteNames{iSite},iBoot);
    fprintf(fid,'#SBATCH --time=12:00:00 \n');
    fprintf(fid,'#SBATCH -N %d \n',1);
    fprintf(fid,'#SBATCH --account=s1688 \n');
    fprintf(fid,'#SBATCH --constraint=hasw \n');
    fprintf(fid,'##SBATCH --qos=long \n');
    fprintf(fid,' \n');  

    fprintf(fid,strcat({'cd '},pwd,{' \n'}));
    fprintf(fid,'./run_benchmark_bootstrap.sh %d %d %d %d \n',iType,...
        iTarg,iSite,iBoot);
    fprintf(fid,'exit 0 \n');
    fprintf(fid,' \n');

    % close file
    fclose(fid);

   end
  end
 end
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
 for iBoot = 1:Nboots
  for iTarg = 1:Ntargs 
   for iSite = 1:Nsites

    % augment job counter
    jobID = jobID + 1;

    % slurm file name
    fslurm = strcat('slurm/benchmark_bootstrap_',num2str(iType),'_',...
        num2str(iTarg),'_',num2str(iSite),'_',num2str(iBoot),'.slurm');

    % dependency name
    if jobID > 1; dname = jname; end

    % job name
    jname = strcat('JOB',num2str(jobID));

    % write job line to file
    if jobID == 1
     fprintf(fid,'%s=$(sbatch %s | cut -f 4 -d'' '') \n',jname,fslurm);
    else
     fprintf(fid,'%s=$(sbatch -d afterany:$%s %s | cut -f 4 -d'' '') \n',...
         jname,dname,fslurm);
    end

    % run it
    fprintf(fid,'echo $%s \n',jname);

   end
  end
 end
end

% close file
fprintf(fid,'squeue -u gnearing \n');
fclose(fid);

%% *** END SCRIPT *********************************************************





