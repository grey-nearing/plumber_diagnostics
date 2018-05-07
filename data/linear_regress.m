%% --- Runtime Environment ------------------------------------------------

clear all; close all; clc;
restoredefaultpath; addpath(genpath(pwd));
fignum = 0;

%% --- Load Data ----------------------------------------------------------

% screen report
fprintf('Loading data ... '); tic;

% load fluxnet data
load('./data/allDailyData.mat');  % fluxnet daily concatenated matrix

% screen report
fprintf('. finished; time = %f \n',toc); 

%% --- Build Linear Models ------------------------------------------------

cols= [2:9,12:23];

for s = 1:Nsites
    
  X = allData(:,cols,s);
  y = allData(:,12,s);
  
  model{s} = stepwiselm([X,y],'interactions');
    
    
end
