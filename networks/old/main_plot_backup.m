clear all
close all
clc
restoredefaultpath
addpath('../matlab_tools');


%% --- Experiment Setup ------------------------------------------

modelNames = [ ...
              {'CABLE.2.0'}
              {'CABLE_2.0_SLI.vxh599_r553'}
              {'CHTESSEL'}
              {'COLASSiB.2.0'}
              {'ISBA_SURFEX_3l.SURFEX7.3'}
              {'ISBA_SURFEX_dif.SURFEX7.3'}
              {'JULES.3.1'}
              {'JULES3.1_altP'}
              {'Mosaic.1'}
              {'NOAH.2.7.1'}
              {'Noah.3.2'}
              {'NOAH.3.3'}
%              {'ORCHIDEE.trunk_r1401'}
%              {'SUMMA.1.0.exp.01.000'}
%              {'SUMMA.1.0.exp.01.001'}
%              {'SUMMA.1.0.exp.01.002'}
%              {'SUMMA.1.0.exp.01.003'}
%              {'SUMMA.1.0.exp.01.004'}
%              {'SUMMA.1.0.exp.01.005'}
%              {'SUMMA.1.0.exp.01.006'}
%              {'SUMMA.1.0.exp.01.007'}
%              {'SUMMA.1.0.exp.01.008'}
%              {'SUMMA.1.0.exp.01.009'}
              {'Manabe_Bucket.2 '}
              {'Penman_Monteith.1'}
];

% site names
siteNames =  [...
             {'Amplero'}
             {'Blodgett'}
             {'Bugac'}
%             {'ElSaler2'}
             {'ElSaler'}
%             {'Espirra'}
             {'FortPeck'}
%             {'Harvard'}
             {'Hesse'}
%             {'Howlandm'}
             {'Howard'}
             {'Hyytiala'}
             {'Kruger'}
             {'Loobos'}
%             {'Merbleue'}
%             {'Mopane'}
             {'Palang'}
             {'Sylvania'}
%             {'Tumba'}
%             {'UniMich'}
];

% temporal lags
%lags = round(logspace(log10(1),log10(48*180),15));
lags = [1,2,6,12,24,48,96,128,256];

% number of sites
Nsites  = length(siteNames);
Nmodels = length(modelNames);
Nlags   = length(lags);

%% --- Load Transfer Results -------------------------------------

% variable names + dimensions
varNames = [{'U'},{'T_a'},{'R_n'},{'R_h'},{'P'},{'Q_e'},{'Q_h'},{'NEE'},{'S_w'},{'SM2'}];
Du = 5;
Dx = 5;
D = Du+Dx;

% init storage
Tm = zeros(D,Dx,Nlags,Nsites,Nmodels);
Tp = zeros(D,Dx,Nlags,Nsites);

% load file data
for iLag = 1
 for iSite = 1:Nsites

  for iModel = 1:Nmodels
   % read model results file
   fname = strcat('results/models/PALS_Normalized_',siteNames{iSite},'_',modelNames{iModel},'_',num2str(iLag),'.txt');
   Tm(:,:,iLag,iSite,iModel) = load_result_file(fname);

   % screen report
   fprintf('finished loading: lag %d/%d - site %d/%d - model %d/%d \n',iLag,Nlags,iSite,Nsites,iModel,Nmodels);
  end

  % read pals results file
  fname = strcat('results/pals/PALS_Normalized_',siteNames{iSite},'_',num2str(iLag),'.txt');
  Tp(:,:,iLag,iSite) = load_result_file(fname);

 end
end

% transfer differences
for iModel = 1:Nmodels
 for iSite = 1:Nsites
  TD(:,:,iSite,iModel) = squeeze(Tp(:,1:2,1,iSite)-Tm(:,1:2,1,iSite,iModel));
  I = find(abs(Tp(:,1,1,iSite       ))<1e-4); TD(I,1,iSite,iModel) = 0/0;
  I = find(abs(Tm(:,1,1,iSite,iModel))<1e-4); TD(I,1,iSite,iModel) = 0/0;
  I = find(abs(Tp(:,2,1,iSite       ))<1e-4); TD(I,2,iSite,iModel) = 0/0;
  I = find(abs(Tm(:,2,1,iSite,iModel))<1e-4); TD(I,2,iSite,iModel) = 0/0;
 end
end

%% --- Load Benchmarking Results ---------------------------------

% nubmer of bootstraps
Nboots = 2;

% benchmark targets
targNames = [{'Qe'},{'Qh'}];%,{'NEE'}];
Ntargs = length(targNames);

% loop through experiments
for iBoot = 1:Nboots
 for iTarg = 1:Ntargs
  for iSite = 1:Nsites

   % load benchmark 
   fname = strcat('../benchmarking/results/Site_all_',targNames{iTarg},'_',siteNames{iSite},'_',num2str(iBoot),'_10000_trainscg.mat');
   load(fname); % stores in variable 'results' 

   % benchmark data
%   Yloo = results.test.Yhat;
   Yloo = results.Yhat;

   % observation data
%   Yobs = results.test.Yobs;
   Yobs = results.Yobs;

   % load model
   fname = strcat('../data/lagged_data/models_',siteNames{iSite},'.mat');
   load(fname); % stores in variable 'model'

   % model data
   Ymod = squeeze(model(:,3+iTarg,:));

   % grandmas
   I = find(isnan(Yloo)); Yloo(I) = []; Yobs(I) = []; Ymod(I,:) = [];
   I = find(isnan(Yobs)); Yloo(I) = []; Yobs(I) = []; Ymod(I,:) = [];
   I = find(Yloo <-9900); Yloo(I) = []; Yobs(I) = []; Ymod(I,:) = [];
   I = find(Yobs <-9900); Yloo(I) = []; Yobs(I) = []; Ymod(I,:) = [];
   I = find(Yloo > 2000); Yloo(I) = []; Yobs(I) = []; Ymod(I,:) = [];
   I = find(Yobs > 2000); Yloo(I) = []; Yobs(I) = []; Ymod(I,:) = [];

   % info bins
   Bw = (max(Yobs)-min(Yobs))*0.02;
   Bhat = (min(Yloo)-Bw):Bw:(max(Yloo)+Bw);
   Bobs = (min(Yobs)-Bw):Bw:(max(Yobs)+Bw);

   % loop throug models
   for iModel = 1:Nmodels

    % segregate model data
    YYmod = Ymod(:,iModel);
    YYobs = Yobs;
    YYloo = Yloo;

    % deal with grandma
    I = find(isnan(YYmod)); YYloo(I) = []; YYobs(I) = []; YYmod(I) = [];
    I = find(YYmod <-9900); YYloo(I) = []; YYobs(I) = []; YYmod(I) = [];
    I = find(YYmod > 2000); YYloo(I) = []; YYobs(I) = []; YYmod(I) = [];

    % info bins
    Bmod = (min(YYmod)-Bw):Bw:(max(YYmod)+Bw);

    % benchmark information
    [Ib(iTarg,iBoot,iSite),Hy(iTarg,iBoot,iSite)] = info(YYobs,YYloo,Bobs,Bhat,1);

    % model information
    [Im(iTarg,iBoot,iSite,iModel),~] = info(YYobs,YYmod,Bobs,Bmod,1);

    % missing information
    MI(iTarg,iBoot,iSite,iModel) = (Ib(iTarg,iBoot,iSite)-Im(iTarg,iBoot,iSite,iModel))./Hy(iTarg,iBoot,iSite);
%    MI(iTarg,iBoot,iSite,iModel) = Im(iTarg,iBoot,iSite,iModel)./Hy(iTarg,iBoot,iSite);

    % screen report
    fprintf('finished calculating: boot %d/%d - targ %d/%d, - site %d/%d - model %d/%d \n',iBoot,Nboots,iTarg,Ntargs,iSite,Nsites,iModel,Nmodels);

   end % model
  end % site
 end % target
end % bootstrap

% average over bootstraps
MI = squeeze(mean(MI,2));

% save results
save('results/1lag_network_results.mat','-v7.3');

%% --- Plotting Tools --------------------------------------------

% plotting tools
figure(1);close(1);figure(1);
cc = plot(randn(7));

Mcolors(1 ,:) = cc(1).Color;
Mcolors(2 ,:) = cc(1).Color;
Mcolors(3 ,:) = cc(2).Color;
Mcolors(4 ,:) = cc(3).Color;
Mcolors(5 ,:) = cc(4).Color;
Mcolors(6 ,:) = cc(4).Color;
Mcolors(7 ,:) = cc(5).Color;
Mcolors(8 ,:) = cc(5).Color;
Mcolors(9 ,:) = cc(6).Color;
Mcolors(10,:) = cc(7).Color;
Mcolors(11,:) = cc(7).Color;
Mcolors(12,:) = cc(7).Color;
Mcolors(13,:) = 0.7*[1,1,1];
Mcolors(14,:) = 0.7*[1,1,1];

Scolors(1 ,:) = cc(1).Color;
Scolors(2 ,:) = cc(5).Color;
Scolors(3 ,:) = cc(1).Color;
Scolors(4 ,:) = cc(5).Color;
Scolors(5 ,:) = cc(1).Color;
Scolors(6 ,:) = cc(3).Color;
Scolors(7 ,:) = cc(4).Color;
Scolors(8 ,:) = cc(5).Color;
Scolors(9 ,:) = cc(4).Color;
Scolors(10,:) = cc(5).Color;
Scolors(11,:) = cc(3).Color;
Scolors(12,:) = cc(2).Color;

close(1);

%i = -1;
%for c = 1:7
% i = i+2;
% colors(i  ,:) = cc(c).Color;
% colors(i+1,:) = cc(c).Color;
%end
%close(1)

%colors(1 ,:) = [1/3,0  ,0];
%colors(2 ,:) = [2/3,0  ,0];
%colors(3 ,:) = [3/3,0  ,0];
%colors(4 ,:) = [0  ,1/3,0];
%colors(5 ,:) = [0  ,2/3,0];
%colors(6 ,:) = [0  ,3/3,0];
%colors(7 ,:) = [0  ,0  ,1/3];
%colors(8 ,:) = [0  ,0  ,2/3];
%colors(9 ,:) = [0  ,0  ,3/3];  
%colors(10,:) = [1  ,1  ,0];
%colors(11,:) = [1  ,0  ,1];
%colors(12,:) = [0  ,1  ,1];

% plot marker shapes
%markers = ['o','s','d','^','p','v','h'];
Mmarkers = ['o','p','o','o','o','p','o','p','o','o','p','v','o','p'];
%Smarkers = ['o','s','d','^','<','>','v','p','h','x','+','*'];
Smarkers = ['o','o','p','p','v','o','o','v','p','s','p','o'];

%% --- Make Plots ------------------------------------------------

make_plots














