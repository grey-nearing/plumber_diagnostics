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
              {'Manabe_Bucket.2 '}
              {'Penman_Monteith.1'}];

% site names
siteNames =  [...
             {'Amplero'}
             {'Blodgett'}
             {'Bugac'}
             {'ElSaler2'}
             {'ElSaler'}
             {'Espirra'}
             {'FortPeck'}
             {'Harvard'}
             {'Hesse'}
             {'Howlandm'}
             {'Howard'}
             {'Hyytiala'}
             {'Kruger'}
             {'Loobos'}
             {'Merbleue'}
             {'Mopane'}
             {'Palang'}
             {'Sylvania'}
             {'Tumba'}
             {'UniMich'}];


%% --- Load Transfer Results -------------------------------------

load('results/main_dpn_results_20.mat');
T = T./H;

Tp = squeeze(        T(:,6:end,:,1    ,:));
Hp = squeeze(        H(:,6:end,:,1    ,:));
Ta = squeeze(        T(:,6:end,:,2    ,:));
Ha = squeeze(        H(:,6:end,:,2    ,:));
Tm = squeeze(permute(T(:,6:end,:,3:end,:),[1,2,3,5,4]));
Hm = squeeze(permute(H(:,6:end,:,3:end,:),[1,2,3,5,4]));

Ns = size(Tm,4);
Nm = size(Tm,5);
Dx = size(Tm,1);
Dy = size(Tm,2);

%% --- Load Benchmarking Results ---------------------------------

% nubmer of bootstraps
Nboots = 1;

% benchmark targets
targNames = [{'Qe'},{'Qh'},{'NEE'}];
Ntargs = length(targNames);

% loop through experiments
for iBoot = 1:Nboots
 for iTarg = 1:Ntargs
  for iSite = 1:Ns

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
   Ymod = squeeze(model(:,8+iTarg,:));

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
   for iModel = 1:Nm

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
    fprintf('finished calculating: boot %d/%d - targ %d/%d, - site %d/%d - model %d/%d \n',iBoot,Nboots,iTarg,Ntargs,iSite,Ns,iModel,Nm);

   end % model
  end % site
 end % target
end % bootstrap

% average over bootstraps
%MI = squeeze(mean(MI,2));
MI = squeeze(MI);
TD = repmat(squeeze(Tp(:,1:3,1,:)),[1,1,1,Nm])-squeeze(Tm(:,1:3,1,:,:));

TDD = TD;
MID = MI;

S = [1,2,3,4,5,7,11,12,13,14,16,17,18];
TD = TDD(:,:,S,:);
MI = MID(:,S,:);
siteNames = siteNames(S);
Ns = size(MI,2);

% save results
save('results/tradeoff_results.mat','-v7.3');

%% --- Make Plots ------------------------------------------------

te_mi_diagrams
clusters_model_groups













