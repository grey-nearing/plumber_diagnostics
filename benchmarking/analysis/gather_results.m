clear all
close all
clc
iFig = 0;
addpath('../../matlab_tools')

% get colors
figure(100);
h = plot(randn(20));
colors = get(h,'Color');
close(100);

% information resolutions
resolutions = [0.01,0.02,0.05,0.10];
resNames = [{'1% info resolution'},{'2% info resolution'},{'5% info resolution'},{'10% info resolution'}];
Nres = length(resolutions);

% site names
sites = [{'Amplero'},{'Blodgett'},{'Bugac'   },{'ElSaler2'},{'ElSaler' },...
         {'Espirra'},{'FortPeck'},{'Harvard' },{'Hesse'   },{'Howlandm'},...
         {'Howard' },{'Hyytiala'},{'Kruger'  },{'Loobos'  },{'Merbleue'},...
         {'Mopane' },{'Palang'  },{'Sylvania'},{'Tumba'   },{'UniMich' }];
Nsites = length(sites);

% Models
modelNames = [{'CABLE 2.0'},{'CABLE 2.0 SLI'},{'CHTESSEL'},{'COLASSiB 2.0'},{'ISBA SURFEX 7.3 3l'},{'ISBA SURFEX 7.3 dif'},{'JULES 3.1'},{'JULES 3.1 (altP)'},{'Manabe Bucket'},{'Mosaic'},{'NOAH 2.7.1'},{'Noah 3.2'},{'NOAH 3.3'},{'Penman Monteith'}];
%          {'ORCHIDEE.trunk_r1401'}
%          {'SUMMA.1.0.exp.01.000'}
%          {'SUMMA.1.0.exp.01.001'}
%          {'SUMMA.1.0.exp.01.002'}
%          {'SUMMA.1.0.exp.01.003'}
%          {'SUMMA.1.0.exp.01.004'}
%          {'SUMMA.1.0.exp.01.005'}
%          {'SUMMA.1.0.exp.01.006'}
%          {'SUMMA.1.0.exp.01.007'}
%          {'SUMMA.1.0.exp.01.008'}
%          {'SUMMA.1.0.exp.01.009'}
Nmodels = length(modelNames);

% stuff we need for the file name
Nepochs = 10000;
trnfctn = 'trainscg';

Du = 5;

% experiment targets
targNames = [{'Qe'},{'Qh'},{'NEE'}];
Ntargs = length(targNames);

% number of boodstraps
Nboots = 1;

% loop through experiments
for iTarg = 1:Ntargs
 for iBoot = 1:Nboots

  % screen report 
  tboot = tic;

  % init storage
  Yloo = []; 
  Yobs = [];
  Ysit = []; 
  Nsite = zeros(Nsites,1);

  for iSite = 1:Nsites

%% --------------------------
   % file name 
   LOOfname = strcat('../results/LOO_all_',targNames{iTarg},'_',sites{iSite},'_',num2str(iBoot),'_',num2str(Nepochs),'_',trnfctn,'.mat');
   SITEfname = strcat('../results/Site_all_',targNames{iTarg},'_',sites{iSite},'_',num2str(iBoot),'_',num2str(Nepochs),'_',trnfctn,'.mat');

   % load file
%   try
    load(LOOfname);  LOO  = results; clear results;
    load(SITEfname); SITE = results; clear results;
%   catch
%    continue
%   end 

   % make sure there was no stupid mistake
   assert(max(abs(LOO.test.Yobs(:)-SITE.Yobs(:)))==0);

   % store all data
   Yloo = [Yloo(:);LOO.test.Yhat(:)];;
   Yobs = [Yobs(:);LOO.test.Yobs(:)];
   assert(length(Yloo) == length(Yobs));
   Nsite(iSite) = length(LOO.test.Yobs(:));

   % store all data
   Ysit = [Ysit(:);SITE.Yhat(:)];;
   assert(length(Ysit) == length(Yobs));

  end % iSite

%% -------------------------- 
  % load models
  Ymod = zeros(length(Ysit),Nmodels)/0;
  
  edex = 0;
  for iSite = 1:Nsites

    % file name
    fname = strcat('../../data/lagged_data/models_',sites{iSite},'.mat');
    load(fname);
    if Nsite(iSite) > 0
     sdex = edex+1; edex = sdex+Nsite(iSite)-1;
     Ymod(sdex:edex,:) = squeeze(model(:,3+Du+iTarg,:)); 
    end

  end % iSite

%% --------------------------
  % remove missing values
  I = find(Yobs<=-990); Yobs(I) = []; Yloo(I) = []; Ysit(I) = []; Ymod(I,:) = []; 
  I = find(Yloo<=-990); Yobs(I) = []; Yloo(I) = []; Ysit(I) = []; Ymod(I,:) = []; 
  I = find(Ysit<=-990); Yobs(I) = []; Yloo(I) = []; Ysit(I) = []; Ymod(I,:) = []; 

  I = find(isnan(Yobs)); Yobs(I) = []; Yloo(I) = []; Ysit(I) = []; Ymod(I,:) = []; 
  I = find(isnan(Yloo)); Yobs(I) = []; Yloo(I) = []; Ysit(I) = []; Ymod(I,:) = []; 
  I = find(isnan(Ysit)); Yobs(I) = []; Yloo(I) = []; Ysit(I) = []; Ymod(I,:) = []; 

%% --------------------------
  % calculate loo statistics
  if isempty(Yloo)
   INFO(1,iTarg,iBoot,1:Nres) = 0./0;
  else
   for r = 1:Nres; tloo = tic;
    Bw = (max(Yobs)-min(Yobs))*resolutions(r);
[min(Yobs),max(Yobs),Bw]
    Bhat = (min(Yloo)-Bw):Bw:max(Yloo+Bw);
    Bobs = (min(Yobs)-Bw):Bw:max(Yobs+Bw);
    [INFO(1,iTarg,iBoot,r),H(iTarg,iBoot,r)] = info(Yloo,Yobs,Bhat,Bobs,2);
    t = toc(tloo); fprintf('Finished loo stats - iTarg = %d/%d - iBoot = %d/%d - iRes = %d/%d - time = %f\n',iTarg,Ntargs,iBoot,Nboots,r,Nres,t);
   end
  end

  % calculate site statistics
  if isempty(Ysit)
   INFO(2,iTarg,iBoot,1:Nres) = 0./0;
  else
   for r = 1:Nres; tsite = tic;
    Bw = (max(Yobs)-min(Yobs))*resolutions(r);
    Bhat = (min(Ysit)-Bw):Bw:max(Ysit+Bw);
    Bobs = (min(Yobs)-Bw):Bw:max(Yobs+Bw);
    INFO(2,iTarg,iBoot,r) = info(Ysit,Yobs,Bhat,Bobs,2);
    t = toc(tsite); fprintf('Finished site stats - iTarg = %d/%d - iBoot = %d/%d - iRes = %d/%d - time = %f\n',iTarg,Ntargs,iBoot,Nboots,r,Nres,t);
   end
  end

  % calculate model statistics
  if iBoot == 1
   for iMod = 1:Nmodels
    if isempty(Ysit)
     INFO(2+iMod,iTarg,iBoot,1:Nres) = 0./0;
    elseif max(abs(diff(Ymod(:,iMod)))) < 1e-6
     INFO(2+iMod,iTarg,iBoot,1:Nres) = 0./0;
    else

     I = find(Ymod(:,iMod)>=-990); %[iMod,length(I),length(Yobs)]
     if length(I)<0.9*length(Yobs)
      INFO(2+iMod,iTarg,iBoot,1:Nres) = 0./0;
      continue; 
     end;

     I = find(~isnan(Ymod(:,iMod))); %[iMod,length(I),length(Yobs)]
     if length(I)<0.9*length(Yobs)
      INFO(2+iMod,iTarg,iBoot,1:Nres) = 0./0;
      continue
     end

     for r = 1:Nres; tmod = tic;
      Bw = (max(Yobs)-min(Yobs))*resolutions(r);
      Bhat = (min(Ymod(I,iMod))-Bw):Bw:max(Ymod(I,iMod)+Bw);
      Bobs = (min(Yobs)-Bw):Bw:max(Yobs+Bw);
      INFO(2+iMod,iTarg,iBoot,r) = info(Ymod(I,iMod),Yobs(I),Bhat,Bobs,2);
      t = toc(tmod) ; fprintf('Finished model stats - iTarg = %d/%d - iMod = %d/%d - iRes = %d/%d - time = %f\n',iTarg,Ntargs,iMod,Nmodels,r,Nres,t);
     end
    end
   end
  end

  % screen report
  t = toc(tboot); fprintf('Finished: iBoot = %d - iTarg = %d - time = %f\n',iBoot,iTarg,t);
  % screen report
  t = toc(tboot); fprintf('Finished: iBoot = %d - iTarg = %d - time = %f\n',iBoot,iTarg,t);

 end % iBoot
end % iTarg

% save results
save('loo_site_results3.mat');

%% ---------------------------------------

make_plots


return

