clear all; close all; clc

iType = 3;
iTarg = 1;
iBoot = 1;
iSite = 1;

% function benchmark_bootstrap(iType,iTarg,iSite,iBoot)

% master routines
addpath('../matlab_tools');

% time the fucntion call
masterTime = tic;

% maximum number of data
Dmax = 1e2*2^11;

% training fractions
Ntrn = 1e2*2.^9;
%Ntrn = 1e2;

% training configuration
Nepochs = 10000;
trnfctn = 'trainscg';
%trnfctn = 'traincgb';
%trnfctn = 'trainbr';

% type names
typeNames = [{'all'},{'dry'},{'wet'}];%

% target names
targNames = [{'Qe'},{'Qh'},{'NEE'}];
targs = [1,2,3];
Ntargs = length(targs);

% site names
siteNames = [{'Amplero'},{'Blodgett'},{'Bugac'   },{'ElSaler2'},{'ElSaler' },...
    {'Espirra'},{'FortPeck'},{'Harvard' },{'Hesse'   },{'Howlandm'},...
    {'Howard' },{'Hyytiala'},{'Kruger'  },{'Loobos'  },{'Merbleue'},...
    {'Mopane' },{'Palang'  },{'Sylvania'},{'Tumba'   },{'UniMich' }];
Nsites = length(siteNames);

% load training/test data from all sites
fprintf('Loading data: Type = %s - Site =  ',typeNames{iType}); loadTime = tic; 
for s = 1:Nsites
    fprintf('%d, ',s);
    
    fname = strcat('../data/lagged_data/',typeNames{iType},'_',siteNames{s},'.mat');
    load(fname);
    
    Ndat = size(data.date,1);
    if iType==3; Nuse = Ndat;
    else Nuse = floor(Dmax/Nsites); end
    
    Idx{s} = randperm(Ndat,Nuse);
    
    if (iSite > 1 && s == 1) || (iSite == 1 && s == 2)
        date = data.date(Idx{s},:);
        forc = data.forc(Idx{s},:,:);
        prog = data.flux(Idx{s},:);
    elseif iSite ~= s
        date = cat(1,date,data.date(Idx{s},:));
        forc = cat(1,forc,data.forc(Idx{s},:,:));
        prog = cat(1,prog,data.flux(Idx{s},:));
    elseif iSite == s
        sdate = data.date;
        sforc = data.forc;
        sprog = data.flux;
    else
        error('check the site loading routine');
    end
    
end; fprintf('/n');

% make sure forcings are not missing
assert(isempty(find(isnan(forc(:)))));
assert(isempty(find(isnan(date(:)))));
assert(isempty(find(isnan(prog(:,targs(iTarg))))));
assert(isempty(find(forc(:)<-990)));
assert(isempty(find(date(:)<-990)));
assert(isempty(find(prog(:,targs(iTarg))<-990)));
assert(isempty(find(isnan(sforc(:)))));
assert(isempty(find(isnan(sdate(:)))));
assert(isempty(find(isnan(sprog(:,targs(iTarg))))));
assert(isempty(find(sforc(:)<-990)));
assert(isempty(find(sdate(:)<-990)));
assert(isempty(find(sprog(:,targs(iTarg))<-990)));

% number of target types
assert(Ntargs <= size(prog,2));

% screen report
t = toc(loadTime); fprintf('finished loading data: t = %f \n',t);

% loo training
fprintf('LOO Training: Type = %d - Targ = %d - Boot = %d - Site = %d \n',iType,iTarg,iBoot,iSite); trainTime = tic; % screen report
X = reshape(forc,[size(forc,1),size(forc,2)*size(forc,3)]);          % extract loo inputs
X = cat(2,date(:,2:3),X);                                            % add dates to inputs
Y = prog(:,targs(iTarg));                                            % extract loo training targets
I = find(isnan(Y)); Y(I) = []; X(I,:) = [];                          % remove missing data
assert(size(X,1)==size(Y,1)); N = size(X,1);                         % number of loo training data available
assert(N>1.5*Ntrn);                                                  % make sure we have enough loo data
I = randperm(N,Ntrn);                                                % select training points
J = 1:N; J(ismember(J,I)) = [];                                      % select test points
[model,results.train] = ann_experiment(X,Y,I,J,trnfctn,Nepochs);     % train loo model
results.train.Isites = Idx;                                          % add site indexes to output stream
t = toc(trainTime); fprintf('finished LOO training: t = %f \n\n',t); % screen report

% test loo on site data
fprintf('LOO Prediction: Type = %d - Targ = %d - Boot = %d - Site = %d \n',iType,iTarg,iBoot,iSite); predTime = tic; % screen report
X = reshape(sforc,[size(sforc,1),size(sforc,2)*size(sforc,3)]);      % extract site inputs
X = cat(2,sdate(:,2:3),X);                                           % add dates to site inputs
Y = sprog(:,targs(iTarg));                                           % extract site targets
I = find(isnan(Y)); Y(I) = []; X(I,:) = [];                          % remove missing data
assert(size(X,1)==size(Y,1)); N = size(X,1);                         % number of site data
I = 1:N; J = [];                                                     % all are test points
results.test = pred_ann(X,Y,model,I,J);                              % predict and stats
t = toc(predTime); fprintf('finished LOO prediction: t = %f \n\n',t);% screen report

% save results
fprintf('Saving LOO results ... '); saveTime = tic; % screen report
results.model = model;
fname = strcat('results/LOO_',typeNames{iType},'_',targNames{iTarg},'_',siteNames{iSite},'_',num2str(iBoot),'_',num2str(Nepochs),'_',trnfctn,'.mat');
save(fname,'results');
t = toc(saveTime); fprintf('finished saving LOO results: t = %f \n',t);
clear results

% site training/test
fprintf('Site Training: Type = %d - Targ = %d - Boot = %d - Site = %d \n',iType,iTarg,iBoot,iSite); trainTime = tic; % screen report
X = reshape(sforc,[size(sforc,1),size(sforc,2)*size(sforc,3)]);      % extract site inputs
X = cat(2,sdate(:,2:3),X);                                           % add dates to inputs
Y = sprog(:,targs(iTarg));                                           % extract site training targets
I = find(isnan(Y)); Y(I) = []; X(I,:) = [];                          % remove missing data
assert(size(X,1)==size(Y,1)); N = size(X,1);                         % number of site training data available
I = 1:N; J = [];                                                     % select training points
[model,results] = ann_experiment(X,Y,I,J,trnfctn,Nepochs);           % train site model
t = toc(trainTime); fprintf('finished site training: t = %f \n\n',t);% screen report

% save results
fprintf('Saving site results ... '); saveTime = tic; % screen report
results.model = model;
fname = strcat('results/Site_',typeNames{iType},'_',targNames{iTarg},'_',siteNames{iSite},'_',num2str(iBoot),'_',num2str(Nepochs),'_',trnfctn,'.mat');
save(fname,'results');
t = toc(saveTime); fprintf('finished saving site results: t = %f \n',t);

%% **** End Experiment ****************************************

% report total time
t = toc(masterTime); fprintf('\nAll Finished - time = %f ',t);
fprintf('\n\n\n----------------------------------------------------\n\n\n');

% report that this matlab job is finished
fname = strcat('reports/benchmark_bootstrap_end_report_',num2str(iType),'_',num2str(iTarg),'_',num2str(iSite),'_',num2str(iBoot),'.end');
fid = fopen(fname,'w'); fprintf(fid,'done'); fclose(fid);

%end
%end
%end
