clear all; close all; clc
restoredefaultpath; addpath('../matlab_tools');
iFig = 0;

% grab some plotting colors
colors = zeros(7,3); hc = plot(rand(7)); 
for i = 1:7; colors(i,:) = hc(i).Color; end
close all;

%% *** Experiment Setup ***************************************************

% target names
targNames = [{'Qe'},{'Qh'},{'NEE'}];

% site names
siteNames = [{'Amplero'},{'Blodgett'},{'Bugac'   },{'ElSaler2'},{'ElSaler' },...
             {'Espirra'},{'FortPeck'},{'Harvard' },{'Hesse'   },{'Howlandm'},...
             {'Howard' },{'Hyytiala'},{'Kruger'  },{'Loobos'  },{'Merbleue'},...
             {'Mopane' },{'Palang'  },{'Sylvania'},{'Tumba'   },{'UniMich' }];
Ns = length(siteNames);

% number of data to use from each site
Nt = 3e4;

% sample sizes
% fracs = logspace(-3,-2,4);
fracs = logspace(-3,log10(0.8),20);
Nf = length(fracs);

% bootstraps
Nb = 5;

% number of mutual info bins
Nbins = 100;

% ann training parameters
trainParms.verbose       = 0;
trainParms.nodes         = 20;
trainParms.trainRatio    = 0.65;
trainParms.max_fail      = 10;
trainParms.epochs        = 1e3;
trainParms.trainFcn      = 'trainscg';
trainParms.performFcn    = 'mse';

%% *** Load Data **********************************************************

% load fluxnet data
fname = strcat('../data/lagged_data/lagged_fluxnet_',siteNames{1},'.mat');
load(fname);
data.flux = data.flux(:,1:3);
data.forc = data.forc(:,:,:);

% load models data
fname = strcat('../data/lagged_data/lagged_models_',siteNames{1},'.mat');
load(fname);
model = model(:,10:12,:);

% dimensions
[~,Nx,Nl] = size(data.forc);
Ny = size(data.flux,2);
Nd = size(data.date,2);
Nm = size(model,3);

% init storage
Iwet = zeros(Nt,Ns)/0;
Idry = zeros(Nt,Ns)/0;
Nwet = zeros(Ns,1)/0;
Ndry = zeros(Ns,1)/0;

date = zeros(Nt,3,Ns);
forc = zeros(Nt,Nx,Nl,Ns);
flux = zeros(Nt,Ny,Ns);
modl = zeros(Nt,Ny,Nm,Ns);

% screen report
fprintf('Loading Data ... \n'); tic;
fprintf(repmat('.',[1,Ns]));
fprintf('\n');

% load training/test data from all sites
for s = 1:Ns
    fprintf('.');
    
    % load fluxnet data
    fname = strcat('../data/lagged_data/lagged_fluxnet_',siteNames{s},'.mat');
    load(fname);
    data.flux = data.flux(:,1:3);
    data.forc = data.forc(:,:,:);
    
    % load models data
    fname = strcat('../data/lagged_data/lagged_models_',siteNames{s},'.mat');
    load(fname);
    model = model(:,10:12,:);
    
    % check dimensions
    assert(length(model) == length(data.date));
    assert(Nt <= length(data.date));
    
    % take random sample
    Idex = randsample(length(data.date),Nt);
    
    % store data in structures
    date(:,:,s)   = data.date(Idex,:);
    forc(:,:,:,s) = data.forc(Idex,:,:);
    flux(:,:,s)   = data.flux(Idex,:);
    modl(:,:,:,s) = model(Idex,:,:);
    
end % s-loop

% error check
assert(isempty(find(isnan(flux),1)))
assert(isempty(find(isnan(forc),1)))
assert(isempty(find(isnan(date),1)))

% screen report
fprintf(' finished; time = %f \n',toc);

%% *** Pre-Processing *****************************************************

% mutual info bins
Bw = zeros(Ny,1);
for y = 1:Ny
    yy = flux(:,y,:);
    Bmin = min(yy(:))-1e-6;
    Bmax = max(yy(:))+1e-6;
    By = linspace(Bmin,Bmax,Nbins);
    Bw(y) = By(2) - By(1);
end

% collapse lagged dimension
Xsite = zeros(Nt,Nx*Nl+1,Ns   )/0;
Ysite = zeros(Nt,Ny     ,Ns   )/0;
Msite = zeros(Nt,Ny     ,Ns,Nm)/0;
for s = 1:Ns
    
    % regressors
    temp = squeeze(forc(:,:,:,s));
    temp = permute(temp,[1,3,2]);
    Xsite(:,2:end,s) = reshape(temp,[Nt,Nl*Nx]);
    Xsite(:,1,s) = date(:,3,s);
    
    % regressands
    Ysite(:,:,s) = flux(:,:,s);
    
    % predictions
    Msite(:,:,s,:) = modl(:,:,:,s);
    
end

% indexes of the precip input (to remove for wet/dry)
pdex = 1+(Nx-1)*Nl+1;
ddex = 1:Nx*Nl+1; ddex(pdex) = [];

% error check
assert(isempty(find(isnan(Xsite(:)),1)));
assert(isempty(find(isnan(Ysite(:)),1)));

% reshape data
Yall = reshape(permute(Ysite,[1,3,2]),[Nt*Ns,Ny]);
Xall = reshape(permute(Xsite,[1,3,2]),[Nt*Ns,Nx*Nl+1]);

%% *** Train Global ANN Models ********************************************

% init sorage - predictions
Zall = zeros(Nt*Ns,Ny,Nf,Nb)./0;

for y = 1:Ny            % loop through output dimensions
    for f = 1:Nf        % loop through data fractions
        for b = 1:Nb    % loop through bootstrap samples
            
            % training/test indexes
            Itrn = randsample(Nt*Ns,round(Nt*Ns*fracs(f)));
            Itst = 1:Nt*Ns; Itst(ismember(Itrn,Itst)) = [];
            
            % segragate trainig test data
            Xtrn = Xall(Itrn,:); Ytrn = Yall(Itrn,y);
            Xtst = Xall(Itst,:); Ytst = Yall(Itst,y);
            
            % skip if there is no data of this type for this site
            if all(isnan(Ytrn)); continue; end
            if all(isnan(Ytst)); continue; end
            
            % remove grandmas
            assert(isempty(find(isnan(Xtrn(:)),1)));
            assert(isempty(find(isnan(Ytrn(:)),1)));
            assert(isempty(find(isnan(Xtst(:)),1)));
            assert(isempty(find(isnan(Ytst(:)),1)));
            
            % screen report
            fprintf('Training ANN: Target = %s - Fraction = %d/%d - Boot = %d/%d - Ndata = %d ...',...
                targNames{y},f,Nf,b,Nb,round(Nt*Ns*fracs(f))); tic;
            
            % train the network on segregated data
            ann{y,f,b} = trainANN(Xtrn,Ytrn,trainParms);
            
            % predict on segregated data
            Zall(:,y,f,b) = ann{y,f}(Xall');
            Ztst = Zall(Itst,y,f,b);
            Ztrn = Zall(Itrn,y,f,b);
            
            % calculate ann-all stats
            temp = calcStats(Ytst,Ztst,Bw); tstStats(y,f,b) = temp.mi;
            temp = calcStats(Ytrn,Ztrn,Bw); trnStats(y,f,b) = temp.mi;
            
            % screen report
            fprintf('. finished; time = %f \n',toc);
            
            % screen report
            fprintf('Training MI = %f\n',trnStats(y,f,b));
            fprintf('Test MI = %f \n\n',tstStats(y,f,b));
            
        end % b-loop
    end % f-loop
    
    % save progress
    fname = './results/progress_convergence.mat';
    save(fname,'-v7.3');
    
end % y-loop

%% *** Make Plots *********************************************************
   
% set up figure
iFig = iFig + 1;
figure(iFig); close(iFig); figure(iFig);
set(gcf,'color','w');

% plotting data: X
X = round(Nt*Ns*fracs);
Mtrn = squeeze(mean(trnStats   ,3));
Strn = squeeze(std( trnStats,[],3));
Mtst = squeeze(mean(tstStats   ,3));
Stst = squeeze(std( tstStats,[],3));
 
% plot data
clear h legLabel
for y = 1:Ny
    h(1,y) = plot(X,Mtrn(y,:),'-o' ,'linewidth',2,'markersize',8,'color',colors(2*y,:)); hold on;
    h(2,y) = plot(X,Mtst(y,:),'--s','linewidth',2,'markersize',8,'color',colors(2*y,:)); hold on;
    legLabel(1,y) = strcat(targNames{y},{' - Training'});
    legLabel(2,y) = strcat(targNames{y},{' - Test'});
end % y-loop

% labels
leg = legend(h(:),legLabel(:),'location','ne');
xlabel('Sample Size','fontsize',20);
ylabel('Mutual Information Ratio','fontsize',20);
title(strcat('Convergence by Sample Size'),'fontsize',22);

% aesthetics
set(gca,'fontsize',16);
grid on;

% save figure
figure(iFig);
set(gcf,'PaperPositionMode','auto')
saveas(gcf,'./figures/FigA1 - convergence.png');

%% *** Save Results *******************************************************

% save progress
fname = './results/progress_convergence.mat';
save(fname,'-v7.3');

% save final statistics only
fname = './results/stats_convergence.mat';
save(fname,'stats','-v7.3');

%% *** END PROGRAM ********************************************************