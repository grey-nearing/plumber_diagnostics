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
fracs = logspace(-3,log10(0.25),20);
Nf = length(fracs);

% bootstraps
Nb = 10;

% number of mutual info bins
Nbins = 100;

% ann training parameters
trainParms.verbose       = 0;
trainParms.nodes         = 20;
trainParms.trainRatio    = 0.75;
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

% dimensions
[~,Nx,Nl] = size(data.forc);
Ny = size(data.flux,2);
Nd = size(data.date,2);

date = zeros(Nt,3,Ns);
forc = zeros(Nt,Nx,Nl,Ns);
flux = zeros(Nt,Ny,Ns);

% screen report
fprintf('Loading Data ... \n'); tic;
fprintf(repmat('.',[1,Ns]));
fprintf('\n');

% load training/test data from all sites
for s = 1:Ns
    fprintf('.');
    
    % load fluxnet data
    fname = strcat('../data/lagged_data/lagged_fluxnet_',...
        siteNames{s},'.mat');
    load(fname);
    data.flux = data.flux(:,1:3);
    data.forc = data.forc(:,:,:);
    
    % check dimensions
    assert(Nt <= length(data.date));
    
    % take random sample
    Idex = randsample(length(data.date),Nt);
    
    % store data in structures
    date(:,:,s)   = data.date(Idex,:);
    forc(:,:,:,s) = data.forc(Idex,:,:);
    flux(:,:,s)   = data.flux(Idex,:);
    
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
for s = 1:Ns
    
    % regressors
    temp = squeeze(forc(:,:,:,s));
    temp = permute(temp,[1,3,2]);
    Xsite(:,2:end,s) = reshape(temp,[Nt,Nl*Nx]);
    Xsite(:,1,s) = date(:,3,s);
    
    % regressands
    Ysite(:,:,s) = flux(:,:,s);
    
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

% % save progress
% fname = './results/progress_convergence.mat';
% load(fname);

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
            fprintf('Training %s ANN: Frac = %d/%d - Boot = %d/%d - Ndata = %d ...',...
                targNames{y},f,Nf,b,Nb,round(Nt*Ns*fracs(f))); tic;
            
            % train the network on segregated data
            ann{y,f,b} = trainANN(Xtrn,Ytrn,trainParms);
            
            % predict on segregated data
            Zall(:,y,f,b) = ann{y,f,b}(Xall');
            Ztst = Zall(Itst,y,f,b);
            Ztrn = Zall(Itrn,y,f,b);
            
            % calculate ann-all stats
            tstStats = calcStats(Ytst,Ztst,Bw); 
            tstMI(y,f,b) = tstStats.mi;
            tstRMSE(y,f,b) = tstStats.rmse;
            tstCORR(y,f,b) = tstStats.r;
            
            trnStats = calcStats(Ytrn,Ztrn,Bw); 
            trnMI(y,f,b) = trnStats.mi;
            trnRMSE(y,f,b) = trnStats.rmse;
            trnCORR(y,f,b) = trnStats.r;
            
            % screen report
            fprintf('. time = %f \n',toc);
            
            % screen report
            fprintf('Training MI = %f\n',trnMI(y,f,b));
            fprintf('Test MI = %f \n\n',tstMI(y,f,b));
            
        end % b-loop
        
        % save progress
        fname = './results/progress_convergence.mat';
        save(fname,'-v7.3');
         
    end % f-loop
end % y-loop

%% *** Make Plots *********************************************************
   
close all

% set up figure
iFig = iFig + 1;
figure(iFig); close(iFig); figure(iFig);
set(gcf,'color','w');
set(gcf,'position',[1925         403        1316         644]);

% plotting data: Y
Mtrn.mi = squeeze(mean(trnMI   ,3));
Strn.mi = squeeze(std( trnMI,[],3));
Mtst.mi = squeeze(mean(tstMI   ,3));
Stst.mi = squeeze(std( tstMI,[],3));

Mtrn.r2 = squeeze(mean(trnCORR   ,3));
Strn.r2 = squeeze(std( trnCORR,[],3));
Mtst.r2 = squeeze(mean(tstCORR   ,3));
Stst.r2 = squeeze(std( tstCORR,[],3));

Mtrn.se = squeeze(mean(trnRMSE   ,3));
Strn.se = squeeze(std( trnRMSE,[],3));
Mtst.se = squeeze(mean(tstRMSE   ,3));
Stst.se = squeeze(std( tstRMSE,[],3));

% plotting data: X
X = round(Nt*Ns*fracs);
X = X(1:size(Mtrn.mi,2));

% plot data
clear h legLabel
for y = 1:Ny
    
    subplot(2,Ny,y)
    
    % plot mi stats
    h(1) = semilogx(X,Mtrn.mi(y,:),'-o' ,'linewidth',2,...
        'markersize',8,'color',colors(2*y,:)); hold on;
    h(2) = semilogx(X,Mtst.mi(y,:),'--s','linewidth',2,...
        'markersize',8,'color',colors(2*y,:)); hold on;
    errorbar(X,Mtrn.mi(y,:),2*Strn.mi(y,:),'color',colors(2*y,:))
    errorbar(X,Mtst.mi(y,:),2*Stst.mi(y,:),'color',colors(2*y,:))
    legLabel(1) = strcat(targNames{y},{' - Training'});
    legLabel(2) = strcat(targNames{y},{' - Test'});
       
    % labels
    leg = legend(h(:),legLabel(:),'location','best');
    xlabel('Sample Size','fontsize',20);
    ylabel('Mutual Information Ratio','fontsize',20);
    title(strcat('Mutual Information Metric'),'fontsize',22);
    
    % aesthetics
    set(gca,'fontsize',16);
    grid on;
    
    subplot(2,Ny,Ny+y)
    
    % plot corr stats
    h(3) = plot(X,Mtrn.r2(y,:),'-o' ,'linewidth',2,...
        'markersize',8,'color',colors(2*y,:)); hold on;
    h(4) = plot(X,Mtst.r2(y,:),'--s','linewidth',2,...
        'markersize',8,'color',colors(2*y,:)); hold on;
    errorbar(X,Mtrn.r2(y,:),2*Strn.r2(y,:),'color',colors(2*y,:))
    errorbar(X,Mtst.r2(y,:),2*Stst.r2(y,:),'color',colors(2*y,:))
    legLabel(1) = strcat(targNames{y},{' - Training'});
    legLabel(2) = strcat(targNames{y},{' - Test'});
    
    % labels
    leg = legend(h(:),legLabel(:),'location','best');
    xlabel('Sample Size','fontsize',20);
    ylabel('Mutual Information Ratio','fontsize',20);
    title(strcat('Correlation Coefficient'),'fontsize',22);
    
    % aesthetics
    set(gca,'fontsize',16);
    grid on;
    
end % y-loop

% save figure
figure(iFig);
set(gcf,'PaperPositionMode','auto')
saveas(gcf,'./figures/FigA1 - convergence.png');

%% *** Save Results *******************************************************

% save progress
fname = './results/progress_convergence.mat';
save(fname,'-v7.3');

% save final statistics only
stats.tstStats.mi = tstMI;
stats.trnStats.mi = trnMI;
stats.tstStats.r2 = tstCORR;
stats.trnStats.r2 = trnCORR;
stats.tstStats.se = tstRMSE;
stats.trnStats.se = trnRMSE;
fname = './results/stats_convergence.mat';
save(fname,'stats','-v7.3');

%% *** END PROGRAM ********************************************************