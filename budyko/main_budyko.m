clear all; close all; clc
iFig = 0;

% grab some plotting colors
h = plot(randn(7));
for i = 1:7; colors(i,:) = h(i).Color; end
close all;

%% *** Experiment Setup ***************************************************

% site names
sites = [{'Amplero'},{'Blodgett'},{'Bugac'   },{'ElSaler2'},{'ElSaler' },...
    {'Espirra'},{'FortPeck'},{'Harvard' },{'Hesse'   },{'Howlandm'},...
    {'Howard' },{'Hyytiala'},{'Kruger'  },{'Loobos'  },{'Merbleue'},...
    {'Mopane' },{'Palang'  },{'Sylvania'},{'Tumba'   },{'UniMich' }];
Nsites = length(sites);

% Models
modelNames = [ ...
    {'CABLE 2.0'}
    {'CABLE 2.0 SLI'}
    {'CHTESSEL'}
    {'COLASSiB 2.0'}
    {'ISBA SURFEX 7.3 3l'}
    {'ISBA SURFEX 7.3 dif'}
    {'JULES 3.1'}
    {'JULES 3.1 (altP)'}
    {'Manabe Bucket'}
    {'Mosaic'}
    {'Noah 2.7.1'}
    {'Noah 3.2'}
    {'Noah 3.3'}
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
    {'Penman Monteith'}];
Nmodels = length(modelNames);

%% *** Load Data **********************************************************

% % screen report
% fprintf('Loading Data ... \n'); tic;
% fprintf(repmat('.',[1,Nsites]));
% fprintf('\n');
% 
% % loop sites
% for s = 1:Nsites
%     fprintf('.');
%     
%     % load pals data
%     fname = strcat('../data/pals_data/extracted/',sites{s},'.txt');
%     pals(:,:,s) = load(fname);
%     
%     % load  data
%     fname = strcat('../data/model_data/extracted/',sites{s},'.mat');
%     data = load(fname);
%     model(:,:,:,s) = data.model;
%     
% end % sites
% 
% % screen report
% fprintf(' finished; time = %f \n',toc);
% 
% all_data.pals = pals;
% all_data.model = model;
% save('all_data.mat','all_data','-v7.3');

load('all_data.mat');
pals = all_data.pals;
model = all_data.model;

%% *** Calculate Indices **************************************************

% screen report
fprintf('Calculating Indices ... \n'); tic;
fprintf(repmat('.',[1,Nsites]));
fprintf('\n');

% init storage
PE  = zeros(1,Nsites)./0;
PP  = zeros(1,Nsites)./0;
DI  = zeros(1,Nsites)./0;
TP  = zeros(1,Nsites)./0;
EE  = zeros(1,Nsites)./0;
EF  = zeros(1,Nsites)./0;
EEm = zeros(Nsites,Nmodels)./0;
EFm = zeros(Nsites,Nmodels)./0;

% aggregate to annual means
for s = 1:Nsites
    fprintf('.');
    
    % find missing values in PALS data
    Iy = find(all(~isnan(squeeze(pals(:,:,s))')))';
    if isempty(Iy); continue; end
    
    % find missing vales in model
    Im1 = find(all(~isnan(squeeze(model(Iy     ,9,:,s))')));
    Im2 = find(all(~isnan(squeeze(model(Iy(Im1),10,:,s))')));
    Iy = Iy(Im1(Im2));
    if length(Iy) < 2e4; continue; end
    
    % dryness index
    PE(s) = mean(model(Iy,10,end,s))/2454000 * 60*30 *length(Iy);
    PP(s) = sum(pals(Iy,9,s));% * 60*30;% *length(Iy);
    DI(s) = PE(s)/PP(s);
    
    % evaporative fraction
    for m = 1:Nmodels
        EEm(s,m) = nanmean(model(Iy,10,m,s))/2454000 * 60*30 *length(Iy);
        EFm(s,m) = EEm(s,m) / PP(s);
    end % models
    EE(s) = mean(pals(Iy,10,s))/2454000 * 60*30 *length(Iy);
    EF(s) = EE(s)/PP(s);
    
    % Turc-Pike
    v = 2;
    TP(s) = (1+DI(s)^-v)^(-1/v);% * PP(y,s);
    
    % regression inputs
    Ta(s) = mean(pals(Iy,5,s));
    Rn(s) = mean(pals(Iy,6,s));
        
end % sites

% screen report
fprintf(' finished; time = %f \n',toc);

%% *** Perform Regression *************************************************

% screen report
fprintf('Training Regressions ... \n'); tic;
fprintf(repmat('.',[1,Nsites]));
fprintf('\n');

clear reg LR
for s = 1:Nsites
    fprintf('.');
    
    % grab loo training/testing data
    S = 1:Nsites; S(S == s) = [];
    di = DI(S); pp = PP(S); pe = PE(S); ta = Ta(S); rn = Rn(S); ef = EF(S); tp = TP(S);
    Xtrain = [tp(:),ta(:),pp(:),rn(:),ones(size(di(:)))];
    Ytrain = ef(:);% - tp(:);
    Xtest = [TP(s),Ta(s),PP(s),Rn(s),1];
    Ytest = EF(s);% - TP(s);
    
    % remove missing values
    Ix = find(any(isnan(Xtrain)));
    Xtrain(Ix,:) = [];
    Ytrain(Ix) = [];
    Iy = find(isnan(Ytrain));
    Xtrain(Iy,:) = [];
    Ytrain(Iy) = [];
    
    % train model
    reg(:,s) = Xtrain\Ytrain;
    
    % predict with trained model
    Ix = find(all(~isnan(Xtest)));
    LR(s) = Xtest * reg(:,s);% + TP(s);
    
end

% screen report
fprintf(' finished; time = %f \n',toc);

%% *** Calculate Bootstrap Significance ***********************************

% total number of data
assert(all(~isnan(DI)));

% number of bootstrap samples
Nboot = 1e4;

% storage
bu = zeros(Nboot,Nmodels)/0;

% loop through models
for m = 1:Nmodels
    
    % true mean
    mu(m) = abs(mean(LR) - mean(EFm(:,m)));
    
    % total data vector
    comb = cat(1,LR',EFm(:,m));
    
    % bootstrap means
    for b = 1:Nboot
        
        % sample with replacement
        S = datasample(1:length(comb),length(comb));
        S1 = comb(S(1:length(comb)/2));
        S2 = comb(S(length(comb)/2+1:end));
        
        % bootstrap mean
        bu(b,m) = abs(mean(S1) - mean(S2));
        
    end % bootstrap
    
    % calculate significance
    pValue(2+m) = length(find(bu(:,m)>mu(m)))/Nboot;
    
end % models

%% *** Calculate RMSEs ****************************************************

% total number of data
I = find(~isnan(DI));
Ni = length(I);

% bootstrap samples
Nboot = 1000;
for b = 1:Nboot
    
    % sample with replcement
    S = datasample(1:Ni,Ni);
    
    % calculte RMSEs
    RMSE.TP(b) = rmse(TP(I(S))',EF(I(S))');
    RMSE.LR(b) = rmse(LR(I(S))',EF(I(S))');
    for m = 1:Nmodels
        efm = squeeze(EFm(:,m));
        RMSE.model(m,b) = rmse(efm(I(S)),EF(I(S))');
    end
    
    % calculte MAEs
    MAE.TP(b) = mae(TP(I(S))',EF(I(S))');
    MAE.LR(b) = mae(LR(I(S))',EF(I(S))');
    for m = 1:Nmodels
        efm = squeeze(EFm(:,m));
        MAE.model(m,b) = mae(efm(I(S)),EF(I(S))');
    end
    
end % Nboot

RMSEs = cat(1,[RMSE.TP;RMSE.LR],RMSE.model(1:end,:));
[~,rank] = sort(mean(RMSEs'));
RMSEs = RMSEs(rank,:);

MAEs = cat(1,[MAE.TP;MAE.LR],MAE.model(1:end,:));
[~,rank] = sort(mean(MAEs'));
MAEs = MAEs(rank,:);

%% *** Plot Budyko Curve **************************************************

% Turc-Pike curve
Xtp = 0:0.1:10.5;
Ytp = (1+Xtp.^-v).^(-1/v);

% initialize figure
iFig = iFig+1; figure(iFig); close(iFig); figure(iFig);
set(gcf,'color','w','position',[1000,894,1050,425]);

% plot turc-pike curve
h1 = plot(Xtp,Ytp,'-k','linewidth',0.5); hold on;
h2 = plot(DI,EF,'s','color',colors(1,:),'markersize',8,'markerfacecolor',colors(1,:));
clear h

% plot model data
for s = 1:Nsites
    h(s,:) = plot(DI(s),squeeze(EFm(s,1:end-1)),'.','color',colors(2,:));
end

% plot regression data
h3 = plot(DI,LR,'o','color',colors(4,:),'markersize',8,'markerfacecolor',colors(4,:));

% plot budyko curve
h0 = plot([0,1],[0,1],'k--');
plot([1,10.5],[1,1],'k--');

% labels
xlabel('Dryness Index','fontsize',16);
ylabel('Evaporative Fraction','fontsize',16);
legend([h1,h2(1),h(1,1),h3(1)],'Turc-Pike','PALS Data','Land Models',...
    'Linear Regression','location','nw');%modelNames(1:end-1)');

% aesthetics
set(gca,'ylim',[0,1.5]);
set(gca,'xlim',[0,10.5]);

% save Budybo Curve figure
figure(1);
set(gcf,'PaperPositionMode','auto')
saveas(gcf,'figures/BudykoCurve.png');

%% *** Plot RMSE Statistics ***********************************************

xnames = [{'Turc-Pike'};{'Benchmark Reg.'};modelNames(1:end)];

iFig = iFig+1; figure(iFig); close(iFig); figure(iFig);
set(gcf,'color','w','position',[440   324   624   474])

% boxplot(RMSEs');
boxplot(MAEs');
set(gca,'xticklabel',xnames(rank));
h = gca;
h.XTickLabelRotation = 60;
grid on;
% ylabel('RMSE of Evap. Frac.','fontsize',20);
ylabel('MAE of Evap. Frac.','fontsize',20);
set(gca,'ylim',[0,0.5]);
set(gca,'fontsize',16)

% save MAE figure
figure(iFig);
set(gcf,'PaperPositionMode','auto')
saveas(gcf,'figures/Figure 3 - BudykoMAE.png');

%% *** End Program ********************************************************



