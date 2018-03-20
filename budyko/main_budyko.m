clear all
close all
clc

% ----------------------------------------------------------------------------

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
              {'NOAH 2.7.1'}
              {'Noah 3.2'}
              {'NOAH 3.3'}
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

% ----------------------------------------------------------------------------

%% load data
%for s = 1:Nsites
%
% % load pals data
% fname = strcat('../data/pals_data/extracted/',sites{s},'.txt');
% pals(:,:,s) = load(fname);
%
% % load  data
% fname = strcat('../data/model_data/extracted/',sites{s},'.mat');
% data = load(fname);
% model(:,:,:,s) = data.model;
% 
% [s/Nsites]
%
%end % sites
%
%all_data.pals = pals;
%all_data.model = model;
%save('all_data.mat','all_data','-v7.3');

load('all_data.mat');
pals = all_data.pals;
model = all_data.model;

% ----------------------------------------------------------------------------

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

 % find missing values in PALS data
 Iy = find(all(~isnan(squeeze(pals(:,:,s))')))';
 if isempty(Iy); continue; end;

 % find missing vales in model
 Im1 = find(all(~isnan(squeeze(model(Iy    ,8,:,s))')));
 Im2 = find(all(~isnan(squeeze(model(Iy(Im1),9,:,s))')));
 Iy = Iy(Im1(Im2));
 if length(Iy) < 2e4; continue; end;

 % dryness index
 PE(s) = mean(model(Iy,9,end,s))/2454000 * 60*30 *length(Iy);
 PP(s) = sum(pals(Iy,8,s));% * 60*30;% *length(Iy);
 DI(s) = PE(s)/PP(s);

 % evaporative fraction
 for m = 1:Nmodels
  EEm(s,m) = nanmean(model(Iy,9,m,s))/2454000 * 60*30 *length(Iy);
  EFm(s,m) = EEm(s,m) / PP(s);
 end % models 
 EE(s) = mean(pals(Iy,9,s))/2454000 * 60*30 *length(Iy);
 EF(s) = EE(s)/PP(s);

 % Turc-Pike
 v = 2;
 TP(s) = (1+DI(s)^-v)^(-1/v);% * PP(y,s);

 % regression inputs
 Ta(s) = mean(pals(Iy,5,s));
 Rn(s) = mean(pals(Iy,6,s));

 [s,length(Iy),length(find(~isnan(model(Iy,9,m,s))))]

end % sites

% ----------------------------------------------------------------------------

clear reg LR
for s = 1:Nsites

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

 % screen report
 fprintf('trained/predicted loo model at site %d of %d \n',s,Nsites);

end

% ----------------------------------------------------------------------------

figure(1);
h = plot(randn(7));
for i = 1:7; colors(i,:) = h(i).Color; end;
close all;

Xtp = 0:0.1:10.5;
Ytp = (1+Xtp.^-v).^(-1/v);

figure(1); close(1); figure(1);
set(gcf,'color','w','position',[1000,894,1050,425]);

h1 = plot(Xtp,Ytp,'-k','linewidth',0.5); hold on;
h2 = plot(DI,EF,'s','color',colors(1,:),'markersize',8,'markerfacecolor',colors(1,:));
clear h
for s = 1:Nsites
 h(s,:) = plot(DI(s),squeeze(EFm(s,1:end-1)),'.','color',colors(2,:));
end
h3 = plot(DI,LR,'o','color',colors(4,:),'markersize',8,'markerfacecolor',colors(4,:));

h0 = plot([0,1],[0,1],'k--');
plot([1,10.5],[1,1],'k--');

xlabel('Dryness Index','fontsize',16);
ylabel('Evaporative Fraction','fontsize',16);
legend([h1,h2(1),h(1,1),h3(1)],'Turc-Pike','PALS Data','Land Models','Linear Regression','location','nw');%modelNames(1:end-1)');
set(gca,'ylim',[0,1.5]);
set(gca,'xlim',[0,10.5]);

pause(1)

% ----------------------------------------------------------------------------

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

% ----------------------------------------------------------------------------

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

end % Nboot

RMSEs = cat(1,[RMSE.TP;RMSE.LR],RMSE.model(1:end,:));
[~,rank] = sort(mean(RMSEs')); 
RMSEs = RMSEs(rank,:);

xnames = [{'Turc-Pike'};{'EF Regression'};modelNames(1:end)];

figure(2); close(2); figure(2);
set(gcf,'color','w');%,'position',[559,755,642,562])

boxplot(RMSEs');
set(gca,'xticklabel',xnames(rank));
h = gca;
h.XTickLabelRotation = 60;
grid on;
ylabel('MAE of Evap. Frac.','fontsize',16);
set(gca,'ylim',[0,0.5]);

% save models figure
figure(1);
set(gcf,'PaperPositionMode','auto')
saveas(gcf,'figures/BudykoCurve.png');


% save models figure
figure(2);
set(gcf,'PaperPositionMode','auto')
saveas(gcf,'figures/BudykoMAE.png');





