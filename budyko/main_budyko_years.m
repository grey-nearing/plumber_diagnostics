%clear all
%close all
%clc

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

%load('all_data.mat');
%pals = all_data.pals;
%model = all_data.model;


% ----------------------------------------------------------------------------

% init storage
PE = zeros(11,Nsites)./0;
PP = zeros(11,Nsites)./0;
DI = zeros(11,Nsites)./0;
TP = zeros(11,Nsites)./0;
EE = zeros(11,Nsites)./0;
EF = zeros(11,Nsites)./0;
EEm = zeros(11,Nsites,Nmodels)./0;
EFm = zeros(11,Nsites,Nmodels)./0;

% aggregate to annual means
for y = 1:11
 for s = 1:Nsites

  % pull indices for this year
  Iy = find(pals(:,1,s) == y+1995);
  if isempty(Iy); continue; end;
  if length(Iy) < 1e4; continue; end;

  % check that model indexes match
  for m = 1:Nmodels
   Iym = find(model(:,1,m,s) == y+1995);
   assert(max(abs(Iym-Iy))==0);
  end

  % dryness index
  PE(y,s) = mean(model(Iy,9,end,s))/2454000 * 60*30 *length(find(~isnan(model(Iy,9,end,s))));
  PP(y,s) = sum(model(Iy,8,1,s));% * 60*30;% *length(Iy);
  DI(y,s) = PE(y,s)/PP(y,s);
  if DI(y,s) > 10
   PE(y,s) = 0/0;
   PP(y,s) = 0/0;
   DI(y,s) = 0/0;
  end

  % evaporative fraction
  for m = 1:Nmodels
   EEm(y,s,m) = nanmean(model(Iy,9,m,s))/2454000 * 60*30 *length(find(~isnan(model(Iy,9,m,s))));
   EFm(y,s,m) = EEm(y,s,m) / PP(y,s);
  end % models 
  EE(y,s) = mean(pals(Iy,9,s))/2454000 * 60*30 *length(Iy);
  EF(y,s) = EE(y,s)/PP(y,s);

  % Turc-Pike
  if DI(y,s) > 0;
   v = 2;
   TP(y,s) = (1+DI(y,s)^-v)^(-1/v);% * PP(y,s);
  else
   DI(y,s) = 0/0;
   PP(y,s) = 0/0;
   PE(y,s) = 0/0;
   EE(y,s) = 0/0;
   EF(y,s) = 0/0;
   TP(y,s) = 0/0;
   EEm(y,s,:) = 0/0;
   EFm(y,s,:) = 0/0;
  end

  % regression inputs
  Ta(y,s) = mean(pals(Iy,5,s));
  Rn(y,s) = mean(pals(Iy,6,s));

  [s,y,length(Iy)]

 end % sites
end % years

% ----------------------------------------------------------------------------

clear reg LR
for s = 1:Nsites

 % grab loo training/testing data
 S = 1:Nsites; S(s) = [];
 di = DI(:,S); ta = Ta(:,S); rn = Rn(:,S); ef = EF(:,S);
 Xtrain = [di(:),ta(:),rn(:),ones(size(di(:)))];  
 Ytrain = ef(:);               
 Xtest = [DI(:,s),Ta(:,s),Rn(:,s),ones(size(Rn(:,s)))]; 
 Ytest = EF(:,s); Ytest = Ytest(:);

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
 LR(:,s) = Xtest * reg(:,s);

 % screen report
 fprintf('trained/predicted loo model at site %d of %d \n',s,Nsites);

end

% ----------------------------------------------------------------------------

figure(1);
h = plot(randn(7));
for i = 1:7; colors(i,:) = h(i).Color; end;
close all;

Xtp = 0:0.1:3.5;
Ytp = (1+Xtp.^-v).^(-1/v);

figure(1); close(1); figure(1);
set(gcf,'color','w','position',[1000,894,1050,425]);

h1 = plot(Xtp,Ytp,'-k','linewidth',0.5); hold on;
h2 = plot(DI,EF,'s','color',colors(1,:),'markersize',8,'markerfacecolor',colors(1,:));
clear h
for s = 1:Nsites
 h(s,:) = plot(DI(:,s),squeeze(EFm(:,s,1:end-1)),'.','color',colors(2,:));
end
h3 = plot(DI,LR,'o','color',colors(4,:),'markersize',8,'markerfacecolor',colors(4,:));

h0 = plot([0,1],[0,1],'k--');
plot([1,3.5],[1,1],'k--');

xlabel('Dryness Index','fontsize',16);
ylabel('Evaporative Fraction','fontsize',16);
legend([h1,h2(1),h(1,1),h3(1)],'Turc-Pike','PALS Data','Land Models','Linear Regression','location','nw');%modelNames(1:end-1)');
set(gca,'ylim',[0,2]);

pause(5)

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
 RMSE.TP(b) = rmse(TP(I(S)),EF(I(S)));
 RMSE.LR(b) = rmse(LR(I(S)),EF(I(S)));
 for m = 1:Nmodels
  efm = squeeze(EFm(:,:,m));
  RMSE.model(m,b) = rmse(TP(I(S)),efm(I(S)));
 end

end % Nboot

RMSEs = cat(1,[RMSE.TP;RMSE.LR],RMSE.model(1:end-1,:));
[~,rank] = sort(mean(RMSEs')); 
RMSEs = RMSEs(rank,:);

xnames = [{'Turc-Pike'};{'EF Regression'};modelNames(1:end-1)];

figure(2); close(2); figure(2);
set(gcf,'color','w');%,'position',[559,755,642,562])

boxplot(RMSEs');
set(gca,'xticklabel',xnames(rank));
h = gca;
h.XTickLabelRotation = 60;
grid on;
ylabel('RMSE Evaporative Fraction','fontsize',16);








