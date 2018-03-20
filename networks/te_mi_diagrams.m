clear all
close all
clc

% raw data
load('results/tradeoff_results.mat');

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
Scolors(2 ,:) = cc(2).Color;
Scolors(3 ,:) = cc(1).Color;
Scolors(4 ,:) = cc(3).Color;
Scolors(5 ,:) = cc(2).Color;
Scolors(6 ,:) = cc(1).Color;
Scolors(7 ,:) = cc(4).Color;
Scolors(8 ,:) = cc(2).Color;
Scolors(9 ,:) = cc(4).Color;
Scolors(10,:) = cc(2).Color;
Scolors(11,:) = cc(4).Color;
Scolors(12,:) = cc(2).Color;
Scolors(13,:) = cc(5).Color;

%Scolors(1 ,:) = cc(1).Color;
%Scolors(2 ,:) = cc(2).Color;
%Scolors(3 ,:) = cc(1).Color;
%Scolors(4 ,:) = cc(2).Color;
%Scolors(5 ,:) = cc(3).Color;
%Scolors(6 ,:) = cc(2).Color;
%Scolors(7 ,:) = cc(1).Color;
%Scolors(8 ,:) = cc(4).Color;
%Scolors(9 ,:) = cc(4).Color;
%Scolors(10,:) = cc(5).Color;
%Scolors(11,:) = cc(2).Color;
%Scolors(12,:) = cc(2).Color;
%Scolors(13,:) = cc(5).Color;
%Scolors(14,:) = cc(2).Color;
%Scolors(15,:) = cc(6).Color;
%Scolors(16,:) = cc(5).Color;
%Scolors(17,:) = cc(2).Color;
%Scolors(18,:) = cc(7).Color;
%Scolors(19,:) = cc(2).Color;
%Scolors(20,:) = cc(4).Color;

close(1);

Mmarkers = ['o','p','o','o','o','p','o','p','o','o','p','v','o','p'];
Smarkers = ['o','o','s','s','o','p','p','^','v','o','+','>','s','<','o','p','h','d','o','s'];

% plot marker shapes

modelLegen = [ ...
              {'CABLE 2.0'}
              {'CABLE 2.0 (alt)'}
              {'CHTESSEL'}
              {'COLASSiB 2.0'}
              {'ISBA SURFEX 7.3'}
              {'ISBA SURFEX 7.3 (alt)'}
              {'JULES 3.1'}
              {'JULES 3.1 (alt)'}
              {'Mosaic'}
              {'NOAH 2.7.1'}
              {'NOAH 3.2'}
              {'NOAH 3.3'}
              {'Manabe Bucket'}
              {'Penman Monteith'}];

% init figure
iFig = 1; figure(iFig); close(iFig); figure(iFig);
set(gcf,'color','w','position',[46,37,1000,1450]);

iFig = 2; figure(iFig); close(iFig); figure(iFig);
set(gcf,'color','w','position',[46,37,1000,1450]);

% loop through subplots
sp = 0;
for iSend = [2:3,8:9];
 for iTarg = 1:Ntargs

  % by model subplot
  figure(1);

  % advance subplot
  sp = sp+1; subplot(4,Ntargs,sp);

  for iModel = 1:Nm
%  for iModel = [1,2,5,6,10,11,12]
   xdata = squeeze(TD(iSend,iTarg,:,iModel));
   ydata = squeeze(MI(iTarg,:,iModel));
   hm(iModel) = plot(xdata,ydata,           ...
        Mmarkers(iModel),'markersize',10, ...
        'color',Mcolors(iModel,:),'linewidth',2);%,'markerfacecolor',Mcolors(iModel,:));
   hold on; 
%   grid on;
  end

  % transfer path (as subplot title)
  tit = strcat(varNames(iSend),{' -> '},varNames(5+iTarg));
  title(tit{1},'fontsize',20);

  % bounds
  xl = 0.12;
  yl = 0.25;
  axis([-xl,xl,0,yl])

  % aesthetics
  xl = xlim;
  yl = ylim;
  plot([0,0],yl,'--k','linewidth',0.5);

  % axis labels
  if iSend == 9; xlabel('transfer diff','fontsize',16); end;
  if iTarg == 1; ylabel('missing info','fontsize',16); end;

  % legend
  if sp == 9; legend(hm,modelLegen); end

  % by site subplot
  figure(2);

  % advance subplot
  subplot(4,Ntargs,sp);

  for iSite = 1:Ns
   xdata = squeeze(TD(iSend,iTarg,iSite,:));
   ydata = squeeze(MI(iTarg,iSite,:));
   hs(iSite) = plot(xdata,ydata,... %repmat(budyko(iSite,2)/budyko(iSite,1),1,Nmodels),           ...
        Smarkers(iSite),'markersize',10, ...
        'color',Scolors(iSite,:),'linewidth',2);%,'markerfacecolor',Scolors(iSite,:));
   hold on; 
  end

  % transfer path (as subplot title)
  tit = strcat(varNames(iSend),{' -> '},varNames(5+iTarg));
  title(tit{1},'fontsize',20);

  % bounds
  xl = 0.12;
  yl = 0.25;
  axis([-xl,xl,0,yl])

  % aesthetics
  xl = xlim;
  yl = ylim;
  plot([0,0],yl,'--k','linewidth',0.5);

  % axis labels
  if iSend == 9; xlabel('transfer diff','fontsize',16); end;
  if iTarg == 1; ylabel('missing info','fontsize',16); end;

  % legend
  if sp == 9; legend(hs,siteNames); end

 end % iSend
end % iTarg

% save models figure
h = figure(1);
set(h,'PaperPositionMode','auto')
saveas(h,'figures/ProcessByModels.png');

% save sites figure
h = figure(2);
set(h,'PaperPositionMode','auto')
saveas(h,'figures/ProcessBySites.png');

% ---------------------------------------------------------------------------------------

figure(3); close(3); figure(3);
set(gcf,'color','w','position',[49,389,1033,391]);

% ------

subplot(2,3,[1,2])
ylabel('mutual info','fontsize',16); 
iSend = 3; iTarg = 1;
tit = strcat(varNames(iSend),{' -> '},varNames(5+iTarg));

for iModel = 1:Nm%[1,2,5,6,10,11,12]
 plot(repmat(iModel,[1,Ns]),squeeze(TD(iSend,iTarg,:,iModel)), ...
   Mmarkers(iModel),'markersize',10, ...
   'color',Mcolors(iModel,:),'linewidth',2);%,'markerfacecolor',Mcolors(iModel,:));
 hold on;
end
plot([0.5,Nm+0.5],[0,0],'k--','linewidth',2)
set(gca,'xtick',1:Nm,'xticklabel',modelLegen);
ylabel('transfer diff','fontsize',16);
ax = gca;
ax.XTickLabelRotation=45;
title(tit{1},'fontsize',20);

subplot(2,3,3)
for iModel = 1:Nm%[1,2,5,6,10,11,12]
 xdata = squeeze(TD(iSend,iTarg,:,iModel));
 ydata = squeeze(MI(iTarg,:,iModel));
 hm(iModel) = plot(xdata,ydata,           ...
      Mmarkers(iModel),'markersize',10, ...
      'color',Mcolors(iModel,:),'linewidth',2);%,'markerfacecolor',Mcolors(iModel,:));
 hold on;
end
legend(modelLegen)
%legend(modelLegen{[1,2,5,6,10,11,12]},'location','nw')
xlabel('transfer diff','fontsize',16); 
title(tit{1},'fontsize',20);
xl = xlim; yl = ylim; plot([0,0],yl,'--k','linewidth',0.5);

% ------

subplot(2,3,[4,5])
iSend = 3; iTarg = 2;
tit = strcat(varNames(iSend),{' -> '},varNames(5+iTarg));

for iModel = 1:Nm%[1,2,5,6,10,11,12]
 plot(repmat(iModel,[1,Ns]),squeeze(TD(iSend,iTarg,:,iModel)), ...
   Mmarkers(iModel),'markersize',10, ...
   'color',Mcolors(iModel,:),'linewidth',2);%,'markerfacecolor',Mcolors(iModel,:));
 hold on;
end
%boxplot(squeeze(TD(iSend,iTarg,:,:)));
%hold on
plot([0.5,Nm+0.5],[0,0],'k--','linewidth',2)
set(gca,'xtick',1:Nm,'xticklabel',modelLegen);
ylabel('transfer diff','fontsize',16);
ax = gca;
ax.XTickLabelRotation=45;
title(tit{1},'fontsize',20);

subplot(2,3,6)
for iModel = 1:Nm%[1,2,5,6,10,11,12]
 xdata = squeeze(TD(iSend,iTarg,:,iModel));
 ydata = squeeze(MI(iTarg,:,iModel));
 hm(iModel) = plot(xdata,ydata,           ...
      Mmarkers(iModel),'markersize',10, ...
      'color',Mcolors(iModel,:),'linewidth',2);%,'markerfacecolor',Mcolors(iModel,:));
 hold on;
end
xlabel('transfer diff','fontsize',16); 
ylabel('mutual info','fontsize',16); 
title(tit{1},'fontsize',20);
xl = xlim; yl = ylim; plot([0,0],yl,'--k','linewidth',0.5);

% save sites figure
figure(3);
set(gcf,'PaperPositionMode','auto')
saveas(gcf,'figures/Rn_Qe_By_Group.png');












