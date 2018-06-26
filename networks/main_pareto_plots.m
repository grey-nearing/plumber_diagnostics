clear all; close all; clc;
iFig = 0;
restoredefaultpath;

%% *** Experiment Setup ***************************************************

modelNames = [ ...
    {'CABLE 2.0'}
    {'CABLE 2.0 (alt)'}
    {'CHTESSEL'}
    {'COLASSiB 2.0'}
    {'ISBA SURFEX 7.3'}
    {'ISBA SURFEX 7.3 (alt)'}
    {'JULES 3.1'}
    {'JULES 3.1 (alt)'}
    {'Mosaic'}
    {'Noah 2.7.1'}
    {'Noah 3.2'}
    {'Noah 3.3'}
    {'Manabe Bucket'}
    {'Penman Monteith'}];

% site names
siteNames =  [{'Amplero'}
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

% variable names
varNames = [{'Ws'},{'Ta'},{'Rsw'},{'Rlw'},{'Hr'},...
    {'PP'},{'Qe'},{'Qh'} ,{'NEE'},{'SM'}];%,{'SM2'}];

%% *** Load Data **********************************************************

% dpn data
fname = strcat('results/dpn_stats.mat');
load(fname); dpn = outdata; clear outdata;

% mi data
fname = '../benchmarking/results/stats_local_anns.mat';
load(fname);

% only us sites with recorded soil moisture
S = [1,2,3,4,5,7,11,12,13,14,16,17,18];
siteNames = siteNames(S);

% gather model-specific data
TD = repmat(dpn.Tp,[1,1,1,14]) - dpn.Tm;
TD = TD(:,:,S,:);
[Nx,Ny,Ns,Nm] = size(TD);

MI = zeros(Ny,Ns,Nm);
for y = 1:3
    for s = 1:Ns
        for m = 1:Nm
            MI(6+y,s,m) = stats(y).site(S(s)).ann.mi - stats(y).site(S(s)).model(m).mi;
        end
    end
end

%% *** Plotting Tools *****************************************************

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

% Scolors(1 ,:) = cc(1).Color;
% Scolors(2 ,:) = cc(2).Color;
% Scolors(3 ,:) = cc(1).Color;
% Scolors(4 ,:) = cc(2).Color;
% Scolors(5 ,:) = cc(3).Color;
% Scolors(6 ,:) = cc(2).Color;
% Scolors(7 ,:) = cc(1).Color;
% Scolors(8 ,:) = cc(4).Color;
% Scolors(9 ,:) = cc(4).Color;
% Scolors(10,:) = cc(5).Color;
% Scolors(11,:) = cc(2).Color;
% Scolors(12,:) = cc(2).Color;
% Scolors(13,:) = cc(5).Color;
% Scolors(14,:) = cc(2).Color;
% Scolors(15,:) = cc(6).Color;
% Scolors(16,:) = cc(5).Color;
% Scolors(17,:) = cc(2).Color;
% Scolors(18,:) = cc(7).Color;
% Scolors(19,:) = cc(2).Color;
% Scolors(20,:) = cc(4).Color;

close(1);

% plot marker shapes
Mmarkers = ['o','p','o','o','o','p','o','p','o','o','p','v','o','p'];
Smarkers = ['o','o','s','s','o','p','p','^','v','o','+','>','s','<',...
            'o','p','h','d','o','s'];
        
%% *** Make Site-Grouped and Model-Grouped Plots **************************

close all

% init figure
iFig = 1; figure(iFig); close(iFig); figure(iFig);
set(gcf,'color','w','position',[318     1   728   804]);

iFig = 2; figure(iFig); close(iFig); figure(iFig);
set(gcf,'color','w','position',[318     1   728   804]);

% loop through subplots
sp = 0;
for x = [2:3,9:10]
    for y = [7:9]
        
        % by model subplot
        figure(1);
        
        % advance subplot
        sp = sp+1; subplot(4,3,sp);
        
        for m = 1:Nm
    
            % gather model-specific data
            xdata = squeeze(TD(x,y,:,m));
            ydata = squeeze(MI(y,:,m));
            
            % plot  
            hm(m) = plot(xdata,ydata,           ...
                Mmarkers(m),'markersize',10, ...
                'color',Mcolors(m,:),'linewidth',2);
            hold on;
            
        end
        
        % transfer path (as subplot title)
        tit = strcat(varNames(x),{' -> '},varNames(y));
        title(tit{1},'fontsize',22);
        
        % bounds
        yl = ylim;
        set(gca,'ylim',[0,yl(2)]);
        plot([0,0],[0,yl(2)],'--k','linewidth',0.5);
        
        % axis labels
        if x == 10; xlabel('\Delta transfer entropy','fontsize',16); end;
        if y == 7; ylabel('missing info','fontsize',16); end;
        
        % legend
%         if sp == 9; legend(hm,modelNames); end
        
        % ------------
        
        % by site subplot
        figure(2);
        
        % advance subplot
        subplot(4,3,sp);
        
        for s = 1:Ns
            
            % gather model-specific data
            xdata = squeeze(TD(x,y,s,:));
            ydata = squeeze(MI(y,s,:));
            
            % plot
            hs(s) = plot(xdata,ydata,...           ...
                Smarkers(s),'markersize',10, ...
                'color',Scolors(s,:),'linewidth',2);
            hold on;
        end
        
        % transfer path (as subplot title)
        tit = strcat(varNames(x),{' -> '},varNames(y));
        title(tit{1},'fontsize',20);
        
        % bounds
        yl = ylim;
        set(gca,'ylim',[0,yl(2)]);
        plot([0,0],[0,yl(2)],'--k','linewidth',0.5);
       
        % axis labels
        if x == 10; xlabel('\Delta transfer entropy','fontsize',16); end;
        if y == 7; ylabel('missing info','fontsize',16); end;
        
        % legend
%         if sp == 9; legend(hs,siteNames); end
        
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

%% *** Rn-Specific Plots **************************************************

figure(3); close(3); figure(3);
set(gcf,'color','w','position',[49         180        1033         600]);

% ------

subplot(2,3,[1,2])
ylabel('mutual info','fontsize',16); 
iSend = 3; iTarg = 7;
tit = strcat(varNames(iSend),{' -> '},varNames(iTarg));

for iModel = 1:Nm%[1,2,5,6,10,11,12]
 plot(repmat(iModel,[1,Ns]),squeeze(TD(iSend,iTarg,:,iModel)), ...
   Mmarkers(iModel),'markersize',10, ...
   'color',Mcolors(iModel,:),'linewidth',2);%,'markerfacecolor',Mcolors(iModel,:));
 hold on;
end
plot([0.5,Nm+0.5],[0,0],'k--','linewidth',2)
set(gca,'xtick',1:Nm,'xticklabel',modelNames);
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
legend(modelNames)
%legend(modelLegen{[1,2,5,6,10,11,12]},'location','nw')
xlabel('transfer diff','fontsize',16); 
title(tit{1},'fontsize',20);
xl = xlim; yl = ylim; plot([0,0],yl,'--k','linewidth',0.5);

% ------

subplot(2,3,[4,5])
iSend = 3; iTarg = 8;
tit = strcat(varNames(iSend),{' -> '},varNames(iTarg));

for iModel = 1:Nm%[1,2,5,6,10,11,12]
 plot(repmat(iModel,[1,Ns]),squeeze(TD(iSend,iTarg,:,iModel)), ...
   Mmarkers(iModel),'markersize',10, ...
   'color',Mcolors(iModel,:),'linewidth',2);%,'markerfacecolor',Mcolors(iModel,:));
 hold on;
end
%boxplot(squeeze(TD(iSend,iTarg,:,:)));
%hold on
plot([0.5,Nm+0.5],[0,0],'k--','linewidth',2)
set(gca,'xtick',1:Nm,'xticklabel',modelNames);
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
saveas(gcf,'figures/Figre - 6 Rn_Qe_By_Group.png');


%% *** End Program ********************************************************










