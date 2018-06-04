clear all; close all; clc

%% *** Initialize Plotting Tools ******************************************

% figure number
iFig = 0;

% get colors
figure(100);
h = plot(randn(20));
colors = get(h,'Color');
close(100);

% model colors
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
Mcolors(14,:) = 0.0*[1,1,1];
close(1)

% model markers
Mmarkers = ['v','v','p','s','s','^','^','*','h','o','o','o','+','+'];

% target names
targNames = [{'Qe'},{'Qh'},{'NEE'}];

% stats names
statNames = [   {'Root Mean Squared Error'},
                {'Mean Bias Error'},
                {'Correlation Coefficient (R)'},
                {'Normalizd Mean Absolute Error'},
                {'5th Percetile Discrepancy'},
                {'95th Percentil Discrepancy'},
                {'2nd Moment discrepancy'},
                {'3rd Moment Discrepancy'},
                {'4th Moment Discrepancy'},
                {'Distribution Overlap'},
                {'Normalized Mutual Information'}];
            
% model names
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
    %              {'ORCHIDEE.trunk_r1401'}
    %              {'SUMMA.1.0.exp.01.000'}
    %              {'SUMMA.1.0.exp.01.001'}
    %              {'SUMMA.1.0.exp.01.002'}
    %              {'SUMMA.1.0.exp.01.003'}
    %              {'SUMMA.1.0.exp.01.004'}
    %              {'SUMMA.1.0.exp.01.005'}
    %              {'SUMMA.1.0.exp.01.006'}
    %              {'SUMMA.1.0.exp.01.007'}
    %              {'SUMMA.1.0.exp.01.008'}
    %              {'SUMMA.1.0.exp.01.009'}
    {'Manabe_Bucket.2 '}
    {'Penman_Monteith.1'}];

%% *** Load Data **********************************************************

% load local results
fname = './results/stats_local_anns.mat';
load(fname);
localStats = stats; clear stats;

% load global results
fname = './results/stats_global_anns.mat';
load(fname);
globalStats = stats; clear stats;

% count number of stats
statsMark = fieldnames(localStats(1).all.ann);
Nstats = length(statsMark);

% other dimensions
Ny = length(globalStats);
Nmodels = length(globalStats(1).all.model);
Nsites = length(localStats(1).site);

%% *** Plot Local Mutual Information Stats (Local, Global, Models) ********

% increment figure number
iFig = iFig+1;

% set up figure
figure(iFig); close(iFig); figure(iFig);
set(gcf,'color','w');
set(gcf,'position',[2121         821          824         433]);

% which stat to use
s = Nstats;

% prepare ann data
annData = zeros(Ny,2);
for y = 1:Ny
    annData(y,1) = localStats(y).all.ann.(statsMark{s});
    annData(y,2) = globalStats(y).all.ann.(statsMark{s});
end % y-loop

% prepare model data
modelData = zeros(Ny,Nmodels);
for y = 1:Ny
    for m = 1:Nmodels
        modelData(y,m) = globalStats(y).all.model(m).(statsMark{s});
    end
end % y-loop


% plot ann data
b = bar(annData); hold on;
b(1).FaceColor = 0.8*ones(3,1);
b(2).FaceColor = 0.3*ones(3,1);

% plot model data
for m = 1:Nmodels
    plot((1:3)-0.125,modelData(:,m),'--',...
        'linewidth',1,...
        'markersize',10,...
        'marker',Mmarkers(m));
end

% aesthetics
set(gca,'fontsize',22);
grid on

% labels
ylabel(statNames{s},'fontsize',24);
set(gca,'xticklabel',targNames);
leg = legend([{'Local ANN'};{'Global ANN'};modelNames],'location','se');
leg.FontSize = 9;
    
%% *** Make Plumber-Type Plots ********************************************

% stats gorupings used by Best et al. (2015)
I1 = [2,3,4,7];
I2 = [5,6];
I3 = [8,9,10];

% stats ratios plotting nonsense
updn = [-1,-1,1,-1,-1,-1,0,-1,-1,1,1];

% calculate rankings
stat = zeros(4,Nsites,Nmodels-2,Ny,Nstats)/0;
rank = zeros(4,Nsites,Nmodels-2,Ny,Nstats)/0;
group = zeros(4,Ny,Nmodels,3);
for m = 1:Nmodels-2
    for y = 1:Ny
        for s = 1:Nstats
            for i = 1:Nsites
                
                % concatenate stats
                stat(:,i,m,y,s)  = [globalStats(y).site(i).ann.(statsMark{s}),...
                    globalStats(y).site(i).model(m).(statsMark{s}),...
                    globalStats(y).site(i).model(Nmodels-1).(statsMark{s}),...
                    globalStats(y).site(i).model(Nmodels).(statsMark{s})];
                
                if updn(s) ~= 0
                    stat(:,i,m,y,s) = -updn(s)*stat(:,i,m,y,s);
                else
                    stat(:,i,m,y,s) = abs(stat(:,i,m,y,s));
                end
                
                % calculate ranks
                [~,rank(:,i,m,y,s)] = sort(stat(:,i,m,y,s));
            end % i-loop
        end % s-loop
        
        % average rankings
        for r = 1:4
            a = squeeze(rank(r,:,m,y,I1)); group(r,y,m,1) = mean(a(:));
            a = squeeze(rank(r,:,m,y,I2)); group(r,y,m,2) = mean(a(:));
            a = squeeze(rank(r,:,m,y,I3)); group(r,y,m,3) = mean(a(:));
        end % r-loop
        
    end % y-loop
end % m-loop

% remove bad data
group(3:4,3,:,:) = 0/0;
group(1:2,[3:6,8:12],:,:) = 0/0;

% plot
for s = 1:Nstats
    
    % set up figure
    iFig = iFig+1;
    figure(iFig); close(iFig); figure(iFig);
    set(gcf,'color','w');
    set(gcf,'position',[2024         131        1154        1246]);
    
    for y = 1:Ny-1
        
        % subplot
        subplot(Ny-1,1,y);
        
        % plot
%         bar(squeeze(group(:,y,1:Nmodels-2,g))')
        bar(squeeze(mean(rank(:,:,1:Nmodels-2,y,s),2))');
        
        % labels
        set(gca,'fontsize',16);
        ylabel('Average Rank','fontsize',20);
        title(strcat(statNames{s},' -- ',targNames{y}),'fontsize',20);
        set(gca,'xticklabel',modelNames);
        xtickangle(gca,60)
        
    end %y-loop
end % g-loop


asdf

mu  = squeeze(nanmean(INFO,3));
sig = squeeze(nanstd(INFO,[],3));

[numgroups,numbars] = size(squeeze(mu(1,:,:)));
groupwidth = min(0.8,numbars/(numbars+1.5));
xx = [];
for i = 1:numbars
 xx = [xx;(1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars)]; 
end
xx = xx(:);

% plot
h = bar(squeeze(mu(2,:,:))); hold on;
for iRes = 1:Nres
% h(iRes).FaceColor = 1-iRes/Nres * ones(3,1);
end
for iMod = 1:Nmodels
 mdat = squeeze(mu(2+iMod,:,:))'; mdat = mdat(:);
 for iTarg = 1:Ntargs
  sdex = (iTarg-1)*Nres + 1;
  edex = sdex + Nres - 1;
  hm(iMod) = plot(xx(sdex:edex),mdat(sdex:edex),'o-',...
      'linewidth',1,...
      'color',colors{iMod});  
 end
end

% labels
set(gca,'xticklabel',targNames,'fontsize',16);%,'fontsize',16);
ylabel('Info Ratio','fontsize',16);
grid on;
title('Out-of-Sample Benchmark','fontsize',18)

% legend
leg = legend([h,hm],[resNames';modelNames']);
set(leg,'fontsize',12)
set(leg,'position',[0.84 0.47 0.15 0.50]);

% save models figure
figure(iFig);
set(gcf,'PaperPositionMode','auto')
saveas(gcf,'figures/NecessaryBenchmark.png');

%% ---------------------------------------


iFig = iFig+1;
figure(iFig); close(iFig); figure(iFig);
set(gcf,'color','w','position',[1400,100,750,1150]);

for iTarg = 1:Ntargs
 subplot(Ntargs,1,iTarg)

 mu  = squeeze(nanmean(INFO(1:2,iTarg,:,:),3));
 sig = squeeze(nanstd(INFO(1:2,iTarg,:,:),[],3));

 loss = squeeze(nanmean((INFO(2,iTarg,:,:)-INFO(1,iTarg,:,:))./INFO(2,iTarg,:,:),3))';
% mu = cat(1,mu,loss);

 [numgroups,numbars] = size(squeeze(mu));
 groupwidth = min(0.8,numbars/(numbars+1.5));
 xx = [];
 for i = 1:numbars
  xx = [xx;(1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars)]; 
 end
 xx = xx';

 % plot
 [ax,h1,h2] = plotyy(xx(:),mu(:),1+xx(2,:),loss,'bar','bar'); hold on;
 h1.FaceColor = 'none';
 h2.FaceColor = 'none';
 h1.EdgeColor = 'none';
 h2.EdgeColor = 'none';
 h3 = bar(cat(1,mu,loss));
 h4 = bar(mu);
 for iRes = 1:Nres
  frac = (iRes-1)/Nres;
  h3(iRes).FaceColor = (1-frac)*[0.850,0.325,0.098];
  frac = (iRes)/Nres;
  h4(iRes).FaceColor = (1-frac) * [1,1,1];
 end
 ax(2).YColor = colors{2};

 if iTarg == 1; leg = legend(h4,resNames,'location','nw'); end
 title(targNames{iTarg},'fontsize',18);
 set(gca,'xtick',[1,2,3],'xticklabel',[{'LOO'},{'in-Sample'},{'fractional loss'}]);
 ax(1).YLabel.String = 'Info Ratio';
 ax(1).YLabel.FontSize = 18;
 ax(2).YLabel.String = 'Info Loss Fraction';
 ax(2).YLabel.FontSize = 18;
 set(gca,'xlim',[0.5,3.5]);
 set(gca,'ylim',[0,0.8],'ytick',0:0.1:0.8);
 grid on

end

% save models figure
figure(iFig);
set(gcf,'PaperPositionMode','auto')
saveas(gcf,'figures/NecessaryVsSufficient.png');

% ---------------------------------------

iFig = iFig+1;
figure(iFig); close(iFig); figure(iFig);
set(gcf,'color','w');
markers = ['o','s','v'];
ax = plot([0,1],[1,1],'w','linewidth',10); hold on;
axis square
set(gcf,'position',[1000,425,1000,850]);
set(gca,'fontsize',12);

for iTarg = Ntargs:-1:1
% subplot(Ntargs,1,iTarg)

 Hmu = squeeze(nanmean(H(          iTarg,:,4),2));
 Umu = squeeze(nanmean(INFO(2     ,iTarg,:,4),3));
 Mmu = squeeze(nanmean(INFO(3:end ,iTarg,:,4),3)); 
  
 for m = 1:Nmodels
  data(1,m) = (Hmu-Umu)   /(Hmu);%-Mmu(m));
  data(2,m) = (Umu-Mmu(m))/(Hmu);%-Mmu(m));
  data(3,m) = (Mmu(m))    /(Hmu);%-Mmu(m));
 end

 plot([0,1/2],[0,sin(pi/3)],'k-','linewidth',3); hold on;
 plot([1,1/2],[0,sin(pi/3)],'k-','linewidth',3);
 plot([0,1  ],[0,0],'k-','linewidth',3);
 for g = 0:0.1:1
  plot([  g,  g/2],[0,g/2*tan(pi/3)],'k--');
  plot([1-g,1-g/2],[0,g/2*tan(pi/3)],'k--');
  plot([g/2,1-g/2],g/2*tan(pi/3)*[1,1],'k--');
  h = text((2-g)/2+0.01,g/2*tan(pi/3)+0.01,num2str(g),'fontsize',12); set(h,'rotation',0);
  if g == 0 || g == 1
   h = text(g,-0.02,num2str(g),'fontsize',12); set(h,'rotation',60); 
   h = text(g/2-0.02,g/2*tan(pi/3)+0.02,num2str(1-g),'fontsize',12); set(h,'rotation',-60)
  else
   h = text(g-0.02,-0.045,num2str(g),'fontsize',12); set(h,'rotation',60); 
   h = text(g/2-0.03,g/2*tan(pi/3)+0.04,num2str(1-g),'fontsize',12); set(h,'rotation',-60)
  end
 end

 for m = Nmodels:-1:1
  h(m) = plot(data(2,m)+data(3,m)/2,data(3,m)*sin(pi/3),'linestyle','none',...
    'color',Mcolors(m,:),...
    'marker',Mmarkers(m),...
    'markersize',15,'linewidth',3);%,...
%    'markerfacecolor',Mcolors(m,:));
 end
end % iTarg
 
leg = legend(h,modelNames,'location','ne');
set(leg,'position',[0.73 0.57 0.19 0.30]);

set(gca,'ycolor','w');
set(gca,'xcolor','w');
set(gca,'ytick',[]);
set(gca,'xtick',0:0.1:1);

h = text(0.3,-0.07,'Info Lost due to Model Error','fontsize',16);
h = text(0.080,0.30,'Info Missing from Forcings','fontsize',16);
set(h,'rotation',60);
h = text(0.73,0.65,'Info Provided by Model','fontsize',16);
set(h,'rotation',-60);
title('Components of Total Entropy in Observations','fontsize',20);

text(0.22,0.25,'Q_h','fontsize',20);
text(0.50,0.43,'Q_e','fontsize',20);
text(0.36,0.32,'NEE','fontsize',20);

% save models figure
figure(iFig);
set(gcf,'PaperPositionMode','auto')
saveas(gcf,'figures/FracMissingInfo.png');
