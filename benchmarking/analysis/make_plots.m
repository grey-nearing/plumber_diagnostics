clear all
close all
clc
iFig = 0;

% get colors
figure(100);
h = plot(randn(20));
colors = get(h,'Color');
close(100);

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
Mcolors(14,:) = 0.0*[1,1,1];
close(1)

Mmarkers = ['o','p','o','o','o','p','o','p','o','o','p','v','o','o'];

% load results
load('loo_site_results3.mat');

a = INFO(15,:,:,:);
b = INFO(11,:,:,:);
c = INFO(12,:,:,:);
INFO(11,:,:,:) = c;
INFO(15,:,:,:) = b;
INFO(12,:,:,:) = a;

a = modelNames{13};
b = modelNames{9};
c = modelNames{10};
modelNames{9} = c;
modelNames{13} = b;
modelNames{10} = a;

% ---------------------------------------

% plot
iFig = iFig+1;
figure(iFig); close(iFig); figure(iFig);
set(gcf,'color','w');
set(gcf,'position',[250 450 1600 815]);

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
 h(iRes).FaceColor = 1-iRes/Nres * ones(3,1);
end
for iMod = 1:Nmodels
 mdat = squeeze(mu(2+iMod,:,:))'; mdat = mdat(:);
 for iTarg = 1:Ntargs
  sdex = (iTarg-1)*Nres + 1;
  edex = sdex + Nres - 1;
  hm(iMod) = plot(xx(sdex:edex),mdat(sdex:edex),'o-','linewidth',1,'color',colors{iMod});  
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

 if iTarg == 1; leg = legend(h4,resNames,'location','nw'); end;
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
