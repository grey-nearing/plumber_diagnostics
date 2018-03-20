clear all
close all
clc

%load('results/all_network_results.mat');
load('results/main_dpn_results_20.mat');
%load('tempsave.mat');
T = T./H;

Tp = squeeze(        T(:,6:end,:,1    ,:));
Hp = squeeze(        H(:,6:end,:,1    ,:));
Ta = squeeze(        T(:,6:end,:,2    ,:));
Ha = squeeze(        H(:,6:end,:,2    ,:));
Tm = squeeze(permute(T(:,6:end,:,3:end,:),[1,2,3,5,4]));
Hm = squeeze(permute(H(:,6:end,:,3:end,:),[1,2,3,5,4]));

Ns = size(Tm,4);
Nm = size(Tm,5);
Dx = size(Tm,1);
Dy = size(Tm,2);

fig = 1;

% plotting tools
figure(fig);close(fig);figure(fig);
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
Scolors(5 ,:) = cc(1).Color;
Scolors(6 ,:) = cc(2).Color;
Scolors(7 ,:) = cc(2).Color;
Scolors(8 ,:) = cc(4).Color;
Scolors(9 ,:) = cc(4).Color;
Scolors(10,:) = cc(5).Color;
Scolors(11,:) = cc(1).Color;
Scolors(12,:) = cc(1).Color;
Scolors(13,:) = cc(5).Color;
Scolors(14,:) = cc(1).Color;
Scolors(15,:) = cc(6).Color;
Scolors(16,:) = cc(5).Color;
Scolors(17,:) = cc(2).Color;
Scolors(18,:) = cc(7).Color;
Scolors(19,:) = cc(2).Color;
Scolors(20,:) = cc(4).Color;

close(fig)

Mmarkers = ['o','p','o','o','o','p','o','p','o','o','p','v','o','p'];
Smarkers = ['o','o','s','o','d','s','d','o','s','o','^','p','s','+','o','d','^','o','p','d'];

%for x = 1:Dx
% for y = 1:Dy
%
%  fig= fig+1; figure(fig); close(fig);figure(fig);
%  set(gcf,'color','w','position',[126,17,1800,1200]);
%
%  for s = 1:Ns
%
%   subplot(4,5,s)
%
%   pdata = squeeze(Tp(x,y,:,s));
%   mdata = squeeze(Tm(x,y,:,s,:));
%
%   semilogx(lags/48,pdata,'k-o','linewidth',2); hold on;
%   for m = 1:Nm
%    semilogx(lags/48,mdata(:,m),'--','color',Mcolors(m,:),'marker',Mmarkers(m),'linewidth',0.5);
%   end
%
%   grid on;
%   if s > 15; xlabel('Lag Time','fontsize',16); end; 
%   ylab = strcat(varNames(x),{' -> '},varNames(5+y));
%   ylabel(ylab,'fontsize',16); 
%   title(siteNames{s},'fontsize',20);
%   %set(gca,'ylim',[0,0.5]);
%
%  end % site
%
% end % y
%end % x

%for m = [1,2,7]
%
% fig= fig+1; figure(fig); close(fig);figure(fig);
% set(gcf,'color','w','position',[126,17,1800,1200]);
%
% sp = 0;
% for y = 1:Dy
%  for x = 1:Dx
%
%   sp = sp +1; subplot(Dx,Dy,sp);
%
%   for s = 1:Ns
%    pdata = squeeze(Tp(x,y,:,s));
%    mdata = squeeze(Tm(x,y,:,s,m));
%    semilogx(lags/48,mdata-pdata,'--','color',Scolors(s,:),'marker',Smarkers(s),'linewidth',1); hold on;
%   end
%
%   grid on;
%%   xlabel('Lag Time','fontsize',16);
%   tit = strcat(varNames(x),{' -> '},varNames(5+y));
%   ylabel(tit,'fontsize',12); 
%%   ylabel(modelNames{m},'fontsize',20);
%  end % site
%
% end % y
%end % x

fig= fig+1; figure(fig); close(fig);figure(fig);
set(gcf,'color','w','position',[126,17,2000,1100]);

xmax = zeros(Dx,Dy);
xmin = ones(Dx,Dy);
sp = 0;
for y = 1:Dy

 for x = 1:Dx

  sp = sp+1; subplot(Dy,Dx,sp);  

  for s = 1:Ns
   l = 1;
   pdata = squeeze(Tp(x,y,l,s));
   mdata = squeeze(Tm(x,y,l,s,:));
   h(s,:) = plot(pdata,mdata,'color',Scolors(s,:),'marker',Smarkers(s)); hold on;
   xmax(x,y) = max([xmax(x,y),pdata,max(mdata)]);
   xmin(x,y) = min([xmin(x,y),pdata,min(mdata)]);
  end
  plot([xmin(x,y)*0.9,xmax(x,y)*1.1],[xmin(x,y)*0.9,xmax(x,y)*1.1],'k--');

  % calculate correlation
  xx = repmat(squeeze(Tp(x,y,l,:)),[1,Nm]);
  yy = squeeze(Tm(x,y,l,:,:)); assert(all(size(xx)==size(yy)));
  I = find(all(~isnan([xx(:),yy(:)]')));
  cc = corrcoef(xx(I),yy(I)); CCpals(x,y) = cc(2);

  grid on;
  path = strcat(varNames(x),{' -> '},varNames(5+y));
  title(path,'fontsize',16);
  if y == 4
   xlabel([strcat({'\rho = '},num2str(round(CCpals(x,y)*100)/100));{'measured'}],'fontsize',12);
  else
   xlabel(strcat('\rho = ',num2str(round(CCpals(x,y)*100)/100)),'fontsize',12);
  end
  if x == 1; ylabel('modeled','fontsize',16); end;
  axis([xmin(x,y)*0.9,xmax(x,y)*1.1,xmin(x,y)*0.9,xmax(x,y)*1.1]);
  axis square
  if y == Dy && x == Dx; legend(h(:,1),siteNames); end;

 end % y
end % x

%fig= fig+1; figure(fig); close(fig);figure(fig);
%set(gcf,'color','w','position',[126,17,2000,2000]);

%xmax = zeros(Dx,Dy);
%sp = 0;
%for y = 1:Dy
% for x = 1:Dx
%
%  sp = sp+1; subplot(Dy,Dx,sp);  
%
%  for s = 1:Ns
%   l = 1;
%   pdata = squeeze(Ta(x,y,l,s));
%   mdata = squeeze(Tm(x,y,l,s,:));
%   h(s,:) = plot(pdata,mdata,'color',Scolors(s,:),'marker',Smarkers(s)); hold on;
%   xmax(x,y) = max([xmax(x,y),max(pdata),max(mdata)]);
%  end
%  plot([0,xmax(x,y)*1.1],[0,xmax(x,y)*1.1],'k--');
%
%  % calculate correlation
%  xx = repmat(squeeze(Ta(x,y,l,:)),[1,Nm]);
%  yy = squeeze(Tm(x,y,l,:,:)); assert(all(size(xx)==size(yy)));
%  I = find(all(~isnan([xx(:),yy(:)]')));
%  cc = corrcoef(xx(I),yy(I)); CCann(x,y) = cc(2);
%
%  grid on;
%  path = strcat(varNames(x),{' -> '},varNames(5+y));
%  title(path,'fontsize',20);
%  axis([0,xmax(x,y)*1.1,0,xmax(x,y)*1.1]);
%  axis square
%  if y == 2 && x == 5; legend(h(:,1),siteNames); end;
%
% end % y
%end % x
%

% save models figure
figure(2);
set(gcf,'PaperPositionMode','auto')
saveas(gcf,'figures/ScattersByPath.png');

