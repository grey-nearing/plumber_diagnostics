clear all
close all
clc
%com.mathworks.services.Prefs.setBooleanPref('EditorGraphicalDebugging',false)
iFig = 0;

% get colors
figure(100);
h = plot(randn(10));
colors = get(h,'Color');
close(100);

% stuff we need for the file name
Nepochs = 10000;
trnfctn = 'trainscg';

% sample sizes
fracs = 100*2.^(0:10);
Nfracs = length(fracs); 

% experiment types
types = [{'all'}];%,{'wet'},{'dry'}];
Ntypes = length(types);

% experiment targets
targs = [{'Qe'},{'Qh'},{'NEE'}];
Ntargs = length(targs);

% number of boodstraps
Nboots = 10;

% loop through experiments
for iType = 1:Ntypes
 for iTarg = 1:Ntargs
  for iBoot = 1:Nboots
   for iFrac = 1:Nfracs

    % file name 
    fname = strcat('../results/',types{iType},'_',targs{iTarg},'_',num2str(fracs(iFrac)),'_',num2str(iBoot),'_',num2str(Nepochs),'_',trnfctn,'.mat');

    % load file
    try
     load(fname);
    catch
     INFO(iType,iTarg,iBoot,iFrac,:,:) = 0/0;
     RMSE(iType,iTarg,iBoot,iFrac,:) = 0/0;
     CORR(iType,iTarg,iBoot,iFrac,:) = 0/0;
     continue
    end 

    for iSamp = 1:3
     for iRes = 1:4
      INFO(iType,iTarg,iBoot,iFrac,iSamp,iRes) = results.stats.info(iSamp,iRes);
      RMSE(iType,iTarg,iBoot,iFrac,iSamp)      = results.stats.rmse(iSamp);
      CORR(iType,iTarg,iBoot,iFrac,iSamp)      = results.stats.corr(iSamp);
     end % iFunc
    end % iRes

    % normalize rmse
    RMSE(iType,iTarg,iBoot,iFrac,1) = 1 - RMSE(iType,iTarg,iBoot,iFrac,1) ./ std(results.Y);
    RMSE(iType,iTarg,iBoot,iFrac,2) = 1 - RMSE(iType,iTarg,iBoot,iFrac,2) ./ std(results.Y(results.Itrn));
    RMSE(iType,iTarg,iBoot,iFrac,3) = 1 - RMSE(iType,iTarg,iBoot,iFrac,3) ./ std(results.Y(results.Itst));

   end % iFrac
  end % iBoot
 end % iTarg
end % iType

% *** Plot Stuff ***********************

iFig = iFig + 1;
figure(iFig); close(iFig); figure(iFig);
%set(gcf,'position',[10,600,2400,600]);
set(gcf,'Renderer','Zbuffer')
set(gcf,'color','w');

for iType = 1:Ntypes
 for iTarg = 1:Ntargs

  % subplot
  subplot(Ntypes,Ntargs,(iType-1)*Ntargs+iTarg);

  % create error bars
  Imu  = squeeze(nanmean(INFO(iType,iTarg,:,:,2:3,:)   ,3));
  Isig = squeeze( nanstd(INFO(iType,iTarg,:,:,2:3,:),[],3));
  Rmu  = squeeze(nanmean(RMSE(iType,iTarg,:,:,2:3)     ,3));
  Rsig = squeeze( nanstd(RMSE(iType,iTarg,:,:,2:3)  ,[],3));
  Cmu  = squeeze(nanmean(CORR(iType,iTarg,:,:,2:3)     ,3));
  Csig = squeeze( nanstd(CORR(iType,iTarg,:,:,2:3)  ,[],3));

  % plot
  for iRes = 1:4
    hi(iRes,1) = errorbar(fracs',Imu(:,1,iRes),Isig(:,1,iRes),'--s','color',colors{iRes}); hold on;
    hi(iRes,2) = errorbar(fracs',Imu(:,2,iRes),Isig(:,2,iRes),':s' ,'color',colors{iRes}); hold on;
  end
  hr(1) = errorbar(fracs',Rmu(:,1),Rsig(:,1),'--s','color',colors{iRes+1}); hold on;
  hr(2) = errorbar(fracs',Rmu(:,2),Rsig(:,2),':s' ,'color',colors{iRes+1}); hold on;
  hc(1) = errorbar(fracs',Cmu(:,1),Csig(:,1),'--s','color',colors{iRes+2}); hold on;
  hc(2) = errorbar(fracs',Cmu(:,2),Csig(:,2),':s' ,'color',colors{iRes+2}); hold on;

  % axis limits
  grid on
  set(gca,'xlim',[0,fracs(end)]);
  set(gca,'ylim',[0,1]);

  % legend
  if iType == 1 && iTarg == 1
   legend('info: 1% resolution - training sample',...
          'info: 1% resolution - test sample',...
          'info: 2% resolution - training sample',...
          'info: 2% resolution - test sample',...
          'info: 5% resolution - training sample',...
          'info: 5% resolution - test sample',...
          'info: 10% resolution - training sample',...
          'info: 10% resolution - test sample',...
          'rmse - training sample',...
          'rmse - test sample',...
          'corr - training sample',...
          'corr - test sample',...
          'location','ne');
  end

  % labels
  xlabel('# Training Points','fontsize',16);
  ylabel('Performance Metric','fontsize',16);
%  title(strcat(targs{iTarg},' - ',types{iType}),'fontsize',18);
  title(strcat(targs{iTarg}),'fontsize',18);

 end
end




