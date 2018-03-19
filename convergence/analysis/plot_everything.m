clear all
close all
clc

Nfracs = 10; % number of sample sizes
Fmin = 1e-3; % miniumum sample size
Fmax = 8e-1; % maximum sample size
fracs = 1e5*logspace(log10(Fmin),log10(Fmax),Nfracs);

models = [ ...
          {'CABLE.2.0'}
          {'CABLE_2.0_SLI.vxh599_r553'}
          {'CHTESSEL'}
          {'COLASSiB.2.0'}
          {'ISBA_SURFEX_3l.SURFEX7.3'}
          {'ISBA_SURFEX_dif.SURFEX7.3'}
          {'JULES.3.1'}
          {'JULES3.1_altP'}
          {'Manabe_Bucket.2 '}
          {'Mosaic.1'}
          {'NOAH.2.7.1'}
          {'Noah.3.2'}
          {'NOAH.3.3'}
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
          {'Penman_Monteith.1'}];
Nmodels = length(models);

% load model data
model = load('../../data/model_data/model_info.mat');
mINFO = model.INFO;

% load ann data
ann = load('ann_info.mat');
aINFO = ann.INFO;

% get colors
figure(100);
h = plot(randn(10));
colors = get(h,'Color');
close(100);

% plot Qe
figure(1); close(1); figure(1);
set(gcf,'color','w');

% create ann error bars
amu  = squeeze(nanmean(aINFO.E,2));
asig = squeeze( nanstd(aINFO.E,[],2));

% create model error bars
mmu  = squeeze(nanmean(mINFO.E,2));
msig = squeeze( nanstd(mINFO.E,[],2));

% separate physical benchmarks
bmu  = mmu(:,[9,14],:);   mmu(:,[9,14],:) = [];
bsig = msig(:,[9,14],:);  msig(:,[9,14],:) = [];

% plot
for iRes = 1:4

 % initiate plot
 clear h
 subplot(2,2,iRes)
 
 % plot anns
 a(iRes,2) = errorbar(fracs',amu(:,2,iRes),asig(:,2,iRes),'--s','color',colors{1},'linewidth',2); hold on;
 a(iRes,1) = errorbar(fracs',amu(:,1,iRes),asig(:,1,iRes),'--o','color',colors{2},'linewidth',2); hold on;
% a(iRes,3) = errorbar(fracs',amu(:,3,iRes),asig(:,3,iRes),':'  ,'color',colors{3},'linewidth',2); hold on;

 % plot models and physical benchmarks
 m(iRes,:) = errorbar(repmat(fracs',1,size(mmu,2)),mmu(:,:,iRes),msig(:,:,iRes),'-','color',colors{4}   ,'linewidth',1);
 b(iRes,:) = errorbar(repmat(fracs',1,size(bmu,2)),bmu(:,:,iRes),bsig(:,:,iRes),'-','color',colors{5}   ,'linewidth',1);

 % axis limits
 grid on
 set(gca,'xlim',[0,0.8]*1e5);
 set(gca,'ylim',[0,0.75]);

 % title
 if iRes == 1; title('Qe - 1% Info Res','fontsize',16); end
 if iRes == 2; title('Qe - 2% Info Res','fontsize',16); end
 if iRes == 3; title('Qe - 5% Info Res','fontsize',16); end
 if iRes == 4; title('Qe - 10% Info Res','fontsize',16); end

 % labels
 xlabel('# Training Points','fontsize',16); 
 ylabel('Info Frac [~]','fontsize',16); 

 % legend
 if iRes == 2; legend([a(2,2),a(2,1),m(2,1),b(2,1)],'ANN Training Data','ANN Test Data','Institutional Models','Physical Benchmarks','location','ne'); end;

end


% plot Qh
figure(2); close(2); figure(2);
set(gcf,'color','w');

% create ann error bars
amu  = squeeze(nanmean(aINFO.H,2));
asig = squeeze( nanstd(aINFO.H,[],2));

% create model error bars
mmu  = squeeze(nanmean(mINFO.H,2));
msig = squeeze( nanstd(mINFO.H,[],2));

% separate physical benchmarks
bmu  = mmu(:,[9,14],:);   mmu(:,[9,14],:) = [];
bsig = msig(:,[9,14],:);  msig(:,[9,14],:) = [];

% plot
for iRes = 1:4

 % initiate plot
 clear h
 subplot(2,2,iRes)
 
 % plot anns
 a(iRes,2) = errorbar(fracs',amu(:,2,iRes),asig(:,2,iRes),'--s','color',colors{1},'linewidth',2); hold on;
 a(iRes,1) = errorbar(fracs',amu(:,1,iRes),asig(:,1,iRes),'--o','color',colors{2},'linewidth',2); hold on;
% a(iRes,3) = errorbar(fracs',amu(:,3,iRes),asig(:,3,iRes),':'  ,'color',colors{3},'linewidth',2); hold on;

 % plot models and physical benchmarks
 m(iRes,:) = errorbar(repmat(fracs',1,size(mmu,2)),mmu(:,:,iRes),msig(:,:,iRes),'-','color',colors{4}   ,'linewidth',1);
 b(iRes,:) = errorbar(repmat(fracs',1,size(bmu,2)),bmu(:,:,iRes),bsig(:,:,iRes),'-','color',colors{5}   ,'linewidth',1);

 % axis limits
 grid on
 set(gca,'xlim',[0,0.8]*1e5);
 set(gca,'ylim',[0,0.75]);

 % title
 if iRes == 1; title('Qh - 1% Info Res','fontsize',16); end
 if iRes == 2; title('Qh - 2% Info Res','fontsize',16); end
 if iRes == 3; title('Qh - 5% Info Res','fontsize',16); end
 if iRes == 4; title('Qh - 10% Info Res','fontsize',16); end

 % labels
 xlabel('# Training Points','fontsize',16); 
 ylabel('Info Frac [~]','fontsize',16); 

 % legend
 if iRes == 1; legend([a(2,2),a(2,1),m(2,1),b(2,1)],'ANN Training Data','ANN Test Data','Institutional Models','Physical Benchmarks','location','ne'); end;

end



