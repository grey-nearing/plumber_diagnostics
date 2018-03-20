clear all
close all
clc

% raw data
load('results/tradeoff_results.mat');

% dimensions
[Nx,Ny,Ns,Nm] = size(TD);

% model groups
Ng = 8;
Imm = zeros(Ng,Nm);
Imm(1,[1,2])      = 1;
Imm(2,[3])        = 1;
Imm(3,[4])        = 1;
Imm(4,[5,6])      = 1;
Imm(5,[7,8])      = 1;
Imm(6,[9])        = 1;
Imm(7,[10,11,12]) = 1;
Imm(8,[13,14])    = 1;

% init storage
Ai = zeros(Nx,Ny)/0;
Si = zeros(Nx,Ny,Ns)/0;
Mi = zeros(Nx,Ny,Ng)/0;
Sc = zeros(Nx,Ny)/0;
Mc = zeros(Nx,Ny)/0;

for x = 1:Nx
 for y = 1:Ny

  % --- all stats ----------------------------------------------------------

  % all cluster center
  td = squeeze(TD(x,y,:,:)); 
  mi = squeeze(MI(y,:,:));   

  % remove any bad site
  Is = find(any(~isnan([td,mi]'))); 
  J = find(all(isnan([td,mi]'))); td(J,:) = []; mi(J,:) = []; 
  assert(isempty(J))

  % remove any bad model
%  Im = find(all(~isnan([td;mi]))); 
  Im = Imm;
  J = find(any(isnan([td;mi]))); td(:,J) = []; mi(:,J) = []; Im(:,J) = [];

  % dimensions
  [Nsites,Nmodels] = size(td);
  [x,y,Nsites,Nmodels]

  % find cluster centers
  Xa = mean(td(:));
  Ya = mean(mi(:));

  % all cluster mean distance
  Ai(x,y) = mean( ((td(:)-Xa).^2 + (mi(:)-Ya).^2).^(1/2) );

  % --- model stats --------------------------------------------------------

  % init storage
  Xm = zeros(Ng,1)/0;
  Ym = zeros(Ng,1)/0;

  % loop through models
  for g = 1:Ng

   % pull model group
   tdg = td(:,find(Im(g,:))); tdg = tdg(:);
   mig = mi(:,find(Im(g,:))); mig = mig(:);

   if isempty(tdg)
    Xm(g) = 0/0;
    Ym(g) = 0/0;
    Mi(x,y,g) = 0/0;
   end

   % find cluster centers
   Xm(g) = mean(tdg);
   Ym(g) = mean(mig);

   % mean distance of each model from their center
   Mi(x,y,g) = mean( ((tdg-Xm(g)).^2 + (mig-Ym(g)).^2).^(1/2) );

  end

  % distance of cluster centers to mean
  Mca(x,y) = mean( ((Xm-Xa).^2 + (Ym-Ya).^2).^(1/2) );

  % mean distance between centers
  Mc(x,y) = 0;
  count = 0;
  for g = 1:Ng
   for gg = 1:g-1
    Mc(x,y) = Mc(x,y) + ((Xm(g)-Xm(gg))^2 + (Ym(g)-Ym(gg))^2)^(1/2);
    count = count + 1;
   end
  end
  Mc(x,y) = Mc(x,y)/count;

  % --- site stats ---------------------------------------------------------

  % init storage
  Xs = zeros(Nsites,1)/0;
  Ys = zeros(Nsites,1)/0;

  % loop through sites
  for s = 1:Nsites

   % find cluster centers
   Xs(s) = mean(td(s,:));
   Ys(s) = mean(mi(s,:));

   % mean distance of each site from their center
   Si(x,y,Is(s)) = mean( ((td(s,:)-Xs(s)).^2 + (mi(s,:)-Ys(s)).^2).^(1/2) );
  end

  % distance of cluster centers to mean
  Sca(x,y) = mean( ((Xs-Xa).^2 + (Ys-Ya).^2).^(1/2) );

  % mean distance between centers
  Sc(x,y) = 0;
  count = 0;
  for s = 1:Nsites
   for ss = 1:s-1
    Sc(x,y) = Sc(x,y) + ((Xs(s)-Xs(ss))^2 + (Ys(s)-Ys(ss))^2)^(1/2);
    count = count + 1;
   end
  end
  Sc(x,y) = Sc(x,y)/count;

 end % y
end % x

%nanmean((repmat(Ai,[1,1,Ng])-Mi)./repmat(Ai,[1,1,Ng]),3)
%nanmean((repmat(Ai,[1,1,Ns])-Si)./repmat(Ai,[1,1,Ns]),3)
%
%(nanmean(Mi,3)-nanmean(Si,3))./nanmean(Si,3)
%(nanmean(Mi,3)-nanmean(Si,3))./Ai
%
%(Mc-Sc)./Ai
%
%Mca-Sca

%  
figure(4);close(4);figure(4);
set(gcf,'color','w');

%bar((nanmean(Mi,3)-nanmean(Si,3))./nanmean(Si,3))
subplot(2,1,1)
bar(nanmean(Mi,3)./Ai);
grid on;
legend(varNames{6:8},'location','sw');
set(gca,'xticklabel',varNames)
title('mean dist to model centers','fontsize',14);
set(gca,'ylim',[0,1]);

subplot(2,1,2)
bar(nanmean(Si,3)./Ai);
grid on;
set(gca,'xticklabel',varNames)
title('mean dist to site centers','fontsize',14);
set(gca,'ylim',[0,1]);

% save sites figure
figure(4);
set(gcf,'PaperPositionMode','auto')
saveas(gcf,'figures/Clusters.png');





