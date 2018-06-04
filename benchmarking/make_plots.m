clear all; close all; clc
iFig = 0;

% grab colors
h = plot(randn(7));
for i = 1:7; colors(i,:) = h(i).Color; end
close all

%% *** Initialize Plotting Tools ******************************************



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
    {'Noah.2.7.1'}
    {'Noah.3.2'}
    {'Noah.3.3'}
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

% prepare ann data
annMI = zeros(Ny,2);
for y = 1:Ny
    annMI(y,1) = localStats(y).all.ann.(statsMark{Nstats});
    annMI(y,2) = globalStats(y).all.ann.(statsMark{Nstats});
end % y-loop

% prepare model data
modelMI = zeros(Ny,Nmodels);
for y = 1:Ny
    for m = 1:Nmodels
        modelMI(y,m) = globalStats(y).all.model(m).(statsMark{Nstats});
    end
end % y-loop

%% *** Calculate Missing Info Fractions ***********************************

% info missing from models
Um = repmat(annMI(:,1),[1,Nmodels]) - modelMI;

% info missing from forcings (bounded)
Ud = 1 - repmat(annMI(:,1),[1,Nmodels]);

% total missing info
Ut = 1 - modelMI;

% fractions
Umf = Um./Ut;
Udf = Ud./Ut;

% screen report
fprintf('--- Fraction of Unexplained Entropy due to Data Error ----------\n');
disp([min(Udf')',max(Udf')'])
fprintf('--- Fraction of Unexplained Entropy due to Model Error ---------\n');
disp([min(Umf')',max(Umf')'])

% HILF 
HILF = (annMI(:,1)-annMI(:,2))./annMI(:,1)

%% *** Plot All Benchmark Results *****************************************

% set up figure
iFig = iFig+1; figure(iFig); close(iFig); figure(iFig);
set(gcf,'color','w');
set(gcf,'position',[2121         821          824         433]);

% plot ann data
b = bar(annMI); hold on;
b(1).FaceColor = 0.8*ones(3,1);
b(2).FaceColor = 0.3*ones(3,1);

% plot model data
for m = 1:Nmodels
    plot((1:3)-0.125,modelMI(:,m),'--',...
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
leg = legend([{'Local Benchmark'};{'Global Benchmark'};modelNames],'location','se');
leg.FontSize = 9;
    
% save figure
figure(iFig);
set(gcf,'PaperPositionMode','auto')
saveas(gcf,'figures/Figure 2 - Benchmark Results.png');

%% *** Make Plumber-Type Plots ********************************************

% % stats gorupings used by Best et al. (2015)
% I1 = [2,3,4,7];
% I2 = [5,6];
% I3 = [8,9,10];
% 
% % stats ratios plotting nonsense
% updn = [-1,-1,1,-1,-1,-1,0,-1,-1,1,1];
% 
% % calculate rankings
% stat = zeros(4,Nsites,Nmodels-2,Ny,Nstats)/0;
% rank = zeros(4,Nsites,Nmodels-2,Ny,Nstats)/0;
% group = zeros(4,Ny,Nmodels,3);
% for m = 1:Nmodels-2
%     for y = 1:Ny
%         for s = 1:Nstats
%             for i = 1:Nsites
%                 
%                 % concatenate stats
%                 stat(:,i,m,y,s)  = [globalStats(y).site(i).ann.(statsMark{s}),...
%                     globalStats(y).site(i).model(m).(statsMark{s}),...
%                     globalStats(y).site(i).model(Nmodels-1).(statsMark{s}),...
%                     globalStats(y).site(i).model(Nmodels).(statsMark{s})];
%                 
%                 if updn(s) ~= 0
%                     stat(:,i,m,y,s) = -updn(s)*stat(:,i,m,y,s);
%                 else
%                     stat(:,i,m,y,s) = abs(stat(:,i,m,y,s));
%                 end
%                 
%                 % calculate ranks
%                 [~,rank(:,i,m,y,s)] = sort(stat(:,i,m,y,s));
%             end % i-loop
%         end % s-loop
%         
%         % average rankings
%         for r = 1:4
%             a = squeeze(rank(r,:,m,y,I1)); group(r,y,m,1) = mean(a(:));
%             a = squeeze(rank(r,:,m,y,I2)); group(r,y,m,2) = mean(a(:));
%             a = squeeze(rank(r,:,m,y,I3)); group(r,y,m,3) = mean(a(:));
%         end % r-loop
%         
%     end % y-loop
% end % m-loop
% 
% % remove where models don't predict NEE
% group(3:4,3,:,:) = 0/0;
% group(1:2,[3:6,8:12],:,:) = 0/0;
% 
% % plot
% for s = 1:Nstats
%     
%     % set up figure
%     iFig = iFig+1; figure(iFig); close(iFig); figure(iFig);
%     set(gcf,'color','w');
%     set(gcf,'position',[2024         131        1154        1246]);
%     
%     for y = 1:Ny-1
%         
%         % subplot
%         subplot(Ny-1,1,y);
%         
%         % plot
% %         bar(squeeze(group(:,y,1:Nmodels-2,g))')
%         bar(squeeze(mean(rank(:,:,1:Nmodels-2,y,s),2))');
%         
%         % labels
%         set(gca,'fontsize',16);
%         ylabel('Average Rank','fontsize',20);
%         title(strcat(statNames{s},' -- ',targNames{y}),'fontsize',20);
%         set(gca,'xticklabel',modelNames);
%         xtickangle(gca,60)
%         
%     end %y-loop
% end % g-loop

%% *** END PROGRAM ********************************************************