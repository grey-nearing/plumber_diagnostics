clear all; close all; clc
restoredefaultpath; addpath('../matlab_tools');

%% *** Experiment Setup ***************************************************

% target names
targNames = [{'Qe'},{'Qh'},{'NEE'}];

%% *** Load Data **********************************************************

% load stats from gap-filled data
fname = './results/stats_local_anns.mat';
load(fname);
gapStats = stats; 
clear stats

% load stats from non-gap-filled data
fname = './results/stats_local_anns_fqc.mat';
load(fname);
nonStats = stats; 
clear stats

% count stats
fn = fieldnames(gapStats(y).all.ann);
Nstats = length(fn);

%% *** Make Plots *********************************************************

% set up figure
iFig = iFig + 1;
figure(iFig); close(iFig); figure(iFig);
set(gcf,'color','w');
% set(gcf,'position',[1925         403        1316         644]);

% calculate values for plotting
for y = 1:length(targNames)
    vals(y,1) = gapStats(y).all.ann.mi;
    vals(y,2) = nonStats(y).all.ann.mi;
    vals(y) = (gapStats(y).all.ann.mi - nonStats(y).all.ann.mi) ./ nonStats(y).all.ann.mi;    
end % y-loop

% plot stuff
bar(vals)

% aesthetics
set(gca,'yticklabels',[{'gap-filled'},{'non-filled'},{'diff. ratio'}])
legend(targNames);

% save figure
figure(iFig);
set(gcf,'PaperPositionMode','auto')
saveas(gcf,'./figures/FigA2 - gap_fill_comparison.png');

%% *** END PROGRAM ********************************************************