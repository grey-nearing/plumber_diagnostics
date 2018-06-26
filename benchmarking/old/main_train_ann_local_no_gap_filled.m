clear all; close all; clc
restoredefaultpath; addpath('../matlab_tools');

%% *** Experiment Setup ***************************************************

% target names
targNames = [{'Qe'},{'Qh'},{'NEE'}];

% site names
siteNames = [{'Amplero'},{'Blodgett'},{'Bugac'   },{'ElSaler2'},{'ElSaler' },...
    {'Espirra'},{'FortPeck'},{'Harvard' },{'Hesse'   },{'Howlandm'},...
    {'Howard' },{'Hyytiala'},{'Kruger'  },{'Loobos'  },{'Merbleue'},...
    {'Mopane' },{'Palang'  },{'Sylvania'},{'Tumba'   },{'UniMich' }];
Ns = length(siteNames);

% number of data to use from each site
Nt = 3e4;

% # k-fold validation segments
kfold = 3;

% number of mutual info bins
Nbins = 100;

% ann training parameters
trainParms.verbose       = 0;
trainParms.nodes         = 20;
trainParms.trainRatio    = 0.65;
trainParms.max_fail      = 100;
trainParms.epochs        = 1e4;
trainParms.trainFcn      = 'trainscg';
trainParms.performFcn    = 'mse';

%% *** Load Data **********************************************************

% load fluxnet data
fname = strcat('../data/lagged_data/lagged_fluxnet_',siteNames{1},'_qfc.mat');
load(fname);
data.flux = data.flux(:,1:3);
data.forc = data.forc(:,:,:);

% load models data
fname = strcat('../data/lagged_data/lagged_models_',siteNames{1},'_qfc.mat');
load(fname);
model = model(:,10:12,:);

% dimensions
[~,Nx,Nl] = size(data.forc);
Ny = size(data.flux,2);
Nd = size(data.date,2);
Nm = size(model,3);

% init storage
Iwet = zeros(Nt,Ns)/0;
Idry = zeros(Nt,Ns)/0;
Nwet = zeros(Ns,1)/0;
Ndry = zeros(Ns,1)/0;

date = zeros(Nt,3,Ns);
forc = zeros(Nt,Nx,Nl,Ns);
flux = zeros(Nt,Ny,Ns);
modl = zeros(Nt,Ny,Nm,Ns);

% screen report
fprintf('Loading Data ... \n'); tic;
fprintf(repmat('.',[1,Ns]));
fprintf('\n');

% load training/test data from all sites
for s = 1:Ns
    fprintf('.');
    
    % load fluxnet data
    fname = strcat('../data/lagged_data/lagged_fluxnet_',siteNames{s},'_qfc.mat');
    load(fname);
    data.flux = data.flux(:,1:3);
    data.forc = data.forc(:,:,:);
    
    % load models data
    fname = strcat('../data/lagged_data/lagged_models_',siteNames{s},'_qfc.mat');
    load(fname);
    model = model(:,10:12,:);
    
    % check dimensions
    assert(length(model) == length(data.date));
    assert(Nt <= length(data.date));
    
    % take random sample
    Idex = randsample(length(data.date),Nt);
    
    % store data in structures
    date(:,:,s)   = data.date(Idex,:);
    forc(:,:,:,s) = data.forc(Idex,:,:);
    flux(:,:,s)   = data.flux(Idex,:);
    modl(:,:,:,s) = model(Idex,:,:);
    
end % s-loop

% error check
assert(isempty(find(isnan(flux),1)))
assert(isempty(find(isnan(forc),1)))
assert(isempty(find(isnan(date),1)))

% screen report
fprintf(' finished; time = %f \n',toc);

%% *** Pre-Processing *****************************************************

% mutual info bins
Bw = zeros(Ny,1);
for y = 1:Ny
    yy = flux(:,y,:);
    Bmin = min(yy(:))-1e-6;
    Bmax = max(yy(:))+1e-6;
    By = linspace(Bmin,Bmax,Nbins);
    Bw(y) = By(2) - By(1);
end

% kfold partitions for site-specific runs
assert(rem(Nt,kfold)==0);
Ntrn = Nt - Nt/kfold;
Ktst = reshape(randperm(Nt),[],kfold);      % testing points
Ktrn = zeros(Ntrn,kfold)/0;                 % training points
for k = 1:kfold
    Ktrn(:,k) = setdiff(1:Nt,Ktst(:,k));
end
Ntst = size(Ktst,1);
assert(Ntrn+Ntst==Nt);

% collapse lagged dimension
Xsite = zeros(Nt,Nx*Nl+1,Ns   )/0;
Ysite = zeros(Nt,Ny     ,Ns   )/0;
Msite = zeros(Nt,Ny     ,Ns,Nm)/0;
for s = 1:Ns
    
    % regressors
    temp = squeeze(forc(:,:,:,s));
    temp = permute(temp,[1,3,2]);
    Xsite(:,2:end,s) = reshape(temp,[Nt,Nl*Nx]);
    Xsite(:,1,s) = date(:,3,s);
    
    % regressands
    Ysite(:,:,s) = flux(:,:,s);
    
    % predictions
    Msite(:,:,s,:) = modl(:,:,:,s);
    
end

% indexes of the precip input (to remove for wet/dry)
pdex = 1+(Nx-1)*Nl+1;
ddex = 1:Nx*Nl+1; ddex(pdex) = [];

% error check
assert(isempty(find(isnan(Xsite(:)),1)));
assert(isempty(find(isnan(Ysite(:)),1)));

%% *** Train K-Fold Site-Specific ANN Models ******************************

% init sorage - predictions
Zsite = zeros(Nt,Ny,Ns,2)./0;

% train/test loop
for s = 1:Ns    % loop through sites
    for y = 1:Ny        % loop through output dimensions
        
        % segregate site data
        Xa = Xsite(:,:,s);
        Ya = Ysite(:,y,s);
        Ma = squeeze(Msite(:,y,s,:));
        Za = zeros(length(Ya),2)/0;
        
        % skip if there is no data of this type for this site
        if all(isnan(Ya)); continue; end
        
        % remove grandmas
        assert(isempty(find(isnan(Xa(:)),1)));
        assert(isempty(find(isnan(Ya(:)),1)));
        
        % train ALL model
        for k = 1:kfold
            
            % screen report
            fprintf('Training Local ANN-ALL: Site = %s; Target = %s; K-fold = %d/%d ...',...
                siteNames{s},targNames{y},k,kfold); tic;
            
            % segregate k-fold training/test
            Xtrn = Xa(Ktrn(:,k),:); Ytrn = Ya(Ktrn(:,k),:);
            Xtst = Xa(Ktst(:,k),:); Ytst = Ya(Ktst(:,k),:);
            
            % train the network on segregated data
            ann{y,s,1} = trainANN(Xtrn,Ytrn,trainParms);
            
            % predict on segregated data
            Za(Ktst(:,k),1) = ann{y,s,1}(Xtst');
            
            % screen report
            fprintf('. finished; time = %f \n',toc);
            
        end % k-loop
        
        %         % train wet/dry model
        %         for k = 1:kfold
        %
        %             % segregate k-fold training/test
        %             Xtrn = Xa(Ktrn(:,k),:); Ytrn = Ya(Ktrn(:,k),:);
        %             Xtst = Xa(Ktst(:,k),:); Ytst = Ya(Ktst(:,k),:);
        %
        %             % find wet training points
        %             Iwtrn = Xtrn(:,pdex) > 0;
        %             Xwtrn = Xtrn(Iwtrn,:); Ywtrn = Ytrn(Iwtrn);
        %
        %             % find wet test points
        %             Iwtst = Xtst(:,pdex) > 0;
        %             Xwtst = Xtst(Iwtst,:); Ywtst = Ytrn(Iwtrn);
        %
        %             % find dry training points
        %             Idtrn = Xtrn(:,pdex) <= 0;
        %             Xdtrn = Xtrn(Idtrn,ddex); Ydtrn = Ytrn(Idtrn);
        %
        %             % find dry test points
        %             Idtst = Xtst(:,pdex) <= 0;
        %             Xdtst = Xtst(Idtst,ddex); Ydtst = Ytrn(Idtrn);
        %
        %             % error checking
        %             assert(sum(Iwtrn)+sum(Idtrn) == Ntrn);
        %             assert(all(~isequal(Iwtrn,Idtrn)));
        %             assert(sum(Iwtst)+sum(Idtst) == Ntst);
        %             assert(all(~isequal(Iwtst,Idtst)));
        %
        %             % train the network on segregated data
        %             fprintf('Training Local ANN-WET: Site = %s; Target = %s; K-fold = %d/%d ...',...
        %                 siteNames{s},targNames{y},k,kfold); tic;
        %             ann{y,s,2} = trainANN(Xwtrn,Ywtrn,trainParms);
        %             Za(Ktst(Iwtst,k),2) = ann{y,s,2}(Xwtst');
        %             fprintf('. finished; time = %f \n',toc);
        %
        %             % screen report
        %             fprintf('Training Local ANN-DRY: Site = %s; Target = %s; K-fold = %d/%d ...',...
        %                 siteNames{s},targNames{y},k,kfold); tic;
        %             ann{y,s,3} = trainANN(Xdtrn,Ydtrn,trainParms);
        %             Za(Ktst(Idtst,k),2) = ann{y,s,3}(Xdtst');
        %             fprintf('. finished; time = %f \n',toc);
        %
        %         end % k-loop
        
        % store result
        Zsite(:,y,s,1) = Za(:,1);
        Zsite(:,y,s,2) = Za(:,2);
        
        % calculate ann-all stats
        stats(y).site(s).ann(1) = calcStats(Ya,Za(:,1),Bw);
        
        %         % calculate ann-wetdry stats
        %         stats(y).site(s).ann(2) = calcStats(Ya,Za(:,2),Bw);
        
        % calculate model stats
        for m = 1:Nm
            stats(y).site(s).model(m) = calcStats(Ya,Ma(:,m),Bw);
        end % m-loop
        
        % calculate improvement fractions
        fn = fieldnames(stats(y).site(s).ann(1));
        ratio = zeros(Nm+1,length(fn))/0;
        updn = [-1,1,-1,-1,-1,0,0,0,1];
        %         ratio(1,1) = -(stats(y).site(s).ann(1).rmse - stats(y).site(s).ann(2).rmse) ...
        %             ./ stats(y).site(s).ann(2).rmse;
        %         ratio(1,2) = (stats(y).site(s).ann(1).mi - stats(y).site(s).ann(2).mi) ...
        %             ./ stats(y).site(s).ann(2).mi;
        for m = 1:Nm
            for f = 1:length(fn)
                ratio(m+1,f) = updn(f)*(stats(y).site(s).ann(1).(fn{f}) - stats(y).site(s).model(m).(fn{f})) ...
                    ./ stats(y).site(s).model(m).(fn{f});
                if updn(f) == 0
                    ratio(m+1,f) = (abs(stats(y).site(s).ann(1).(fn{f})) - abs(stats(y).site(s).model(m).(fn{f}))) ...
                        ./ abs(stats(y).site(s).model(m).(fn{f}));
                end
            end
        end % m-loop
        
        % print to screen
        fprintf('--- Site-Specific Improvement Results --------------------------- \n')
        disp(ratio);
        fprintf('----------------------------------------------------------------- \n\n')
        
    end % y-loop
    
% save progress
fname = './results/progress_local_anns_qfc.mat';
save(fname,'-v7.3');

end % s-loop

%% *** Site-Specific Model Evaluation *************************************

% all-site data
Yall = reshape(permute(Ysite,[1,3,2]),[Nt*Ns,Ny]);

Zall = zeros(Nt*Ns,Ny,2)/0;
for z = 1:2
    Zall(:,:,z) = reshape(permute(squeeze(Zsite(:,:,:,z)),[1,3,2]),[Nt*Ns,Ny]);
end % z-loop

Mall = zeros(Nt*Ns,Ny,Nm)/0;
for m = 1:Nm
    Mall(:,:,m) = reshape(permute(squeeze(Msite(:,:,:,m)),[1,3,2]),[Nt*Ns,Ny]);
end % m-loop

% global statistics
for y = 1:Ny
    
    % calculate ann-all stats
    stats(y).all.ann(1) = calcStats(Yall(:,y),Zall(:,y,1));
    
    %     % calculate ann-wetdry stats
    %     stats(y).all.ann(2) = calcStats(Yall(:,y),Zall(:,y,2));
    
    % calculate model stats
    for m = 1:Nm
        stats(y).all.model(m) = calcStats(Yall(:,y),Mall(:,y,m));
    end % m-loop
    
    %     ratio(1,1) = -(stats(y).site(s).ann(1).rmse - stats(y).site(s).ann(2).rmse) ...
    %         ./ stats(y).site(s).ann(2).rmse;
    %     ratio(1,2) = (stats(y).site(s).ann(1).mi - stats(y).site(s).ann(2).mi) ...
    %         ./ stats(y).site(s).ann(2).mi;
    for m = 1:Nm
        for f = 1:length(fn)
            ratio(m+1,f) = updn(f)*(stats(y).site(s).ann(1).(fn{f}) - stats(y).site(s).model(m).(fn{f})) ...
                ./ stats(y).site(s).model(m).(fn{f});
            if updn(f) == 0
                ratio(m+1,f) = (abs(stats(y).site(s).ann(1).(fn{f})) - abs(stats(y).site(s).model(m).(fn{f}))) ...
                    ./ abs(stats(y).site(s).model(m).(fn{f}));
            end
        end
    end % m-loop
    % print to screen
    fprintf('--- Site-Specific Improvement Results --------------------------- \n')
    disp(ratio);
    fprintf('----------------------------------------------------------------- \n\n')
    
end % y-loop

%% *** Save Results *******************************************************

% save progress
fname = './results/progress_local_anns_qfc.mat';
save(fname,'-v7.3');

% save final statistics only
fname = './results/stats_local_anns_qfc.mat';
save(fname,'stats','-v7.3');

%% *** END PROGRAM ********************************************************