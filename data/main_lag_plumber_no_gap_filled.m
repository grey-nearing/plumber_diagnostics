clear all; close all; clc
restoredefaultpath; addpath(genpath(pwd));

%% *** Model & Site Names *************************************************

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
    {'ORCHIDEE.trunk_r1401'}
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

% site names
siteNames = [{'Amplero'},{'Blodgett'},{'Bugac'   },{'ElSaler2'},{'ElSaler' },...
    {'Espirra'},{'FortPeck'},{'Harvard' },{'Hesse'   },{'Howlandm'},...
    {'Howard' },{'Hyytiala'},{'Kruger'  },{'Loobos'  },{'Merbleue'},...
    {'Mopane' },{'Palang'  },{'Sylvania'},{'Tumba'   },{'UniMich' }];

% specified parameters
Nsites  = length(siteNames);    % number of sites
Nmodels = length(modelNames);   % number of models

%% *** Experiment Dimensions **********************************************

% number of models & sites
Nsites = length(siteNames);
Nmodels = length(modelNames);

% lag parameters
Ltimes      = 4;            % timestep lags
Ldays       = 2;            % day lags
Lweeks      = 2;            % weeks lags
Lmonths     = 2;            % month lags

% derived parameters
Linpt = Ltimes+Ldays*48+Lweeks*48*7+Lmonths*48*30;  % total number of lagged inputs
Louts = Ltimes+Ldays+Lweeks+Lmonths;                % total number of lagged outputs

% data dimensions
data = load('pals_data/extracted/Amplero.txt');
[Tmax,D] = size(data);
Dx = 6; Dy = 5;
assert(D == Dx+Dy+3);

%% *** Gather Data *******************************************************

% loop through sites
for s = 1:Nsites
    
    % init storage
    X = zeros(Tmax,Dx,Louts)./0;
    Y = zeros(Tmax,Dy      )./0;
    D = zeros(Tmax,3       )./0;
    
    Xfqc = zeros(Tmax,1);
    Yfqc = zeros(Tmax,Dy);
    
    % screen report
    fprintf('Site %d of %d ...',s,Nsites); tic;
    
    % load forcing and flux data
    data = load(strcat('pals_data/extracted/',siteNames{s},'.txt'));
    date = data(:,1:3);
    forc = data(:,3+(1:Dx));
    prog = data(:,3+Dx+(1:Dy));
    
    % load forcing and flux quality control files
    data = load(strcat('pals_data/extracted/',siteNames{s},'_fqc.txt'));
    datefqc = data(:,1:3); 
    forcfqc = data(:,3+(1:Dx));
    progfqc = data(:,3+Dx+(1:Dy));
    
    % load model data
    model = load(strcat('model_data/extracted/',siteNames{s},'.mat'));
    model = model.model;

    % check that the file sizes are correct and aligned
    Tfile = size(data,1);
    Tmodl = size(model,1);
    assert(Tfile == Tmax);
    assert(Tmodl == Tmax);

    I = find(model(:,1,1)>-9990);
    assert(max(abs(model(I,1,1)-date(I,1)))==0)
    assert(max(abs(model(I,2,1)-date(I,2)))==0)
    assert(max(abs(model(I,3,1)-date(I,3)))==0)
    
    % store site inputs
    for t = Linpt:Tfile
        
        % store dates
        D(t,:) = date(t,:);
        
        % store timestep lags
        for l = 1:Ltimes
            X(t,:,l) = forc(t-l+1,:);
        end
        
        % store day lags
        for l = 1:Ldays
            edex = t-(Ltimes+48*(l-1));
            sdex = t-(Ltimes+48*l)+1;
            dat = forc(sdex:edex,:); dat = dat(:);
            X(t,:,Ltimes+l) = mean(forc(sdex:edex,:));
            Xfqc(t) = max(max(forcfqc(sdex:edex,:)));
        end
        
        % store weeks lags
        for l = 1:Lweeks
            edex = t-Ltimes-48*Ldays-(48*7*(l-1));
            sdex = t-Ltimes-48*Ldays-(48*7*l)+1;
            dat = forc(sdex:edex,:); dat = dat(:);
            X(t,:,Ltimes+Ldays+l) = mean(forc(sdex:edex,:));
            Xfqc(t) = max(Xfqc(t),max(max(forcfqc(sdex:edex,:))));
        end
        
        % store month lags
        for l = 1:Lmonths
            edex = t-Ltimes-48*Ldays-48*7*Lweeks-(48*30*(l-1));
            sdex = t-Ltimes-48*Ldays-48*7*Lweeks-(48*30*l)+1;
            dat = forc(sdex:edex,:); dat = dat(:);
            X(t,:,Ltimes+Ldays+Lweeks+l) = mean(forc(sdex:edex,:));
            Xfqc(t) = max(Xfqc(t),max(max(forcfqc(sdex:edex,:))));
        end
        
        % store site targets
        Y(t,:) = prog(t,:);
        Yfqc(t,:) = progfqc(t,:);
        
    end
        
    % beginning missing values
    X(1:Linpt-1,:,:) = [];
    Y(1:Linpt-1,:) = [];
    D(1:Linpt-1,:) = [];
    model(1:Linpt-1,:,:) = [];
    
    % remove any lines with missing data
    XX = reshape(X,[size(X,1),size(X,2)*size(X,3)]);
    I = find(any(isnan(XX')));
    X(I,:,:) = []; Y(I,:) = []; D(I,:) = []; model(I,:,:) = [];
    
    % remove any lines with any gap-filled data
    I = find(any(Xfqc > 0));
    X(I,:,:) = []; Y(I,:) = []; D(I,:) = []; model(I,:,:) = [];

    I = find(any(Yfqc > 0));
    X(I,:,:) = []; Y(I,:) = []; D(I,:) = []; model(I,:,:) = [];
    
    % remove any lines with missing data
    XX = reshape(X,[size(X,1),size(X,2)*size(X,3)]);
    I = find(any(XX'<-9990));
    X(I,:,:) = []; Y(I,:) = []; D(I,:) = []; model(I,:,:) = [];
    
    % save lagged data
    clear data
    data.date = D;
    data.forc = X;
    data.flux = Y;
    fname = strcat('lagged_data/lagged_fluxnet_',siteNames{s},'_fqc.mat');
    save(fname,'data');
    
    % save corresponding model data
    fname = strcat('lagged_data/lagged_models_',siteNames{s},'_fqc.mat');
    save(fname,'model');
    
    % screen report
    t = toc; fprintf('Finished: time = %f \n',t);
           
end

%% *** END SCRIPT *********************************************************


