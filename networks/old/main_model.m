clear all
close all
clc
restoredefaultpath
addpath('../matlab_tools');

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

% site names
siteNames =  [{'Amplero'}
             {'Blodgett'}
             {'Bugac'}
             {'ElSaler2'}
             {'ElSaler'}
             {'Espirra'}
             {'FortPeck'}
             {'Harvard'}
             {'Hesse'}
             {'Howlandm'}
             {'Howard'}
             {'Hyytiala'}
             {'Kruger'}
             {'Loobos'}
             {'Merbleue'}
             {'Mopane'}
             {'Palang'}
             {'Sylvania'}
             {'Tumba'}
             {'UniMich'}];

% number of sites
Ns = length(siteNames);
Nm = length(modelNames);

% variable names
varNames = [{'Ws'},{'Ta'},{'Rn'},{'Rh'},{'PP'},{'Qe'},{'Qh'},{'NEE'},{'SM1'},{'SM2'}];
Du = 5; % number of forcing data
D = 10;

% number of histogram bins
Nb = 30;

% temporal lags
%lags = round(logspace(log10(1),log10(48*180),15));
lags = [1,2,6,12,24,48,96,128,256];
%lags = [1,48];
Nl = length(lags);

% loop through sites
for s = 1:Ns

 % screen report
 fprintf('Working on site %d of %d (%s) ...',s,Ns,siteNames{s}); tsite = tic;

 % load data
 fname = strcat('../data/model_data/extracted/',siteNames{s},'.mat');
 load(strcat(fname));
 data = model; clear model;
 dates = data(:,1:3,1);
 data(:,1:3,:) = [];  % remove dates

 % find start and end dates of actual data
 If = find(~isnan(dates(:,1)),1,'first');
 Il = find(~isnan(dates(:,1)),1,'last');
 dates = dates(If:Il,:); data = data(If:Il,:,:);

 % find discontinuous missing data 
 imiss = isempty(find(isnan(dates(:)),1,'first'));
 count = 0;
 while ~imiss
  If = find(isnan(dates(:,1)),1,'first');
  if If>0.5*length(dates(:,1));
   dates(If:end,:) = []; data(If:end,:,:) = []; 
  else
   dates(1:If,:) = []; data(1:If,:,:) = []; 
  end

  count = count+1;
  if count > 1; error('strange patterns of missing data: (%d,%d)',If,Il); end;
  imiss = isempty(find(isnan(dates(:)),1,'first'));
 end

 % data dimensions
 [Nt,Dz,M] = size(data);
 assert(D==Dz);
 assert(M==Nm);

 % screen report
 if count == 0 
  fprintf('. Ndata = %d -- Any Missing? No. ',Nt);
 else
  fprintf('. Ndata = %d -- Any Missing? Yes!! ',Nt);
 end

 % loop through models
 for m = 1:Nm

  % separate data for individual model
  mdata = squeeze(data(:,:,m));

  % histogram bin sizes
  for d = 1:D
   Bmin = min(mdata(:,d))-1e-6;
   Bmax = max(mdata(:,d))+1e-6;
   B(:,d) = linspace(Bmin,Bmax,Nb); 
  end

  % process networks at the different scales
  for l = 1:Nl

   % init storage
   T = zeros(D)./0;
   H = zeros(D)./0;
   S = zeros(D)./0;

   % input data loop at site/model
   for x = 1:D

    % skip precip
%    if x == Du; continue; end;

    % grab raw data
    XXX = mdata(:,x); XXX(abs(XXX-9999)<1) = 0/0;

    % find start and end dates of actual data
    Ifx = find(~isnan(XXX),1,'first');
    Ilx = find(~isnan(XXX),1,'last');
    XXX = XXX(Ifx:Ilx,:); 

    % deal with grandmas
    try
     XXX = grandma_smoothing(XXX,lags(l));
    catch
     continue
    end

    % if data is constant, skip
    if isempty(XXX); continue; end;
    if max(abs(diff(XXX)))==0; continue; end;

    % target data loop at site/model
    for y = Du+1:D

     % don't target self
     if x==y; continue; end;

     % raw target data
     YY = mdata(:,y); YY(abs(YY-9999)<1) = 0/0;
     YY = YY(Ifx:Ilx,:);

     % find start and end dates of actual data
     Ify = find(~isnan(YY),1,'first');
     Ily = find(~isnan(YY),1,'last');
     XX = XXX(Ify:Ily,:); 
     YY = YY(Ify:Ily,:); 

     % deal with grandma
     try
      YY = grandma_smoothing(YY,lags(l));
     catch
      continue
     end
    
     % if data is constant, skip
     if isempty(XX); continue; end;
     if max(abs(diff(XX)))==0; continue; end;
     if isempty(YY); continue; end;
     if max(abs(diff(YY)))==0; continue; end;
 
     % deal with lags
     Xw = window_average(XX,lags(l));
     Yw = window_average(YY,lags(l));
     assert(numel(Xw)==numel(Yw));
     if numel(Xw)<500; continue; end;

     % time shift
     Xt = Xw(1:end-1,:); 
     Yt = Yw(1:end-1,:); 
     Ys = Yw(2:end,:); 

     % make sure we have dealt with missing values
     assert(isempty(find(isnan(Xt),1,'first')));
     assert(isempty(find(isnan(Yt),1,'first')));
     assert(isempty(find(isnan(Ys),1,'first')));

     % do the actaul calculations
     [T(x,y),H(x,y),S(x,y)] = transfer_entropy_window_average(Ys,Xt,Yt,B(:,x),B(:,y));

     % screen report
     fprintf('\n    Site: %s - Model: %s - Lag = %d - %s -> %s = %f - Ndata = %d',siteNames{s},modelNames{m},lags(l),varNames{x},varNames{y},T(x,y)./H(x,y),length(YY));

    end % x
   end % y

   % save in Circos format
   fname = strcat('results/models/MODEL_',siteNames{s},'_',modelNames{m},'_',num2str(l),'.txt');
   saveCircos(T,varNames,Du,D,fname);

   % save in Circos format
   TT = T./H; 
   fname = strcat('results/models/MODEL_Normalized_',siteNames{s},'_',modelNames{m},'_',num2str(l),'.txt');
   saveCircos(TT,varNames,Du,D,fname);

  end % lags

  % save progress
  save('results/models/saveModelProgress.mat');

 end % model

 % screen report
 t = toc(tsite); fprintf('\n');
 fprintf('Finished Site %s - time = %f\n',siteNames{s},t);

end




