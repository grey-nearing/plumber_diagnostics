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
varNames = [{'Ws'},{'Ta'},{'Rn'},{'Rh'},{'PP'},{'Qe'},{'Qh'},{'NEE'},{'SM'}];%,{'SM2'}];
Du = 5; % number of forcing data
Dx = length(varNames);

% number of histogram bins
Nb = 30;

% temporal lags
lags = [1,2,6,12,24,48,96];
Nl = length(lags);

% init storage
T = zeros(Dx,Dx,Nl,Nm+2,Ns)./0;
H = zeros(Dx,Dx,Nl,Nm+2,Ns)./0;
S = zeros(Dx,Dx,Nl,Nm+2,Ns)./0;

% load data
for s = 1:Ns; tsite = tic;

 % screen report
 fprintf('Working on site %d of %d (%s) ...',s,Ns,siteNames{s}); tsite = tic;

 % load model data
 fname = strcat('../data/model_data/extracted/',siteNames{s},'.mat');
 load(strcat(fname));
 Mdates = model(:,1:3,1);
 model(:,1:3,:) = [];  
 model(:,end,:) = [];
 model(model<=-999) = 0/0;
 assert(size(model,3) == Nm);
 assert(size(model,2) == Dx);

 % load pals data
 fname = strcat('../data/pals_data/extracted/',siteNames{s},'.txt');
 pals = load(strcat(fname));
 Pdates = pals(:,1:3);
 pals(:,1:3) = [];  
 pals(:,end) = [];
 pals(pals<=-999) = 0/0;
 assert(size(pals,2) == Dx);

 % find start and end dates of model data
 Ifm = find(~isnan(Mdates(:,1)),1,'first');
 Ilm = find(~isnan(Mdates(:,1)),1,'last');
 Mdates = Mdates(Ifm:Ilm,:); model = model(Ifm:Ilm,:,:);

 % find start and end dates of pals data
 Ifp = find(~isnan(Pdates(:,1)),1,'first'); assert(Ifm == Ifp);
 Ilp = find(~isnan(Pdates(:,1)),1,'last');  assert(Ilm == Ilp);
 Pdates = Pdates(Ifp:Ilp,:); pals = pals(Ifp:Ilp,:,:);

 % make sure the files have the same dates
 for d = 1:3
  assert(max(abs(Mdates(:,d)-Pdates(:,d)))==0);
 end
 dates = Pdates; clear Mdates Pdates

 % missing years at Espirra 
 if s == 6 
  dates(52609:end,:)   = [];
  pals( 52609:end,:)   = [];
  model(52609:end,:,:) = [];
 end

 %  soil indexes
 I = isnan(pals(:,9))';
 first = find([1,diff(I)]);
 runs = diff([first,length(I)+1]);
 [~,A] = sort(runs,'descend');
 soilstart = [];
 soilend = [];
 for a = 1:length(runs)
  if ~isnan(pals(first(A(a)),9))
   soilstart = first(A(a));
   soilend   = first(A(a)) + runs(A(a))-1;
   break;
  end
 end
 assert(all(~isnan(pals(soilstart:soilend,9))));
 Nsoil = soilend-soilstart+1;

 if Nsoil > 17520
  dates = dates(soilstart:soilend,:);
  pals = pals(soilstart:soilend,:);
  model = model(soilstart:soilend,:,:);
 else
  Nsoil = 0;
  startsoil = []; 
  endsoil = [];
  pals(:,Dx) = [];
  model(:,Dx,:) = [];
 end

 % deal with grandma
 for d = 1:size(pals,2)
  [pals(:,d),Nmissing] = grandma_smoothing(pals(:,d),0.001);
  if all(isnan(pals(:,d))); fprintf('bad pals vector from grandma smoothing: %s. \n',varNames{d}); end;
  for m = 1:Nm
%   if m == 3 
%    [model(1:end-46,d,m),~] = grandma_smoothing(model(1:end-47,d,m),0.001,1);
%    model(end-47:end,d,m) = 0/0;
%   else
    [model(:,d,m),~] = grandma_smoothing(model(:,d,m),0.001,0);
%   end
  end
 end

% for d = Dx
%
%  if Nsoil > 17520
%   [pals(soilstart:soilend,d),~] = grandma_smoothing(pals(soilstart:soilend,d),0.01);
%   if all(isnan(pals(soilstart:soilend,d))); fprintf('bad soils vector from grandma smoothing: %s. \n',varNames{d}); end;
%
%   for m = 1:Nm
%    [model(soilstart:soilend,d,m),~] = grandma_smoothing(model(soilstart:soilend,d,m),0.01,0);
%   end
%
%   if soilstart > 1
%    pals(1:soilstart-1,d) = 0/0;
%    model(1:soilstart-1,d,:) = 0/0;
%   end
%   if soilend < size(pals,1)
%    pals(soilend+1:end,d) = 0/0;
%    model(soilend+1:end,d,:) = 0/0;
%   end
%
%  else
%   pals(:,d) = 0/0;
%   model(:,d,:) = 0/0;
%   Nsoil = 0;
%  end
%
% end

 % data dimensions
 assert(all(size(squeeze(model(:,:,1)))==size(pals)));
 assert(isempty(find(isnan(pals))));
 [Nt,Dp] = size(pals);

 % screen report
 fprintf('. Ndata = %d -- Nmissing = %d -- Nsoil = %d \n',Nt,Nmissing,Nsoil);

 % loop through lags
 for l = 1:Nl

  % lag data 
  clear pals_lag;
  for x = 1:Dp
   pals_lag(:,:,x) = window_average(pals(:,x),lags(l)); 
   if size(pals_lag,1) < 1000; break; end;
  end
  if size(pals_lag,1) < 1000; continue; end;

  % loop through dpn pathways
  for y = Du+1:Dp

   % ann training input data 
   Xtrn = zeros((size(pals_lag,1)-1)*size(pals_lag,2),Dp)./0;
   for x = 1:Du
    xx = pals_lag(2:end,:,x);
    Xtrn(:,x) = xx(:);
   end 
   for x = Du+1:Dp
    xx = pals_lag(1:end-1,:,x);
    Xtrn(:,x) = xx(:);
   end 

   % ann training target
   Ytrn = pals_lag(2:end,:,y);
   Ytrn = Ytrn(:);

   % check grandma
   assert(isempty(find(isnan(Xtrn),1,'first')));
   assert(isempty(find(isnan(Ytrn),1,'first')));
   assert(isempty(find(Xtrn<=-999,1,'first')));
   assert(isempty(find(Ytrn<=-999,1,'first')));

   % train ann
   Yann = ann_train_pred(Xtrn,Ytrn,1:size(Xtrn,1));

   for x = 1:Dp

    % perform on pals - no ann
    Bmin = min(min(Ytrn),min(Xtrn(:,y)))-1e-6; Bmax = max(max(Ytrn),max(Xtrn(:,y)))+1e-6; By = linspace(Bmin,Bmax,Nb); 
    Bmin = min(Xtrn(:,x))-1e-6;                Bmax = max(Xtrn(:,x))+1e-6;                Bx = linspace(Bmin,Bmax,Nb); 
    if x ~= y
     [T(x,y,l,1,s),H(x,y,l,1,s),~] = transfer_entropy(Ytrn,Xtrn(:,x),Xtrn(:,y),Bx,By);
    else
     Yrnd = rand(size(Xtrn(:,y))) * (max(Xtrn(:,y))-min(Xtrn(:,y))) + min(Xtrn(:,y));
     [T(x,y,l,1,s),H(x,y,l,1,s),~] = transfer_entropy(Ytrn,Xtrn(:,x),Yrnd,Bx,By);
    end

    % perform on pals - with ann
    Bmin = min(min(Yann),min(Xtrn(:,y)))-1e-6; Bmax = max(max(Yann),max(Xtrn(:,y)))+1e-6; By = linspace(Bmin,Bmax,Nb); 
    Bmin = min(Xtrn(:,x))-1e-6;                Bmax = max(Xtrn(:,x))+1e-6;                Bx = linspace(Bmin,Bmax,Nb); 
    if x ~= y
     [T(x,y,l,2,s),H(x,y,l,2,s),~] = transfer_entropy(Yann,Xtrn(:,x),Xtrn(:,y),Bx,By);
    else
     Yrnd = rand(size(Xtrn(:,y))) * (max(Xtrn(:,y))-min(Xtrn(:,y))) + min(Xtrn(:,y));
     [T(x,y,l,2,s),H(x,y,l,2,s),~] = transfer_entropy(Yann,Xtrn(:,x),Yrnd,Bx,By);
    end

    % perform on models
    for m = 1:Nm; tmodel = tic;

     I = find(all(~isnan(model(:,[x,y],m))'));
     if length(I) < 0.9*length(pals(:,x)); continue; end;

     % deal with grandma
     [model(:,x,m),~] = grandma_smoothing(model(:,x,m),0.001);
     [model(:,y,m),~] = grandma_smoothing(model(:,y,m),0.001);

     % lag data 
     X = window_average(model(:,x,m),lags(l)); if numel(X) < 1000; continue; end;
     Y = window_average(model(:,y,m),lags(l)); if numel(Y) < 1000; continue; end;

     % get
     Ys = Y(2:end,:);   Ys = Ys(:);
     Yt = Y(1:end-1,:); Yt = Yt(:);
     if x <= Du 
      Xt = X(2:end,:); 
     else
      Xt = X(1:end-1,:); 
     end;               Xt = Xt(:);

     % check grandma
     assert(isempty(find(isnan(Xt),1,'first')));
     assert(isempty(find(isnan(Yt),1,'first')));
     assert(isempty(find(isnan(Ys),1,'first')));

     % calculate
     Bmin = min(min(Ys),min(Yt))-1e-6; Bmax = max(max(Ys),max(Yt))+1e-6; By = linspace(Bmin,Bmax,Nb); 
     Bmin = min(Xt)-1e-6;              Bmax = max(Xt)+1e-6;              Bx = linspace(Bmin,Bmax,Nb); 
     if x ~= y
      [T(x,y,l,2+m,s),H(x,y,l,2+m,s),~] = transfer_entropy(Ys,Xt,Yt,Bx,By);
     else
      Yrnd = rand(size(Yt)) * (max(Yt)-min(Yt)) + min(Yt);
      [T(x,y,l,2+m,s),H(x,y,l,2+m,s),~] = transfer_entropy(Ys,Xt,Yrnd,Bx,By);
     end

     % screen report
     t = toc(tmodel); 
     fprintf(' ---- M(%d/%d) L(%d/%d); %s -> %s;  %f / %f / %f; time = %f \n',m,Nm,l,Nl,varNames{x},varNames{y},...
        T(x,y,l,2+m,s),T(x,y,l,1,s),T(x,y,l,2,s),t);

    end % m

   end % x
  end % y
 end % l

 % screen report
 t = toc(tsite); fprintf('\n');
 fprintf('Finished Site %s - time = %f\n',siteNames{s},t);

 % save progress
 fname = strcat('results/main_dpn_results_',num2str(s),'.mat');
 save(fname);

end % s

%% save models figure
%figure(1);
%set(gcf,'PaperPositionMode','auto')
%saveas(gcf,'figures/ProcessByModels.png');





