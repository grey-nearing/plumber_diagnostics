function convergence_bootstrap(iType,iTarg)

% master routines
addpath('../matlab_tools');

% time the fucntion call
masterTime = tic;

% maximum number of data
Dmax = 1e2*2^11;

% training configuration
Nepochs = 10000;
trnfctn = 'trainscg';
%trnfctn = 'traincgb';
%trnfctn = 'trainbr';

% site names
siteNames = [{'Amplero'},{'Blodgett'},{'Bugac'   },{'ElSaler2'},{'ElSaler' },...
             {'Espirra'},{'FortPeck'},{'Harvard' },{'Hesse'   },{'Howlandm'},...
             {'Howard' },{'Hyytiala'},{'Kruger'  },{'Loobos'  },{'Merbleue'},...
             {'Mopane' },{'Palang'  },{'Sylvania'},{'Tumba'   },{'UniMich' }];
Nsites = length(siteNames);

% type names
typeNames = [{'all'}];%,{'dry'},{'wet'}];

% target names
targNames = [{'Qe'},{'Qh'},{'NEE'}];
targs = [1,2,3];
Ntargs = length(targs);

% training fractions
fracs = 1e2*2.^(0:10);
Nfracs = length(fracs);
if iType == 3; fracs = round(fracs/1.5); end;

% number of bootstraps
Nboots = 10;

% load training/test data from all sites
fprintf('Loading data: Type = %s - Site =  ',typeNames{iType}); tic; % screen report
for iSite = 1:Nsites
 fprintf('%d, ',iSite);

 fname = strcat('../data/lagged_data/',typeNames{iType},'_',siteNames{iSite},'.mat');
 load(fname);

 Ndat = size(data.date,1);
 if iType==3; Nuse = Ndat; 
 else Nuse = floor(Dmax/Nsites); end
% Nuse = min(Ndat,floor(Dmax/Nsites));
% assert(Ndat>Dmax/Nsites);

 Idx{iSite} = randperm(Ndat,Nuse);
 
 if iSite == 1 
  date = data.date(Idx{iSite},:); 
  forc = data.forc(Idx{iSite},:,:); 
  prog = data.flux(Idx{iSite},:);
 else 
  date = cat(1,date,data.date(Idx{iSite},:));
  forc = cat(1,forc,data.forc(Idx{iSite},:,:));
  prog = cat(1,prog,data.flux(Idx{iSite},:));
 end

end
t = toc; fprintf('finished: %f seconds. \n',t);

% make sure forcings are not missing
assert(isempty(find(isnan(forc(:)))));
assert(isempty(find(forc(:)<-9900)));
assert(isempty(find(prog(:)<-9900)));

% number of target types
assert(Ntargs <= size(prog,2));

% loop through bootstrap
for iBoot = 5:Nboots
 for iFrac = 1:Nfracs

  % screen report
  fprintf('Type = %d - Targ = %d - Boot = %d - Frac = %d \n',iType,iTarg,iBoot,iFrac); repTime = tic;

  % extract inputs
  X = reshape(forc,[size(forc,1),size(forc,2)*size(forc,3)]);
  X = cat(2,date(:,2:3),X);

  % extract targets
  Y = prog(:,targs(iTarg));

  % remove missing data
  I = find(isnan(Y));
  Y(I) = []; X(I,:) = [];

  % check sizes
  assert(size(X,1)==size(Y,1));
  N = size(X,1);
%  [N,fracs(end)]
  assert(N>1.5*fracs(end));

  % select training points
  I = randperm(N,fracs(iFrac)); 

  % select test points
  J = 1:N; J(ismember(J,I)) = [];  
  assert(length(J) >= length(I)/3);

  % perfom experiment
  results = ann_experiment(X,Y,I,J,trnfctn,Nepochs);

  % add site indexes to output stream
  results.Isites = Idx;

  % save results
  fprintf('Saving results ... '); tic; % screen report
  fname = strcat('results/',typeNames{iType},'_',targNames{iTarg},'_',num2str(fracs(iFrac)),'_',num2str(iBoot),'_',num2str(Nepochs),'_',trnfctn,'.mat');
  save(fname,'results');
  t = toc; fprintf('finished: %f seconds. \n',t);

  % screen report
  t = toc(repTime); fprintf('\nAll Finished - time = %f ',t); 
  fprintf('\n\n\n----------------------------------------------------\n\n\n');
 
 end % iFrac
end % iBoot

%% **** End Experiment ****************************************

% report total time
t = toc(masterTime); fprintf('Total time = %f',t);

% report that this matlab job is finished
fname = strcat('reports/convergence_bootstrap_end_report_',num2str(iType),'_',num2str(iTarg),'.end');
fid = fopen(fname,'w'); fprintf(fid,'done'); fclose(fid);


