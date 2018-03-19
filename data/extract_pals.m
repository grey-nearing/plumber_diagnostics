clear all
close all
clc

%% *** Site Data *************************************

% maximum number of points per site
Tmax = 227904; 

% site names
siteNames = [{'Amplero'},{'Blodgett'},{'Bugac'   },{'ElSaler2'},{'ElSaler' },...
             {'Espirra'},{'FortPeck'},{'Harvard' },{'Hesse'   },{'Howlandm'},...
             {'Howard' },{'Hyytiala'},{'Kruger'  },{'Loobos'  },{'Merbleue'},...
             {'Mopane' },{'Palang'  },{'Sylvania'},{'Tumba'   },{'UniMich' }];

siteAbbrv = [{'IT-Amp'},{'US-Blo'},{'HU-Bug'},{'ES-ES2'},{'ES-ES1'},...
             {'PT-Esp'},{'US-FPe'},{'US-Ha1'},{'FR-Hes'},{'US-Ho1'},...
             {'AU-How'},{'FI-Hyy'},{'ZA-Kru'},{'NL-Loo'},{'CA-Mer'},...
             {'BW-Ma1'},{'ID-Pag'},{'US-Syv'},{'AU-Tum'},{'US-UMB'}];

assert(length(siteNames)==length(siteAbbrv));

% dimensions
Nsites  = length(siteNames);

%% *** Loop Through Sites to Extgact Data *************

for s = 1:Nsites

 % screen report
 fprintf('Site %d of %d ... ',s,Nsites); tic;

 % netcdf file name for met and obs variables at site
 mname = strcat('pals_data/netcdf/',siteNames{s},'Fluxnet.1.4_met.nc');
 oname = strcat('pals_data/netcdf/',siteNames{s},'Fluxnet.1.4_flux.nc');

 % extract lat/lon from met file
 mlat = ncread(mname,'latitude');
 mlon = ncread(mname,'longitude');
 olat = ncread(oname,'latitude');
 olon = ncread(oname,'longitude');
 assert(mlat == olat); lat(s) = mlat;
 assert(mlon == olon); lon(s) = mlon;

 % read forcing data
 UU = ncread(mname,'Wind');
 Ta = ncread(mname,'Tair');
 Lw = ncread(mname,'LWdown');
 Sw = ncread(mname,'SWdown');
 Rh = ncread(mname,'Qair');
 PP = ncread(mname,'Rainf');
 Qe = ncread(oname,'Qle');
 Qh = ncread(oname,'Qh');
 NEE= ncread(oname,'NEE');

 % init storage
 date = zeros(Tmax,3)./0;
 prog = zeros(Tmax,5)./0;
 forc = zeros(Tmax,5)./0;

 % read spreadsheet for the rest of obs
 edex = 0;  % index looper
 D = 77;    % number of columns in spreadsheet

 % loop through years
 for y = 1994:2006 % loop through all possible years

  Nt = 365*48;
  if rem(y,4) == 0; Nt = Nt + 48; end

  % augment index counters
  sdex = edex+1;
  edex = sdex+Nt-1;
  dex = sdex:edex;

  try
   % read the file
   fname = strcat('pals_data/spreadsheets/',siteAbbrv{s},'.',num2str(y),'.synth.hourly.coreplusquality.csv');
   fid = fopen(fname);
   header = textscan(fid,'%s',D,'delimiter',',');
   data   = textscan(fid,'%f','delimiter',',');
   fclose(fid);
  catch
   continue
  end % try

  % check for correct number of entries in file
  Ny = size(data{1},1)/D;
  assert(Ny==Nt);
  data = reshape(data{1},D,Nt)';

  % assign dates
  date(dex,1) = y;                   date(edex,1) = y+1;
  date(dex,2) = floor((1:Nt)./48)+1; date(edex,2) = 1;
  date(dex,3) = rem(1:Nt,48)/2;      date(edex,3) = 0;

  prog(dex,1) = data(:,8);  % Qe
  prog(dex,2) = data(:,9);  % Qh
  prog(dex,3) = data(:,5);  % NEE
  prog(dex,4) = data(:,16); % SM1
  prog(dex,5) = data(:,17); % SM2

  forc(dex,1) = data(:,18); % Uu
  forc(dex,2) = data(:,11); % Ta
  forc(dex,3) = data(:,21); % Rn
  forc(dex,4) = data(:,36); % Rh
  forc(dex,5) = data(:,15); % Pp

 end % year

 % check that lengths match netcdf file
 assert(length(Qe) ==length(find(~isnan(prog(:,1)))));
 assert(length(Qh) ==length(find(~isnan(prog(:,2)))));
 assert(length(NEE)==length(find(~isnan(prog(:,3)))));
 assert(length(UU) ==length(find(~isnan(forc(:,1)))));
 assert(length(Ta) ==length(find(~isnan(forc(:,2)))));
 assert(length(Rh) ==length(find(~isnan(forc(:,3)))));
 assert(length(Lw) ==length(find(~isnan(forc(:,4)))));
 assert(length(PP) ==length(find(~isnan(forc(:,5)))));

 % check if spreadsheet matches netcdf file
 I = find(~isnan(prog(:,1)));
 if ~(max(abs(squeeze(Qe)    - prog(I,1)))<1e-6);
  fprintf('\n   Qe(%f,%f) ...',max(squeeze(Qe) - prog(I,1)),min(squeeze(Qe) - prog(I,1)));
  prog(I,1) = Qe;
 end
 if ~(max(abs(squeeze(Qh)    - prog(I,2)))<1e-6);
  fprintf('\n   Qh(%f,%f) ...',max(squeeze(Qh) - prog(I,2)),min(squeeze(Qh) - prog(I,2)));
  prog(I,2) = Qh;
 end
 if ~(max(abs(squeeze(NEE)   - prog(I,3)))<1e-6);
  fprintf('\n   NEE(%f,%f) ...',max(squeeze(NEE) - prog(I,3)),min(squeeze(NEE) - prog(I,3)));
  prog(I,3) = NEE;
 end
 if ~(max(abs(squeeze(UU)    - forc(I,1)))<1e-3);
  fprintf('\n   Uu(%f,%f) ...',max(squeeze(UU) - forc(I,1)),min(squeeze(UU) - forc(I,1)));
  forc(I,1) = UU;
 end
 if ~(max(abs(squeeze(Ta)    - forc(I,2)-273.15))<1e-3);
  fprintf('\n   Ta(%f,%f) ...',max(squeeze(Ta) - forc(I,2)),min(squeeze(Ta) - forc(I,2)));
  forc(I,2) = Ta;
 end
 if ~(max(abs(squeeze(Lw+Sw) - forc(I,3)))<1e-6);
  fprintf('\n   Rn(%f,%f) ...',max(squeeze(Lw+Sw) - forc(I,3)),min(squeeze(Lw+Sw) - forc(I,3)));
  forc(I,3) = Lw+Sw;
 end
 if ~(max(abs(squeeze(Rh)    - forc(I,4)))<1e-6);
  fprintf('\n   Rh(%f,%f) ...',max(squeeze(Rh) - forc(I,4)),min(squeeze(Rh) - forc(I,4)));
  forc(I,4) = Rh;
 end
 if ~(max(abs(squeeze(PP)    - forc(I,5)./1800))<1e-6);
  fprintf('\n   PP(%f,%f) ...',max(squeeze(PP) - forc(I,5)./1800),min(squeeze(PP) - forc(I,5)./1800));
  forc(I,5) = PP*1800;
 end

% % check file length
% assert(isempty(find(date(:,1)<0))); 
% assert(isempty(find(isnan(date(dex,1))))); 

 % write total data file
 outdata = [date,forc,prog];
 fname = strcat('pals_data/extracted/',siteNames{s},'.txt');
 save(fname,'outdata','-ascii') 

 % screen report
 t=toc; fprintf('finished: time = %f \n',t);

end % sites
 
%% *** Print Lat/Lon File *************************************

% print lat/lons of all sites in a single file
fid = fopen('LatLon.txt','w');
for s = 1:Nsites
 fprintf(fid,'%f %f \n',lat(s),lon(s));
end
fclose(fid);





