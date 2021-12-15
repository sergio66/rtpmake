% fill_era5_daily.m
% 
% S.Machado - tweaks through the years
%

%{
lev = read_netcdf_lls('/asl/models/era5/2014/01/20140102_lev_test.nc');
sfc = read_netcdf_lls('/asl/models/era5/2014/01/20140102_sfc.nc');
%}

function [prof, head, pattr] = fill_era5_daily(prof, head, pattr);

addpath /asl/matlib/aslutil
addpath /asl/packages/time

%{
if nargin ~= nargout
    error(['>>> ERROR: mismatch between fill_era5_daily inputs and ' ...
           'outputs.\n\tUse either [p,h]=fill_era5_daily(p,h) or ' ...
           '[p,h,pa]=fill_era5_daily(p,h,pa) (preferred)\n\tTerminating'], '\n');
end
%}

addpath /asl/matlib/aslutil
addpath /asl/packages/time

% Location of grib files
fhdr = '/asl/models/era5/';   %% from Oct 2019 - 

ename = '';  % This should be placed outside a rtp file loop
mtime = tai2dnum(prof.rtime);

% Get a cell array of ecmwf grib files for each time
% I think this will be BROKEN if using datetime above!!
% enames = get_ecmwf_enames(mtime,prof.rtime);
enames = get_era5_daily_enames(mtime);

% Get a cell array of era grib files for each time
% Round to get 4 forecast hours per day
rmtime = round(mtime*8)/8;
timestr = datestr(rmtime,'yyyymmddhh');
ystr = timestr(:,1:4);
mstr = timestr(:,5:6);
dstr = timestr(:,7:8);
hstr = timestr(:,9:10);
yearindex = str2num(ystr);
dayindex = str2num(dstr);
hourindex = str2num(hstr);
  figure(1); hist(hourindex,25);

%enames = [ystr mstr dstr];
%enames = cellstr(enames);
%[u_enames, ~, ic] = unique(enames);
%n = length(u_enames); % Generally 2 names for 1 day's worth of data
ic = ones(size(mtime));
n = 1;

for i=1:n
   %fn = fullfile(fhdr,u_enames{i}(1:4),u_enames{i}(5:6),u_enames{i});
   %fn_lev = [fn '_lev.nc'];
   %fn_sfc = [fn '_sfc.nc'];
   fn     = enames.fn;
   fn_lev = enames.lev;
   fn_sfc = enames.sfc;
   fn_2m  = enames.twm;
   fprintf(1,'era5 daily file name fn = %s \n',[fn '_*.nc where * = sfc or lev_test or 2meter'])   

   % Does the netcdf files exist? dont worry too much about fn_2m
   if exist(fn_sfc,'file') == 0 || exist(fn_lev,'file') == 0 
      disp(['Netcdf grib files missing for root: ' fn])
      break % Go to next partition
   end
   % If the filename has changed, re-load F   
   if ~strcmp(ename,fn) 
      clear F  % Probably not needed
      %disp('New ERA5 monthly file'); %for debugging
      if ~exist(fn_2m,'file')
        disp('did not find 2m file ....')
        F(1) = grib_interpolate_era(fn_sfc,fn_lev,1);
        F(2) = grib_interpolate_era(fn_sfc,fn_lev,2);
        F(3) = grib_interpolate_era(fn_sfc,fn_lev,3);
        F(4) = grib_interpolate_era(fn_sfc,fn_lev,4);
        F(5) = grib_interpolate_era(fn_sfc,fn_lev,5);
        F(6) = grib_interpolate_era(fn_sfc,fn_lev,6);
        F(7) = grib_interpolate_era(fn_sfc,fn_lev,7);
        F(8) = grib_interpolate_era(fn_sfc,fn_lev,8);
        ename = fn;
      else
        disp('yay, did find 2m file ....')
        F(1) = grib_interpolate_era5(fn_sfc,fn_lev,fn_2m,1);
        F(2) = grib_interpolate_era5(fn_sfc,fn_lev,fn_2m,2);
        F(3) = grib_interpolate_era5(fn_sfc,fn_lev,fn_2m,3);
        F(4) = grib_interpolate_era5(fn_sfc,fn_lev,fn_2m,4);
        F(5) = grib_interpolate_era5(fn_sfc,fn_lev,fn_2m,5);
        F(6) = grib_interpolate_era5(fn_sfc,fn_lev,fn_2m,6);
        F(7) = grib_interpolate_era5(fn_sfc,fn_lev,fn_2m,7);
        F(8) = grib_interpolate_era5(fn_sfc,fn_lev,fn_2m,8);
        ename = fn;
      end
     end   

   % Fill rtp fields
   m = find( ic == i );  % indices of first era file
   %   fhi = 0;   % this was new on Jul 20, 2015!
   u_hour = unique(hourindex);
   nn = length(u_hour);
   % Only loop over hours needed
   for jj = 1:nn
      % index for this hour (1:8);  u_hour = [0 3 6 9 12 15 18 21]
      fhi = (u_hour(jj)/3) + 1;
      l = find( hourindex == u_hour(jj));
      k = intersect(l,m);
      % sfhi(k,:) = fhi;   % Debug, showed that fhi changes properly
      if k > 0         
         % Assume rtp lat/lon are +-180??  Need to be 0-360 for grib interpolation
         rlat = prof.rlat(k);
         rlon = prof.rlon(k);
         rlon(rlon<0) = rlon(rlon<0) + 360;

         if exist(fn_2m)
           try
             %% 2m air and dewpoint temperatures
             prof.d2m(k)   = F(fhi).d2m.ig(rlat,rlon);
             prof.t2m(k)   = F(fhi).t2m.ig(rlat,rlon);
           end
         end

         prof.spres(k)   = F(fhi).sp.ig(rlat,rlon);
         prof.stemp(k)   = F(fhi).skt.ig(rlat,rlon);
         wind_v          = F(fhi).v10.ig(rlat,rlon);
         wind_u          = F(fhi).u10.ig(rlat,rlon);
         prof.wspeed(k)  = sqrt(wind_u.^2 + wind_v.^2);
         prof.wsource(k) = mod(atan2(single(wind_u), single(wind_v)) * 180/pi,360);
         prof.tcc(k)     = F(fhi).tcc.ig(rlat,rlon);
         ci_udef = 1;
         prof.udef(ci_udef,k) = F(fhi).ci.ig(rlat,rlon);

         % Estimate model grid centers used
         gdlat = abs(nanmean(diff(F(fhi).h_latitude)));  % lat spacing
         gdlon = abs(nanmean(diff(F(fhi).h_longitude))); % lon spacing
         prof.plat(k) = floor(rlat/gdlat)*gdlat + gdlat/2;
         prof.plon(k) = floor(rlon/gdlon)*gdlon + gdlon/2;

         % F(fhi).tcwv.ig  % Total column water?  Use this instead of ours?
         % F(fhi).msl.ig   % Not in rtp for now
         % Hybrid parameters
         % levid = 1 is top of atmosphere
         % b are the sortedd level IDs   
         [b,j]=sort(F(fhi).levid);
         for l=1:length(F(fhi).levid)
            prof.ptemp(l,k) = F(fhi).t(j(l)).ig(rlat,rlon);
            prof.gas_1(l,k) = F(fhi).q(j(l)).ig(rlat,rlon);
            prof.gas_3(l,k) = F(fhi).o3(j(l)).ig(rlat,rlon);
            prof.cc(l,k)    = F(fhi).cc(j(l)).ig(rlat,rlon);
            prof.clwc(l,k)  = F(fhi).clwc(j(l)).ig(rlat,rlon);
            prof.ciwc(l,k)  = F(fhi).ciwc(j(l)).ig(rlat,rlon);
         end
% Only want pressure levels in grib file, in order
         xtemp = p137_ecmwf(prof.spres(k));  % all 137 pressure levels         
         prof.plevs(:,k) = xtemp(b,:);  % subset to ones in grib file
         prof.nlevs(k) = length(F(fhi).levid);

%         xtemp = double(era5plevs);
%         prof.plevs(1:37,k) = xtemp * ones(1,length(k));  % subset to ones in grib file
%         prof.nlevs(k) = 37;

      end  % k loop  LLS
   end
end

% Header info
head.ptype = 0;
head.ngas = 2;
head.glist = [1; 3];
head.gunit = [21; 21];
head.pmin = min( prof.plevs(1,:) );
head.pmax = max( prof.plevs(end,:) );

% Setting attributes needs work...
% pattr = set_attr(pattr,'profiles','ECMWF','profiles');

% I think this is needed to avoid negatives in SARTA?
min_H2O_gg = 3.1E-7;  % 0.5 pppm
min_O3_gg = 1.6E-8;   % 0.01 ppm
% Find/replace bad mixing ratios
if isfield(prof,'gas_1')
  ibad = find(prof.gas_1 <= 0);
  nbad = length(ibad);
  if (nbad > 0)
    prof.gas_1(ibad) = min_H2O_gg;
   %say(['Replaced ' int2str(nbad) ' negative/zero H2O mixing ratios'])
  end
end

if isfield(prof,'gas_3')
ibad = find(prof.gas_3 <= 0);
nbad = length(ibad);
  if (nbad > 0)
    prof.gas_3(ibad) = min_O3_gg;
    %say(['Replaced ' int2str(nbad) ' negative/zero O3 mixing ratios'])
  end
end

%  fix any cloud frac
if isfield(prof,'tcc')
  ibad = find(prof.tcc > 1);
  nbad = length(ibad);
  if (nbad > 0)
    prof.tcc(ibad) = 1;
    % say(['Replaced ' int2str(nbad) ' TCC > 1 fields'])
  end
end
if isfield(prof,'tcc')
  ibad = find(prof.tcc < 0);
  nbad = length(ibad);
  if (nbad > 0)
    prof.tcc(ibad) = 0;
    % say(['Replaced ' int2str(nbad) ' TCC > 1 fields'])
  end
end
if isfield(prof,'cc')
  ibad = find(prof.cc > 1);
  nbad = length(ibad);
  if (nbad > 0)
    prof.cc(ibad) = 1;
    %say(['Replaced ' int2str(nbad) ' CC > 1 fields'])
  end
end
if isfield(prof,'cc')
  ibad = find(prof.cc < 0);
  nbad = length(ibad);
  if (nbad > 0)
    prof.cc(ibad) = 0;
     say(['Replaced ' int2str(nbad) ' CC > 1 fields'])
  end
end
if isfield(prof,'ciwc')
  ibad = find(prof.ciwc < 0);
  nbad = length(ibad);
  if (nbad > 0)
    prof.ciwc(ibad) = 0;
    % say(['Replaced ' int2str(nbad) ' CIWC > 1 fields'])
  end
end
if isfield(prof,'clwcc')
  ibad = find(prof.clwcc < 0);
  nbad = length(ibad);
  if (nbad > 0)
    prof.clwc(ibad) = 0;
    % say(['Replaced ' int2str(nbad) ' CLWC > 1 fields'])
  end
end

switch nargin
  case 2
    fprintf(2, ['>>> WARNING: fill_ecmwf now sets model attribute in ' ...
                'pattr.\n\tUpdate calls to fill_ecmwf to include pattr. ' ...
                'i.e. [p,h,pa] = fill_ecmwf(p,h,pa)\n'])
  case 3
    % set an attribute string to let the rtp know what we have done
    pattr = set_attr(pattr,'model','ecmwf');
  case 4
    if iDoneOne < 0
      disp('OOPS : could not read in a single ECMWF match file!')
    end
end

