function [prof, head, pattr, iOLR] = fill_era5_monthly(prof, head, pattr, iOLR)

% fill_era5_monthly.m, copied from fill_era
% 
% % S.Machado - tweaks through the years
%

%{
sfc = read_netcdf_lls('/asl/models/era5_avg/2002/2002-08_sfc.nc');
lev = read_netcdf_lls('/asl/models/era5_avg/2002/2002-08_lev.nc');
era5plevs =  lev.level;
save /home/sergio/MATLABCODE/RTPMAKE/CLUST_RTPMAKE/CLUSTMAKE_ERA5/era5plevs.mat era5plevs
NOTE era5plevs are from 1 to 1000 mb ie automatically already sorted from minimum to maximum
%}

if nargin < 3
  iOLR = -1;
end

% Check args in and out to see if they conform to the new API and
% aren't split between old and new styles
if nargin ~= nargout
    error(['>>> ERROR: mismatch between fill_era5_monthly inputs and ' ...
           'outputs.\n\tUse either [p,h]=fill_era5_monthly(p,h) or ' ...
           '[p,h,pa]=fill_era5_monthly(p,h,pa) (preferred)\n\tTerminating or also include iOLR'], '\n')
end

addpath /asl/matlib/aslutil
addpath /asl/packages/time

% Location of grib files
fhdr = '/asl/models/era5_avg/';         %% from July 2021
fhdr = '/asl/models/era5_monthly/';     %% from Jan 2024

ename = '';  % This should be placed outside a rtp file loop
mtime = tai2dnum(prof.rtime);

% Get a cell array of ecmwf grib files for each time
% I think this will be BROKEN if using datetime above!!
% enames = get_ecmwf_enames(mtime,profin.rtime);
%enames = get_era5_monthly_enames(mtime);
enames = get_era5_monthly_enames(mtime,1,fhdr);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get a cell array of era grib files for each time
% Round to get 8 forecast hours per day
rmtime = round(mtime*8)/8;
timestr = datestr(rmtime,'yyyymmddhh');
ystr = timestr(:,1:4);
mstr = timestr(:,5:6);
dstr = timestr(:,7:8);
hstr = timestr(:,9:10);
yearindex = str2num(ystr);
monthindex = str2num(mstr);
dayindex = str2num(dstr);
hourindex = str2num(hstr);

%% hmmm ... oct 30, 2022
[boo,~] = size(hstr);
for ii = 1 : boo
  hourindex(ii) = str2num(hstr(ii,:));  
end
figure(1); plot(yearindex);
figure(2); plot(monthindex);
figure(3); plot(dayindex);
figure(4); plot(hourindex);
figure(5); hist(hourindex,25); 

%{ 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% enames = [ystr mstr dstr];
% enames = cellstr(enames);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%} 

xnames = enames.fn;
[u_enames, ~, ic] = unique(xnames);

n = length(u_enames); % Generally 2 names for 1 day's worth of data

% ic = ones(size(mtime));
% n = 1;

fprintf(1,'need to read in %3i ERA5 monthly files \n',n)

load /home/sergio/MATLABCODE/RTPMAKE/CLUST_RTPMAKE/CLUSTMAKE_ERA5/era5plevs.mat

for i=1:n
   %fn = fullfile(fhdr,u_enames{i}(1:4),u_enames{i}(5:6),u_enames{i});
   %fn_lev = [fn '_lev.nc'];
   %fn_sfc = [fn '_sfc.nc'];

   fn     = u_enames{i};
   fn_lev = [fn '_lev.nc'];
   fn_sfc = [fn '_sfc.nc'];
   fn_OLR = [fn '_rad.nc'];
   fn_2m  = [fn '_2meter.nc'];
     fn_2m = fn_sfc;

%   kapoo = findstr(fn_OLR,'_avg');
%   kamoo = fn_OLR(kapoo+4:kapoo+9);
%   fn_OLR = strrep(fn_OLR,kamoo,'/INCOMING/')
%   %% fn_OLR = '/asl/models/era5_avg/INCOMING/2002-09_rad.nc';
%   %% keyboard_nowindow

   fprintf(1,'era5 monthly file name fn = %s \n',[fn '_*.nc where * = sfc or lev or 2meter'])    
   fprintf(1,'  eg fn_sfc = %s fn_lev = %s fn_OLR = %s  \n',fn_sfc,fn_lev,fn_OLR)

   % Does the netcdf files exist?
   if exist(fn_sfc,'file') == 0 || exist(fn_lev,'file') == 0 
      disp(['Netcdf grib files missing for root: ' fn])
      break % Go to next partition
   end
   % If the filename has changed, re-load F   
   if ~strcmp(ename,fn) 
      clear F  % Probably not needed
      if ~exist(fn_2m,'file') & ~exist(fn_OLR,'file')
        disp('did not find 2m file OR olr ....')
        %disp('New ERA5 monthly file'); %for debugging
        F(1) = grib_interpolate_era(fn_sfc,fn_lev,1);
        F(2) = grib_interpolate_era(fn_sfc,fn_lev,2);
        F(3) = grib_interpolate_era(fn_sfc,fn_lev,3);
        F(4) = grib_interpolate_era(fn_sfc,fn_lev,4);
        F(5) = grib_interpolate_era(fn_sfc,fn_lev,5);
        F(6) = grib_interpolate_era(fn_sfc,fn_lev,6);
        F(7) = grib_interpolate_era(fn_sfc,fn_lev,7);
        F(8) = grib_interpolate_era(fn_sfc,fn_lev,8);
        ename = fn;
      elseif exist(fn_2m,'file') & ~exist(fn_OLR,'file')
        disp('yay, did find 2m file but not OLR ....')
        F(1) = grib_interpolate_era5(fn_sfc,fn_lev,fn_2m,fn_OLR,1,+1);
        F(2) = grib_interpolate_era5(fn_sfc,fn_lev,fn_2m,fn_OLR,2,+1);
        F(3) = grib_interpolate_era5(fn_sfc,fn_lev,fn_2m,fn_OLR,3,+1);
        F(4) = grib_interpolate_era5(fn_sfc,fn_lev,fn_2m,fn_OLR,4,+1);
        F(5) = grib_interpolate_era5(fn_sfc,fn_lev,fn_2m,fn_OLR,5,+1);
        F(6) = grib_interpolate_era5(fn_sfc,fn_lev,fn_2m,fn_OLR,6,+1);
        F(7) = grib_interpolate_era5(fn_sfc,fn_lev,fn_2m,fn_OLR,7,+1);
        F(8) = grib_interpolate_era5(fn_sfc,fn_lev,fn_2m,fn_OLR,8,+1);
        ename = fn;
      elseif ~exist(fn_2m,'file') & exist(fn_OLR,'file')
        disp('yay, did NOT find 2m file but found OLR ....')
        F(1) = grib_interpolate_era5(fn_sfc,fn_lev,fn_2m,fn_OLR,1,-1);
        F(2) = grib_interpolate_era5(fn_sfc,fn_lev,fn_2m,fn_OLR,2,-1);
        F(3) = grib_interpolate_era5(fn_sfc,fn_lev,fn_2m,fn_OLR,3,-1);
        F(4) = grib_interpolate_era5(fn_sfc,fn_lev,fn_2m,fn_OLR,4,-1);
        F(5) = grib_interpolate_era5(fn_sfc,fn_lev,fn_2m,fn_OLR,5,-1);
        F(6) = grib_interpolate_era5(fn_sfc,fn_lev,fn_2m,fn_OLR,6,-1);
        F(7) = grib_interpolate_era5(fn_sfc,fn_lev,fn_2m,fn_OLR,7,-1);
        F(8) = grib_interpolate_era5(fn_sfc,fn_lev,fn_2m,fn_OLR,8,-1);
        ename = fn;
      elseif exist(fn_2m,'file') & exist(fn_OLR,'file')
        disp('yay, found 2m file and found OLR ....')
        F(1) = grib_interpolate_era5(fn_sfc,fn_lev,fn_2m,fn_OLR,1,0);
        F(2) = grib_interpolate_era5(fn_sfc,fn_lev,fn_2m,fn_OLR,2,0);
        F(3) = grib_interpolate_era5(fn_sfc,fn_lev,fn_2m,fn_OLR,3,0);
        F(4) = grib_interpolate_era5(fn_sfc,fn_lev,fn_2m,fn_OLR,4,0);
        F(5) = grib_interpolate_era5(fn_sfc,fn_lev,fn_2m,fn_OLR,5,0);
        F(6) = grib_interpolate_era5(fn_sfc,fn_lev,fn_2m,fn_OLR,6,0);
        F(7) = grib_interpolate_era5(fn_sfc,fn_lev,fn_2m,fn_OLR,7,0);
        F(8) = grib_interpolate_era5(fn_sfc,fn_lev,fn_2m,fn_OLR,8,0);
        ename = fn;
      end
   end   

   % Fill rtp fields
   m = find( ic == i );  % indices of first era file
   %   fhi = 0;   % this was new on Jul 20, 2015!
   u_hour = unique(hourindex);
   nn = length(u_hour);
   fprintf(1,'  ERA file %3i of %3i --- %6i of %6i entries contained inside %3i unique hours \n',i,n,length(m),length(prof.rtime),nn);

   % Only loop over hours needed

   for jj = 1:nn
      % index for this hour (1:8);  u_hour = [0 3 6 9 12 15 18 21]
      fhi = floor(u_hour(jj)/3) + 1;
%% fhi = floor([0:23]/3) + 1
%% fhi =      1     1     1     2     2     2     3     3     3     4     4     4     5     5     5     6     6     6     7     7     7     8     8     8
      l = find( hourindex == u_hour(jj));
      k = intersect(l,m);
      % sfhi(k,:) = fhi;   % Debug, showed that fhi changes properly
      if k > 0         
         % Assume rtp lat/lon are +-180??  Need to be 0-360 for grib interpolation
         rlat = prof.rlat(k);
         rlon = prof.rlon(k);
         rlon(rlon<0) = rlon(rlon<0) + 360;

         try
           %% 2m air and dewpoint temperatures
           prof.d2m(k)   = F(fhi).d2m.ig(rlat,rlon);
           prof.t2m(k)   = F(fhi).t2m.ig(rlat,rlon);
         end

         if iOLR > 0
           try
             prof.olr(k)     = F(fhi).olr.ig(rlat,rlon);
             prof.olr_clr(k) = F(fhi).olr_clr.ig(rlat,rlon);
             prof.ilr(k)     = F(fhi).ilr.ig(rlat,rlon);
             prof.ilr_clr(k) = F(fhi).ilr_clr.ig(rlat,rlon);
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
         % NOTE era5plevs are from 1 to 1000 mb ie automatically already sorted from minimum to maximum   
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
         %xtemp = p60_ecmwf(prof.spres(k));  % all 137 pressure levels         
         %prof.plevs(:,k) = xtemp(b,:);  % subset to ones in grib file
         %prof.nlevs(k) = length(F(fhi).levid);

         % NOTE era5plevs are from 1 to 1000 mb ie automatically already sorted from minimum to maximum   
         xtemp = double(era5plevs);          
         xtemp = double(sort(F(fhi).levid));  %% THIS IS EQUIVALENT
         prof.plevs(1:37,k) = xtemp * ones(1,length(k));  % subset to ones in grib file
         prof.nlevs(k) = 37;

      end  % k loop  LLS
   end
end
prof.nlevs = int32(prof.nlevs);

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
%    say(['Replaced ' int2str(nbad) ' negative/zero H2O mixing ratios'])
  end
end
%
if isfield(prof,'gas_3')
  ibad = find(prof.gas_3 <= 0);
  nbad = length(ibad);
  if (nbad > 0)
    prof.gas_3(ibad) = min_O3_gg;
%    say(['Replaced ' int2str(nbad) ' negative/zero O3 mixing ratios'])
  end
end
%  fix any cloud frac
if isfield(prof,'tcc')
  ibad = find(prof.tcc > 1);
  nbad = length(ibad);
  if (nbad > 0)
    prof.tcc(ibad) = 1;
%    say(['Replaced ' int2str(nbad) ' TCC > 1 fields'])
  end
end

switch nargin
  case 2
    fprintf(2, ['>>> WARNING: fill_era5_monthly now sets model attribute in ' ...
                'pattr.\n\tUpdate calls to fill_era5_monthly to include pattr. ' ...
                'i.e. [p,h,pa] = fill_era5_monthly(p,h,pa)\n'])
  case 3
    % set an attribute string to let the rtp know what we have done
    pattr = set_attr(pattr,'model','era5_monthly');
end

