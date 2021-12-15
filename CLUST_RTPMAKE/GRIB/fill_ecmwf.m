% fill_ecmwf.m
% 
% L. Strow, 22 Oct 2014
% S.Machado - tweaks through the years
%
% Modify to include era?

function [prof, head, pattr, iDoneOne, iaDone] = fill_ecmwf(prof, head, pattr);

profin = prof;
headin = head;

%% carry this index around
headin.nchan = headin.nchan+1;
headin.ichan = [headin.ichan ; max(headin.ichan)+1];
if isfield(headin,'vchan')
  headin.vchan = [headin.vchan; max(headin.vchan)+1];
end
profin.robs1(headin.nchan,:) = 1:length(profin.rlat);
if isfield(profin,'calflag')
  profin.calflag(headin.nchan,:) = 1:length(profin.rlat);
end

clear prof
prof = struct;

%% 2017/03/19 removed iDoneOne as an input argument since it is always set to 1 a few lines below
%% if nargin == 3
%%  iDoneOne = -1;
%% end

% Check args in and out to see if they conform to the new API and
% aren't split between old and new styles

%{
if nargin ~= nargout
    error(['>>> ERROR: mismatch between fill_ecmwf inputs and ' ...
           'outputs.\n\tUse either [p,h]=fill_ecmwf(p,h) or ' ...
           '[p,h,pa]=fill_ecmwf(p,h,pa) (preferred)\n\tTerminating'], '\n');
end
%}

addpath /asl/matlib/aslutil
addpath /asl/packages/time

iaDone = zeros(size(profin.rlat));

% Location of grib files
fhdr = '/asl/data/ecmwf_nc/';  %% till Sept 2019
fhdr = '/asl/data/ecmwf/';     %% from Oct 2019 - 
fhdr = '/asl/models/ecmwf/';   %% from Oct 2019 - 

ename = '';  % This should be placed outside a rtp file loop

mtime = tai2dnum(profin.rtime);

% Get a cell array of ecmwf grib files for each time
% I think this will be BROKEN if using datetime above!!
% enames = get_ecmwf_enames(mtime,profin.rtime);
enames = get_ecmwf_enames(mtime);

% Find the unique grib files and indices that go with them
[u_enames, ia, ic] = unique(enames);
%u_enames
n = length(u_enames);

skipped = 0;

% Loop over unique grib file names
iDoneOne = -1;
for i = 1:n
  % Build file name from parts
  fne = ['UAD' u_enames{i} '001'];
  e_mth_year = datestr(mtime(ia(i)),'yyyymm');
  fn = fullfile(fhdr,e_mth_year(1:4),e_mth_year(5:6),fne);
  fprintf(1,'ecmwf file name fn = %s \n',fn);

  %% this ptime can mess up eg when looking at fovs on 12/31, about 23-24 hrs
  %% (since then it looks for files from next year, which I could have symbolicly linked in
  %% and pretended they are from 12/01 of this year)
  %%
  %% or never have symbolic links to files of next year, so the "close to next year" FOVS will not be matched
  %% to ECMWF files ==> "CONTINUE" statement from line 75 below takes into effect, and you
  %% won't have this problem ==> ptime will be fine

  addpath /home/sergio/MATLABCODE/TIME
  %% /asl/data/ecmwf_nc/2013/07/UAD07021200070221001
  %% /asl/data/ecmwf_nc/YYYY/MM/UADMMDD1200MMDD21001  
  slash = findstr(fn,'/');
  zahyear  = str2num(fn(slash(4)+1:slash(4)+4));
  zahmonth = str2num(fn(slash(5)+1:slash(5)+2));
  zahday   = str2num(fn(slash(6)+6:slash(6)+7));
  zahhr    = str2num(fn(slash(6)+8:slash(6)+9));
  zahhr    = str2num(fn(slash(6)+16:slash(6)+17));  
  ptime    = utc2taiSergio(zahyear,zahmonth,zahday,zahhr);

  % Actually read grib1, grib2 .nc files
  fn_s = [fn '-1.nc'];
  fn_h = [fn '-2.nc'];
  % Do the netcdf files exist?
  if exist(fn_s) == 0 | exist(fn_h) == 0 
    disp(['Netcdf grib files missing for root: ' fn])
    k = find( ic == i );  % indices of first partition (of n total)
    skipped = skipped + length(k);
    fprintf(1,'skipping to next time partition, will have to ignore %4i FOVS \n',length(k))
    continue % Go to next partition  ORIGNALLY WAS BREAK SO EXITED THE LOOP COMPLETELY
   end
   
  % If the filename has changed, re-load F   
  %keyboard
   if ~strcmp(ename,fn) 
      clear F  % Probably not needed
      disp('New file')
      F = grib_interpolate(fn_s,fn_h);
      ename = fn;
   end   
  % Fill rtp fields
   k = find( ic == i );  % indices of first partition (of n total)
   if ~isfield(headin,'ngas')
     headin.ngas = 0;
   end
   [head,prof] = subset_rtp(headin,profin,[],[],k);
   iaDone(k) = 1;

  % Assume rtp lat/lon are +-180??  Need to be 0-360 for grib interpolation
  %   rlat = prof.rlat(k);
  %   rlon = prof.rlon(k);
   rlat = prof.rlat;
   rlon = prof.rlon;
   rlon(rlon<0) = rlon(rlon<0) + 360;

   prof.ptime   = ones(1,length(k)) * ptime;
   prof.sst     = F.sst.ig(rlat,rlon);
   prof.spres   = F.sp.ig(rlat,rlon);
   prof.stemp   = F.skt.ig(rlat,rlon);
   wind_v          = F.v10.ig(rlat,rlon);
   wind_u          = F.u10.ig(rlat,rlon);
   prof.wspeed  = sqrt(wind_u.^2 + wind_v.^2);   
   prof.wsource = mod(atan2(single(wind_u), single(wind_v)) * 180/pi,360);
   prof.tcc   = F.tcc.ig(rlat,rlon);

   %% comment this out Oct 2020
   %ci_udef = 1;
   %prof.udef(ci_udef,:) = F.ci.ig(rlat,rlon);

   % Estimate model grid centers used
   gdlat = abs(nanmean(diff(F.h_latitude)));  % lat spacing
   gdlon = abs(nanmean(diff(F.h_longitude))); % lon spacing
   prof.plat = floor(rlat/gdlat)*gdlat + gdlat/2;
   prof.plon = floor(rlon/gdlon)*gdlon + gdlon/2;

  % F.tcwv.ig  % Total column water?  Use this instead of ours (Sergio?)?
  % F.msl.ig   % Not in rtp for now

  % Hybrid parameters
  % levid = 1 is top of atmosphere
  % b are the sortedd level IDs   
  %   prof.nlevs = ones(1,length(k))*length(F.levid);
   [b,j]=sort(F.levid);
   for l=1:length(F.levid)
      prof.ptemp(l,:) = F.t(j(l)).ig(rlat,rlon);
      prof.gas_1(l,:) = F.q(j(l)).ig(rlat,rlon);
      prof.gas_3(l,:) = F.o3(j(l)).ig(rlat,rlon);
      prof.cc(l,:)    = F.cc(j(l)).ig(rlat,rlon);
      prof.clwc(l,:)  = F.clwc(j(l)).ig(rlat,rlon);
      prof.ciwc(l,:)  = F.ciwc(j(l)).ig(rlat,rlon);
   end
  % Only want pressure levels in grib file, in order
  % Is this a 91 or 137 level forecast?
  % Note: On June 25, 2013 ECMWF moved to 137 levels, and they selected
  % 91 of these to send to us!!  Need to map them correctly!
  % Need to check, how many levels in Sept. 1, 2092 ECMWF?
   max_lev = max(F.levid);
   if max_lev > 91
      xtemp = p137_ecmwf(prof.spres);
   elseif max_lev == 91
      xtemp = p91_ecmwf(prof.spres);
   elseif max_lev < 91
      xtemp = p60_ecmwf(prof.spres);
   end
   prof.plevs(:,:) = xtemp(b,:);  % subset to ones in grib file
   prof.nlevs(1:length(k)) = length(F.levid);

   if iDoneOne < 0
     iDoneOne = +1;
     prof_all = prof;
     head_all = headin;
     if mod(head_all.pfields,2) == 0
       %% pfields = 2 or 4 so obs and.or calc, no prof
       head_all.pfields = head_all.pfields + 1;  %% add on profile
     end
     head_all.ptype = 0;                       %% levels     
     head_all.ngas = 2;
     head_all.glist = [1; 3];
     head_all.gunit = [21; 21];
     blah = [i n length(prof.stemp) length(profin.rlat)];
     fprintf(1,'loop %3i of %3i -- adding on %4i profiles  .. orig input total was %4i \n',blah);     
   else
     blah = [i n length(prof.stemp) length(prof_all.stemp) length(prof.stemp)+length(prof_all.stemp) length(profin.rlat)];
     fprintf(1,'loop %3i of %3i -- adding on %4i profiles to existing %4i for total %4i .. orig input total was %4i \n',blah);
     [head_all,prof_all] = cat_rtp(head_all,prof_all,head_all,prof);
   end
   plot((prof_all.rtime-prof_all.ptime)/60/60); title('\deltaT (rtime-ptime)');ylabel('hours');pause(0.1)
   
end

if exist('prof_all')
  if iDoneOne < 0
    error('hmm iDoneOne < 0 but there is prof_all???')
  end
  prof = prof_all;
  prof.nlevs = int32(prof.nlevs);
  fprintf(1,'done with loops ... final profiles %4i .. orig FOVS were %4i \n',length(prof.stemp),length(profin.rlat));
else
  if iDoneOne > 0
    error('hmm iDoneOne > 0 but there is no prof_all???')
  end
  fprintf('oops ... could not do anything!!!!')
  prof = profin;
  skipped = length(prof.rlat);
end  
if skipped > 0
  fprintf(1,'    skipped %4i \n',skipped)
end  
disp(' ')

% Header info
if exist('prof_all')
  if mod(head.pfields,2) == 0
    %% pfields = 2 or 4 so obs and.or calc, no prof
    head.pfields = head.pfields + 1;  %% add on profile
  end
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
%      say(['Replaced ' int2str(nbad) ' CC > 1 fields'])
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
else
  head = headin;
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% now re-sort
robsx = prof.robs1(head.nchan,:);
[Y,I] = sort(robsx);
if sum(abs(Y-I)) ~= 0
  index = I;
  disp('oops indices got messed while reading different fn files .... resorting')
  fn = fieldnames(prof);
  for ii = 1 : length(fn)
    str = ['junkjunk = prof.' fn{ii} ';'];
    eval(str);
    [m,n] = size(junkjunk);
    if m == 1
      junkjunk = junkjunk(index);
    else
      junkjunk = junkjunk(:,index);
    end
    str = ['prof.' fn{ii} ' =  junkjunk;'];
    eval(str);
  end
else
  disp('hmm no need to re-sort at end')
end

%% finally get rid of last point
head.nchan = head.nchan-1;
head.ichan = head.ichan(1:head.nchan);
if isfield(head,'vchan')
  head.vchan = head.vchan(1:head.nchan);
end  
prof.robs1 = prof.robs1(1:head.nchan,:);
if isfield(profin,'calflag')
  prof.calflag = prof.calflag(1:head.nchan,:);
end  
