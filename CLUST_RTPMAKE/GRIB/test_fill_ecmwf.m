% fill_ecmwf.m
% 
% L. Strow, 22 Oct 2014
%
% Modify to include era?

function [prof, head, pattr] = fill_ecmwf(prof, head, pattr);

% Check args in and out to see if they conform to the new API and
% aren't split between old and new styles
if nargin ~= nargout
    error(['>>> ERROR: mismatch between fill_ecmwf inputs and ' ...
           'outputs.\n\tUse either [p,h]=fill_ecmwf(p,h) or ' ...
           '[p,h,pa]=fill_ecmwf(p,h,pa) (preferred)\n\tTerminating'], '\n');
end

addpath /asl/matlib/aslutil
addpath /asl/packages/time

% Location of grib files
fhdr = '/asl/data/ecmwf_nc/';

ename = '';  % This should be placed outside a rtp file loop

mtime = tai2dnum(prof.rtime);

% Get a cell array of ecmwf grib files for each time
% I think this will be BROKEN if using datetime above!!
enames = get_ecmwf_enames(mtime);

% Find the unique grib files and indices that go with them
[u_enames, ia, ic] = unique(enames);
n = length(u_enames);

keyboard

% Loop over unique grib file names
for i = 1:n
% Build file name from parts
   fne = ['UAD' u_enames{i} '001'];
   e_mth_year = datestr(mtime(ia(i)),'yyyymm');
   fn = fullfile(fhdr,e_mth_year(1:4),e_mth_year(5:6),fne);
% Actually read grib1, grib2 .nc files
   fn_s = [fn '-1.nc']
   fn_h = [fn '-2.nc'];
% Do the netcdf files exist?
   if exist(fn_s) == 0 | exist(fn_h) == 0 
      disp(['Netcdf grib files missing for root: ' fn])
      break % Go to next partition
   end
% If the filename has changed, re-load F   
%keyboard
   if ~strcmp(ename,fn) 
      clear F  % Probably not needed
      disp('New file')
      ename = fn;
   end   
% Fill rtp fields
   k = find( ic == i );  % indices of first partition (of n total)
% Assume rtp lat/lon are +-180??  Need to be 0-360 for grib interpolation
   rlat = prof.rlat(k);
   rlon = prof.rlon(k);
   rlon(rlon<0) = rlon(rlon<0) + 360;


% F.tcwv.ig  % Total column water?  Use this instead of ours (Sergio?)?
% F.msl.ig   % Not in rtp for now
end

