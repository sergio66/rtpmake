function F = grib_interpolate_era5(fn_s,fn_h,fn_2m,fn_olr,hindex,iWhich);
% 
% Inputs: fn_s, fn_h, fn_2m, fn_olr, hindex, iWhich 
%         Netcdf files containing grib1 and grib2 
%         data respectively.  This code assumes the grib1 data
%         is surface data, and the grib2 data is hybrid ecmwf data.
%         hindex is the hour index (4 total) in the era file

% Output: F, a structure with housekeeping data (lat, lon of grid) and
%         interpolants for interpolating to arbitrary lat/lon positions.
%
% Presently the interpolant is set to be linear, you can also 
% select nearest-neighbor, cubic, etc.  I may make this an option
% in the future.  Minor changes are needed to use this for ERA, etc.
%
% L. Strow, June 11, 2014

if nargin <= 3
  error('grib_interpolate_era5 needs 4 arguments, not 3')
end
if hindex < 1 | hindex > 8
  error('need hindex == 1 .. 8 for grib_interpolate_era5.m');
end

F.s_longitude = ncread(fn_s,'longitude');
  %% longitude is 0 to 360
F.s_latitude  = ncread(fn_s,'latitude');
  %% latitude is +90 to -90
F.s_time      = ncread(fn_s,'time',[hindex],[1]);
F.s_mtime     = datenum(1900,0,0,double(F.s_time),0,0);

[X,Y] = ndgrid(F.s_latitude,F.s_longitude);
iX = flipud(X); iY = flipud(Y);

F.sp.ig   = griddedInterpolant(iX,iY,flipud(single(ncread(fn_s,'sp',[1 1 hindex],[Inf Inf 1]))'/100),'linear');
F.skt.ig  = griddedInterpolant(iX,iY,flipud(single(ncread(fn_s,'skt',[1 1 hindex],[Inf Inf 1]))'),'linear');
F.v10.ig  = griddedInterpolant(iX,iY,flipud(single(ncread(fn_s,'v10',[1 1 hindex],[Inf Inf 1]))'),'linear');
F.u10.ig  = griddedInterpolant(iX,iY,flipud(single(ncread(fn_s,'u10',[1 1 hindex],[Inf Inf 1]))'),'linear');
F.tcc.ig  = griddedInterpolant(iX,iY,flipud(single(ncread(fn_s,'tcc',[1 1 hindex],[Inf Inf 1]))'),'linear');

%% for when the 2m stuff is available
if iWhich == 0 | iWhich == +1
  try
    %% 2m air and dewpoint temperatures
    F.t2m.ig  = griddedInterpolant(iX,iY,flipud(single(ncread(fn_2m,'t2m',[1 1 hindex],[Inf Inf 1]))'),'linear');
    F.d2m.ig  = griddedInterpolant(iX,iY,flipud(single(ncread(fn_2m,'d2m',[1 1 hindex],[Inf Inf 1]))'),'linear');
  end
end

%% for when the OLR stuff is available
if iWhich == 0 | iWhich == -1
%% for when the OLR stuff is available
  try
    %% olr and olr_clr, ilr and ilr_clr
    F.olr.ig     = griddedInterpolant(iX,iY,flipud(single(ncread(fn_olr,'mtnlwrf',[1 1 hindex],[Inf Inf 1]))'),'linear');
    F.olr_clr.ig = griddedInterpolant(iX,iY,flipud(single(ncread(fn_olr,'mtnlwrfcs',[1 1 hindex],[Inf Inf 1]))'),'linear');
    F.ilr.ig     = griddedInterpolant(iX,iY,flipud(single(ncread(fn_olr,'msnlwrf',[1 1 hindex],[Inf Inf 1]))'),'linear');
    F.ilr_clr.ig = griddedInterpolant(iX,iY,flipud(single(ncread(fn_olr,'msnlwrfcs',[1 1 hindex],[Inf Inf 1]))'),'linear');
  end
end

%% from ~/MATLABCODE/matlib/rtp_prod2/grib/grib_interpolate_era.m
%% need to check here for existence of ci vs si and read accordingly
try
    F.ci.ig   = griddedInterpolant(iX,iY,flipud(single(ncread(fn_s,'ci',[1 1 hindex],[Inf Inf 1]))'),'linear');
catch
    F.ci.ig   = griddedInterpolant(iX,iY,flipud(single(ncread(fn_s,'siconc',[1 1 hindex],[Inf Inf 1]))'),'linear');
end

%F.tcwv.ig = griddedInterpolant(iX,iY,flipud(single(ncread(fn_s,'tcwv',[1 1 hindex],[Inf Inf 1]))'),'linear');
%F.msl.ig  = griddedInterpolant(iX,iY,flipud(single(ncread(fn_s,'msl',[1 1 hindex],[Inf Inf 1]))'),'linear');

F.h_longitude = ncread(fn_h,'longitude');
F.h_latitude  = ncread(fn_h,'latitude');
F.levid       = ncread(fn_h,'level');
F.h_time      = ncread(fn_h,'time',[hindex],[1]);
F.h_mtime     = datenum(1900,0,0,double(F.h_time),0,0);

[X,Y] = ndgrid(F.h_latitude,F.h_longitude);
iX = flipud(X); iY = flipud(Y);

t = permute(single(ncread(fn_h,'t',[1 1 1 hindex],[Inf Inf Inf 1])),[2,1,3]);
for i=1:length(F.levid)
   F.t(i).ig = griddedInterpolant(iX,iY,flipud(squeeze(t(:,:,i))),'linear');
end
clear t

ciwc = permute(single(ncread(fn_h,'ciwc',[1 1 1 hindex],[Inf Inf Inf 1])),[2,1,3]);
for i=1:length(F.levid)
   F.ciwc(i).ig = griddedInterpolant(iX,iY,flipud(squeeze(ciwc(:,:,i))),'linear');
end
clear ciwc

cc = permute(single(ncread(fn_h,'cc',[1 1 1 hindex],[Inf Inf Inf 1])),[2,1,3]);
for i=1:length(F.levid)
   F.cc(i).ig = griddedInterpolant(iX,iY,flipud(squeeze(cc(:,:,i))),'linear');
end
clear cc

q = permute(single(ncread(fn_h,'q',[1 1 1 hindex],[Inf Inf Inf 1])),[2,1,3]);
for i=1:length(F.levid)
   F.q(i).ig = griddedInterpolant(iX,iY,flipud(squeeze(q(:,:,i))),'linear');   
end
clear q

o3 = permute(single(ncread(fn_h,'o3',[1 1 1 hindex],[Inf Inf Inf 1])),[2,1,3]);
for i=1:length(F.levid)
   F.o3(i).ig = griddedInterpolant(iX,iY,flipud(squeeze(o3(:,:,i))),'linear');
end
clear o3

clwc = permute(single(ncread(fn_h,'clwc',[1 1 1 hindex],[Inf Inf Inf 1])),[2,1,3]);
for i=1:length(F.levid)
   F.clwc(i).ig = griddedInterpolant(iX,iY,flipud(squeeze(clwc(:,:,i))),'linear');
end
clear clwc

% whos
