function F = grib_interpolateERA5sfc(fn_s,hindex);

% quick version of grib_interpolate_era5x.m
% 
% Inputs: fn_s, hindex
%         Netcdf files containing grib1 and grib2 
%         data respectively.  This code assumes the grib1 data
%         is surface data, and the grib2 data is hybrid ecmwf data.
%         hindex is the hour index (8 total) in the era5 file

% Output: F, a structure with housekeeping data (lat, lon of grid) and
%         interpolants for interpolating to arbitrary lat/lon positions.
%
% Presently the interpolant is set to be linear, you can also 
% select nearest-neighbor, cubic, etc.  I may make this an option
% in the future.  Minor changes are needed to use this for ERA, etc.
%
% L. Strow, June 11, 2014

F.s_longitude = ncread(fn_s,'longitude');
F.s_latitude  = ncread(fn_s,'latitude');
%F.s_time      = ncread(fn_s,'time',[hindex],[1]);
%F.s_mtime     = datenum(1900,0,0,double(F.s_time),0,0);

F.s_longitude = wrapTo180(F.s_longitude);

[X,Y] = ndgrid(F.s_latitude,F.s_longitude);
iX = flipud(X); iY = flipud(Y);

F.sp.ig       = griddedInterpolant(iX,iY,flipud(single(ncread(fn_s,'sp', [1 1 hindex],[Inf Inf 1]))'/100),'linear');
F.skt.ig      = griddedInterpolant(iX,iY,flipud(single(ncread(fn_s,'skt',[1 1 hindex],[Inf Inf 1]))'),'linear');
F.v10.ig      = griddedInterpolant(iX,iY,flipud(single(ncread(fn_s,'v10',[1 1 hindex],[Inf Inf 1]))'),'linear');
F.u10.ig      = griddedInterpolant(iX,iY,flipud(single(ncread(fn_s,'u10',[1 1 hindex],[Inf Inf 1]))'),'linear');
F.tcc.ig      = griddedInterpolant(iX,iY,flipud(single(ncread(fn_s,'tcc',[1 1 hindex],[Inf Inf 1]))'),'linear');
F.pblh_nwp.ig = griddedInterpolant(iX,iY,flipud(single(ncread(fn_s,'blh',[1 1 hindex],[Inf Inf 1]))'),'linear');
F.t2m.ig      = griddedInterpolant(iX,iY,flipud(single(ncread(fn_s,'t2m',[1 1 hindex],[Inf Inf 1]))'),'linear');

