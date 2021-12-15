function [sst] = readecmwfsst_nearest(filename)

% improved from /asl/matlab/gribtools/readecmwfsst_nearest.m

% function [sst] = readecmwfsst_nearest(filename)
%
% Routine to read in a 60 or 91 level ECMWF file and return SST
% for the closest grid points to the specified (lat,lon) locations.
%
% Input:
%    filename : (string) complete ECMWG GRIB file name
%    lat : (1 x nprof) latitudes (degrees -90 to +90)
%    lon : (1 x nprof) longitude (degrees, either 0 to 360 or -180 to 180)
%
% Output:
%    sst  : (1 x nprof) sea surface temperature (actually ECMWF "SKT")
%
% Note: uses external routines: p60_ecmwf.m, p91_ecmwf.m, readgrib_inv.m,
%    readgrib_rec.m, as well as the "wgrib" program.
%

% Created: 10 Apr 2007, Scott Hannon - based on readecmwf_nearest.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Edit this section as needed

addpath /asl/matlab/gribtools  % for readgrib_inv.m, readgrib_rec.m

[lat,lon] = meshgrid(-90:0.25:90, 0:0.25:359.75);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check the GRIB file exists
d = dir(filename);
if (length(d) ~= 1)
   disp(['did not find GRIB file=' filename]);
   return
end

%filename = keepNfiles(filename,2);
filename = keepNfiles(filename,1);

% Determine 60 or 91 level based on file size.
% Note: typical file sizes are 193E+6 for 60 lev and 1160E+6 for 91 lev.
if (d.bytes < 200E+6)
   % Number of ECMWF hybrid levels
   nlev = 60;
   % Number of ECMWF latitude points
   nlat = 361;  % -90:0.50:90
   nlon = 720;  %   0:0.50:359.50
else
   % Number of ECMWF hybrid levels
   nlev = 91;
   % Number of ECMWF latitude points
   nlat = 721;  % -90:0.25:90
   nlon = 1440; %   0:0.25:359.75
end


%%%%%%%%%%%%%%%%%%%
% Check lat and lon
%%%%%%%%%%%%%%%%%%%

nprof = length(lat);
if (length(lon) ~= nprof)
   disp('Error: lon and lat are different sizes!');
   return
end

% Latitude must be between -90 (south pole) to +90 (north pole)
ii = find(lat < -90 | lat > 90);
if (length(ii) > 0)
   disp('Error: latitude out of range!')
   ii
   lat(ii)
   return
end

% Note: longitude can be either 0 to 360 or -180 to 180 
ii = find(lon < -180 | lon > 360);
if (length(ii) > 0)
   disp('Error: longitude out of range!')
   ii
   lon(ii)
   return
end

% Convert any negative longitudes to positive equivalent
xlon =lon;
ii = find( lon < 0 );
xlon(ii) = 360 + lon(ii);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convert (lat,lon) to fractional indices into a 2-D ECMWF grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

halfnlatp1 = round( 0.5*(nlat-1) + 1 );
nlonp1 = nlon + 1;
latres = round( (nlat-1)/180 );
lonres = round( nlon/360 );

% Convert latitude
glat = halfnlatp1 - latres*lat;
ii = find(glat < 1);
glat(ii) = 1;
ii = find(glat > nlat); % impossible except for tiny precision errors
glat(ii) = nlat;

% Convert longitude
glon = 1 + lonres*xlon;
ii = find(glon < 1);
glon(ii) = 1;
ii = find(glon >= nlonp1);  % Note: nlonp1 is 360=0 degrees
glon(ii) = 1;

 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find the single nearest 2-D grid point
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lon grid
iglon = floor( glon );
dg = glon - iglon;
ii = find( dg > 0.5);
iglon(ii) = iglon(ii) + 1;
ii = find(iglon == nlonp1);  % non-existant grid nlonp1 = grid 1 (0 deg)
iglon(ii) = 1;

% Lat grid
% Note: max(glat) corresponds to min(lat) and vice versa
iglat = floor( glat );
dg = glat - iglat;
ii = find( dg > 0.5);
iglat(ii) = iglat(ii) + 1;
clear ii dg


%%%%%%%%%%%%%%%%%%%
% 1-D ECMWF indices
%%%%%%%%%%%%%%%%%%%
% Note: in MATLAB, a 2-D (nrow, ncol) matrix is equivalent to a 1-D vector
% (nrow*ncol), with index translation 1-D index=irow + nrow*(icol-1).
%%%
% Old code
% i1D = iglat + nlat*(iglon - 1);
%%%
% C switches rows/columns compared to FORTRAN?
i1D = iglon + nlon*(iglat - 1);
clear iglat iglon

%%%
%% All needed grid points (without repeats)
%i1Dneed = unique(i1D);
%nneed = length(i1Dneed);
%
%% Create a lookup table to translate output index into index in "i1Dneed"
%i1Dtable = zeros(361*720,1);
%i1Dtable(i1Dneed) = 1:nneed;
%
%% Indices for each output profile in "i1Dneed"
%indo = i1Dtable(i1D);
%clear i1D i1Dtable
%%%

%%%
%spres=d(i1Dneed)/100;  % Divide by 100 to convert Pa to mb
%prof.spres = spres(indo);
%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine record numbers for profile parameters in GRIB file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get an inventory of the GRIB file
[rec,param,level] = readgrib_inv(filename);


% Parameter "SKT" skin temperature (K)
% Note: "sfc"
iparam = strcmp('SKT',param);
ii = find(iparam == 1);
if (length(ii) ~= 1)
   disp('did not find SKT in GRIB inventory');
   return
end
irec_SKT = ii;


junk = readgrib_rec(filename,irec_SKT);
sst = junk(i1D);

%%% end of function %%%


%%%%%%%%%%%%%%%%%%%%
delete(filename)