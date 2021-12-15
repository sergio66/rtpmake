function [era] = load_era5_monthly(date1,date2,lonbin,latbin)


% NB: longitude +ve is East of Greenwich Meridian
% acceptable variable names: (lev) longitude latitude, level, time,
%   cc, o3, ciwc, clwc, t, q, 
% (sfc) skt, v10, u10, tcc, sp, siconc, time, longitude, latitude, time
% grids are 1440 x 721.
% era latitude  = [-90:0.25:90];
% era longitude = [0:0.25:359.75];
% Joel time period: Jan 2003 thru Dec. 2017
%
% selected airs.tile: lonbin 27, latbin 47 ([],[-50,-48]) (lat:+38.5  lon: -50.0)
% rtv.lon(47,26:27,1)= -50   -48
% rtv.lat(47:48,27,1)= 38.5000 41.2500
% -> era.tile lonbin: 521    latbin:207
% Eastings are -ve WEST of GM, +ve EAST of GM.
% AIRS.tiles to test: (Eastings first Northings second!)
% 1. seasonal lonbin 26:27, latbin 46:47 (
% 2. clear    lonbin 54, latbin 24
% 3. cloudy   lonbin:67, latbin 35
airs.test_tiles = [27, 47;   % seasonal (lon:-50:-55.0 lat:+35.75:+38.50)
                   54, 24;   % clear
		   67, 35;   % cloudy (lon:+150:+155, lat:+5.5:+8.25)
		   64, 34;   % Equator (lonbin, latbin) (N,:) N=1,2,3
		   36,  3;   %
		   36, 62];  %



addpath /asl/matlib/time
addpath /home/chepplew/myLib/matlib
addpath /home/strow/Matlab/Math                  % Math_tsfit_lin_robust.m
addpath /home/chepplew/myLib/matlib/math 

% redirect matlab stdout to file
%foH = fopen('/home/chepplew/logs/matlab/stdout.log','w')
foH = 1;

switch nargin
  case 0
    % Define date range to match the retrieval dataset.
    sdate = '2002/09/08';
    edate = '2020/10/29';
    D1 = datenum(sdate,'yyyy/mm/dd');
    D2 = datenum(edate,'yyyy/mm/dd');
  case 2
    % check validity
    try
      D1 = datenum(date1,'yyyy/mm/dd');
      D2 = datenum(date2,'yyyy/mm/dd');
    catch
      error('Incorrect Date Format')
      return
    end
  case 4
    % check lonBin latBin
    if(~ismember(lonbin,[1:72]) | ~ismember(latbin,[1:64]) )
      error('latbin or lonbin out of range'); return
    end 
    try
      D1 = datenum(date1,'yyyy/mm/dd');
      D2 = datenum(date2,'yyyy/mm/dd');
    catch
      error('Incorrect Date Format')
      return
    end
  otherwise
    error('invalid number of arguments [0 or 2]');
end

% ERA5 tile boundaries to find inda, indn(BEFORE longitude shift as loaded):
era.lonB  = [0:0.25:359.75];       % Before
era.lonBc = [-180:0.25:180.00];    % AFTER circshift
era.lonC  = wrapTo180(era.lonB);
era.latB  = [90:-0.25:-90];

% airs tile boundaries: (72 x 64 lon x lat tiles)
load('/home/motteler/shome/airs_tiling/latB64.mat','latB2');  % latB2 x 65
airs.latB2 = latB2; clear latB2;
airs.lonB2 = [-180:5:180];

% Specified AIRS.tile region to load:
if(nargin ~= 4)
  iwb = 3;
  lonbin = airs.test_tiles(iwb,1);
  latbin = airs.test_tiles(iwb,2);
end
  disp(['lonbin: ' num2str(lonbin) ' latbin: ' num2str(latbin)]);

% Check specified tile for levels ingest:
%[en, elonbin] = deal(521);
%[et, elatbin] = deal(207);
% era5.tiles required to cover airs.tile lonB:[-50,-48] latB:[38.5,41.25].
%[en, elonbin] = deal(514:532);
%[et, elatbin] = deal(196:207);
% era5.tiles required to cover airs.tile lonB:[85:90]. latB:[-24.75:-22.0]
inda = find(era.latB >= airs.latB2(latbin)  & era.latB <= airs.latB2(latbin+1) );
if(airs.lonB2(lonbin) <= 0 & airs.lonB2(lonbin+1) <= 0)
  indn = find(era.lonB >= 360+airs.lonB2(lonbin) & ...
              era.lonB <= 360+airs.lonB2(lonbin+1) );
end
if(airs.lonB2(lonbin) >= 0 & airs.lonB2(lonbin+1) >= 0)
  indn = find(era.lonB >= airs.lonB2(lonbin) & ...
              era.lonB <= airs.lonB2(lonbin+1) );
end
%  indn = find(era.lonBc >= airs.lonB2(lonbin) & era.lonBc <= airs.lonB2(lonbin+1));
  
% used to size output arrays
iinda = [1:length(inda)];
iindn = [1:length(indn)];

%{
% sanity mapping checks
airs.latB2(latbin:latbin+1)
era.lat(inda)
airs.lonB2(lonbin:lonbin+1)
era.lon(indn)

addpath /asl/matlib/aslutil
figure; simplemap(rot90(era.skt(:,:,1)))       % rot90 since long x lat
xlat = rot90(repmat([89.5:-1:-89.5],360,1));
xlon = repmat([-179.5:1:179.5],180,1);
figure;simplemap(xlat, xlon, (squeeze( airs.skt_a(:,:,1))) );caxis([220 320]);

hf = pcolor(squeeze(airs.skt_d(:,:,1)) );  set(hf, 'EdgeColor', 'none'); colorbar;


%}
% check valid variable names of interest:
allvars = {};

disp('locating era5 data')
d.home = '/asl/models/era5_avg/';
d.sdir = [];
d.ldir = [];
for iyr = 2002:2020
  d.sdir  = [d.sdir; dir([d.home '/' num2str(iyr) '/*_sfc.nc'])];  
  d.ldir  = [d.ldir; dir([d.home '/' num2str(iyr) '/*_lev.nc'])];  

end

disp('Loading levels fields')
era.ednum = [];
era.tz    = single(zeros(length(indn),length(inda),37,length(d.ldir)));
era.q     = single(zeros(length(indn),length(inda),37,length(d.ldir)));
era.o3    = single(zeros(length(indn),length(inda),37,length(d.ldir)));
era.clwc  = single(zeros(length(indn),length(inda),37,length(d.ldir)));
era.ciwc  = single(zeros(length(indn),length(inda),37,length(d.ldir)));
era.cc    = single(zeros(length(indn),length(inda),37,length(d.ldir)));
k = 1;
for i = 1:length(d.ldir)
% i = 110;
  fname = [d.ldir(i).folder '/' d.ldir(i).name] ;

  if(k==1) 
    era.plevs = ncread(fname,'level');
    era.lat   = ncread(fname,'latitude'); 
    era.lon   = ncread(fname,'longitude');
  end
  ntime = ncread(fname,'time');
  % ERA5 monthly means are centered middle month
  junk = datetime(double(ntime(4))*3600, 'ConvertFrom', 'epochtime', ...
         'Epoch', '15-Jan-1900');
  era.ednum = [era.ednum; datenum(junk)];    

  junk = ncread(fname,'t');
  era.tz(iindn,iinda,:,k)   = nanmean(junk(indn,inda,:,:),4);
  junk = ncread(fname,'q');
  era.q(iindn,iinda,:,k)    = nanmean(junk(indn,inda,:,:),4);
  junk = ncread(fname,'o3');
  era.o3(iindn,iinda,:,k)   = nanmean(junk(indn,inda,:,:),4);
  junk = ncread(fname,'clwc');
  era.clwc(iindn,iinda,:,k) = nanmean(junk(indn,inda,:,:),4);
  junk = ncread(fname,'ciwc');
  era.ciwc(iindn,iinda,:,k) = nanmean(junk(indn,inda,:,:),4);
  junk = ncread(fname,'cc');
  era.cc(iindn,iinda,:,k)   = nanmean(junk(indn,inda,:,:),4);
  
  
  k = k+1;    
  fprintf(1,'.')
end
fprintf(1,'\n')

disp('Loading surface fields')
era.sdnum = [];
era.skt   = single(zeros(length(indn),length(inda),length(d.sdir))); 
%era.lat   = []; era.lon = [];
k = 1;

% latitude runs from [90:0.25:-90]. longitude from [0:0.25:359.75]
for i = 1:length(d.sdir)
% i = 110;
  fname = [d.sdir(i).folder '/' d.sdir(i).name];
  ntime = ncread(fname,'time');

  % ERA5 monthly means are centered middle month
  junk = datetime(double(ntime(4))*3600, 'ConvertFrom', 'epochtime', ...
         'Epoch', '15-Jan-1900');
  era.sdnum = [era.sdnum; datenum(junk)];    
  
  %if (i == 1) 
  %  era.lat  = [era.lat; ncread(fname, 'latitude')];
  %  era.lon  = [era.lon; ncread(fname, 'longitude')];
  %end
  junk = ncread(fname,'skt');
  era.skt(iindn,iinda,k)  = nanmean(junk(indn, inda,:),3);
  
  k = k+1;   
  fprintf(1,'.')
end

whos *dnum lat lon skt

era.inda   = inda;
era.indn   = indn;
era.lonbin = lonbin;
era.latbin = latbin;

%era




% ------- END of FUNCTION ------------

%{
% translate the field longitudes from 0:360 to -180:180
era.skt = circshift(era.skt,720,1);
era.lon = era.lon - 180;

% Subset time period if requested:
switch nargin
  case 0
    indt  = find(era.fdnum);
  case 2
    indt  = find(era.sdnum >= D1 & era.sdnum <= D2);
end


 Pick selected tile: lat:+38.5  lon: -50.0
% find(era.lat > 38 & era.lat < 39)     -> 207
% find(era.lon > -50.5 & era.lon < -50) -> 520
% Lat: -22.4, Lon: 87.25
% find(era.lat > -22.5 & era.lat < -22) ->   450
% find(era.lon > 87 & era.lon < 87.5)   ->  1070
figure; hold on;
 plot(datetime(era.fdnum,'convertfrom','datenum'),squeeze(era.skt(1070,450,:)),'.-')






% Get trend at original grid - too slow: need to parallelize)
linr  = zeros(1440,721,1);
rerr  = zeros(1440,721,1);
anom  = zeros(1440,721,length(indt));
xdat  = era.sdnum(indt);
warning('off','all')
parfor i = 1:1440
  for j = 1:721
    ydat = squeeze(era.skt(i,j,indt)) - nanmean(era.skt(i,j,indt),3);
    try
      [B, Be, Bstats] = Math_tsfit_lin_robust(xdat-xdat(1)+1,ydat,2);
      %[B, Be, stats] = Math_tsfit_lin_robust(xdat_daily-xdat_daily(1)+1, ydat_daily,1);

      %yhat = B(2)*(days(ndt1:ndt2)-days(ndt1))./365 + B(1);
      linr(i, j,:) = B(2);
% Lag-1 uncertainty estimates
      lxc    = xcorr(Bstats.resid,1,'coeff');
      lag1   = lxc(1);
      lag1_c = sqrt( (1 + lag1) ./ ( 1 - lag1));
      rerr(i,j) = Bstats.se(2).*lag1;
% to compare with SM anom retrievals add back linear term
      dr = (xdat - xdat(1)+1)./365. *B(2);
      anom(i,j,:) = Bstats.resid + dr;
    catch ME
      fprintf(1, '%s\n', ME.message);
      continue; 
    end
  end
  if(~mod(i,20)) fprintf(1,'.'); end
end
warning('on','all')

% ---------------------------
% fit and anomaly levels field

warning('off','all')
olinfit = zeros(1440,721,37,1);
anom    = [];
fiterr  = zeros(1440,721,37,1);
indt    = find(era.ednum);
xdat    = era.ednum(indt);

i = indn(1)
  j = inda(1)
    for ilev = 1:37
      ydat = squeeze(era.tlev(i,j,ilev,indt)) - nanmean(era.tlev(i,j,ilev,indt),4);
      try
        [B, Be, Bstats] = Math_tsfit_lin_robust(xdat-xdat(1)+1,ydat,2);
        olinfit(i, j, ilev,:) = B(2);
        % to compare with SM anom retrievals add back linear term
        dr = (xdat - xdat(1)+1)./365. *B(2);
        anom(i,j,ilev,:) = Bstats.resid + dr;
      catch ME
        fprintf(1, '%s\n', ME.message);
        continue; 
      end
    end
% end
% end
warning('on','all')



% test tile with olinfit > 0.31 K/yr (nr. Severny Island Ru:
% plot(datetime(era.fdnum,'convertfrom','datenum'),squeeze(era.skt(955,74,:)),'.-');grid on

% Interpolate era to airs tiles
% 1. load airs.tiled grid:
gridfile = '/home/motteler/repos/airs_tiling/latB64.mat';
load(gridfile, 'latB2');
% build longitude array
lonB2 = -180:5:180;
% convert from boundaries to centers
tclat = (latB2(2:end) + latB2(1:end-1) ) / 2;
tclon = (lonB2(2:end) + lonB2(1:end-1) ) / 2;

% grid of tile center interpolation points
[aX,aY] = meshgrid(tclon, tclat);

% era geo are centers, to meshgrid.
eclat = (era.lat(2:end) + era.lat(1:end-1) ) / 2;
eclon = (era.lon(2:end) + era.lon(1:end-1) ) / 2;
eclon = [eclon; 179.875];
[eX,eY] = meshgrid(eclon, eclat);

% Interpolate ERA to airs.tiles and add fields to returning structure.


%tdata = nanmean(era.skt,3)';
tdata = era.skt(:,1:720,99)';
tdata = olinfit(:,1:720)';
%datetime(fdnum(99),'convertfrom','datenum')
interp_td = interp2(eX, eY, tdata, aX, aY, 'linear');
era.skt2  = interp_td;

% global map
addpath /asl/matlib/maps          % aslmap 
addpath /home/strow/Matlab/Extra/
addpath /asl/matlib/plotutils
load llsmap5

[SZ, SS, SFG] = smoothn(dlr);
 
mopts.color = 'k';
mopts.title = 'ERA5 month.mn skt (arb timestep) K';
mopts.caxis = [-Inf Inf];
mopts.caxis = [-0.25 0.25];
mopts.cmap  = llsmap5;
mopts.titlesize = 14;

fign = 2;
tdata = squeeze(era.skt(:,:,1));
xlat = flipud(double(era.lat));
xlon = double(era.lon);
 %fh = equal_area_map(fign, gis.latB2, gis.lonB2, SZ, txtstr);
 fh = aslmap(fign, xlat, xlon, rot90(tdata(1:end-1,1:end-1)), ...
            [-90 90], [-180 180], mopts);

mopts.title='ERA 2010.Oct mean skt (K)';
  fh = aslmap(fign, double(lat), double(elon), double(Y(1:720,:)), [-90 90], [-180 180], mopts);

fign=3
mopts.title = 'ERA5 month SKT lin.rate K/yr';
mopts.caxis = [-0.3 0.3];
tdata = rot90(linr(1:end-1,1:end-1));
fh = aslmap(fign, xlat, xlon, tdata, [-90 90], [-180 180], mopts);

fclose(foH)
fign=5;
lonB3 = lonB2+180;
fh=aslmap(fign,latB2,lonB3,interp_td,[-90 90],[-180 180],mopts);


% ANOMALY
figure; hold on;
for i=1:length(indn);
  for j = 1:length(inda)
    plot(datetime(era.sdnum(indt),'convertfrom','datenum'), ...
         squeeze(anom(indn(i),inda(j),:)),'.-');
  end
end
  title(['ERA5 SKT all tiles anom bin: ' ...
   sprintf('%2d,%2d',lonBin,latBin)] )
  ylabel('temp (K)');grid on;

junk = squeeze(nanmean(anom(indn,inda,:),[1 2]));
figure;plot(datetime(era.sdnum(indt),'convertfrom','datenum'),junk,'.-')
  grid on; ylabel('temp (K)'); 
  title(['ERA5 skt mean anom. bin: ' ...
         sprintf('%2d,%2d',lonBin,latBin)] );





 junk = squeeze(anom(i,j,:,:));
 hf=pcolor(datetime(era.ednum,'convertfrom','datenum'),era.plevs,junk);
 set(hf, 'EdgeColor', 'none')
set(gca,'Yscale','log');set(gca,'YDir','reverse')
colormap(llsmap5);  colorbar
title('ERA.5 T month mean anom (K)')
caxis([-5 5]);  ylim([1 1000])

%}
