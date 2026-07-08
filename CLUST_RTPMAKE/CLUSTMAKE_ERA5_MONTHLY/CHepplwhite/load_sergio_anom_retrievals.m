function [a] = load_sergio_anom_retrievals()

% Sergio's anomaly retrievals
% sub-dirs /[01..64]/ -> /[01..72]/
%
% oem.finalrates(1:64)
% 1,2,3,4 = CO2 stemp cld1 cld2
%         next 20.      = WV(z)
%         next 20       = T(z)
%         final 20.     = O3(z)

addpath /home/chepplew/myLib/matlib
addpath /home/strow/Matlab/Math                  % Math_tsfit_lin_robust.m
addpath /home/chepplew/myLib/matlib/math


% Home of Sergio's anomaly retrievals
d.home = ['/home/sergio/MATLABCODE/oem_pkg_run/' ...
          'AIRS_gridded_TILES_Anomaly_ROSES2020/OutputAnomaly_OBS/'];

% Hardwire start date of sequence
sdate = '2002/10/20';
sdnum = datenum(sdate,'yyyy/mm/dd');

% allocate order sorted file index
osnum = cell(64,72);
ixnum = cell(64,72);

% get listing: NB: loads in order 1, 10, 100:199 2,20 200:299, 3,30, 300:387, 
% 31:99 but with missing steps.
for i = 1:64

  %for j = 1:72 
  for j = 1:72
    d.dir = dir([d.home sprintf('%02d/%02d',i,j) '/anomtest_timestep*.mat']);
    %d.dir  =  [d.dir; dir([d.home sprintf('%02d/%02d',i,j)  ...
    %           '/anomtest_timestep*.mat']) ];  %end

    inum = [];
    for k = 1:length(d.dir)
      junk  = strsplit(d.dir(k).name, {'step','.'});
      inum(k)  = str2num(junk{2});
    end
    [ixnum{i,j},osnum{i,j}] = sort(inum);
    % ixnum records actual time step

  end 
  fprintf(1,'.') 
end

% figure;hold on; for i=1:64 for j=1:72 plot(j,length(ixnum{i,j}),'.'); end; end
% for stemp use: x.oem.finalrates(2)

% sample tile: 'Seasonal' lonbin 27, latbin 47 (lat:+38.5  lon:-50.0 ).i=47;j=27

linfit = zeros(64,72);
for i = 1:64 
  for j=1:72
    d.dir = dir([d.home sprintf('%02d/%02d',i,j) '/anomtest_timestep*.mat']);
    n = 1;
    rate_skt = [];
    for k = osnum{i,j}
      %disp(d.dir(k).name);
      x = load([d.dir(k).folder '/' d.dir(k).name]);
      rate_skt(n) = x.oem.finalrates(2);
      n = n + 1;
    end

    % sample datenums
    sdnums = [sdnum + ixnum{i,j}*16];
    try
      [B, Be, stats] = Math_tsfit_lin_robust(sdnums-sdnums(1)+1,rate_skt,1);
      linfit(i,j) = B(2);
    catch ME
      fprintf(1, '%s\n', ME.message);
      continue; 
    end
  end
  fprintf(1,'.')
end  


%
% Do global map at airs.tiled grid
% 1. load airs.tiled grid:
gridfile = '/home/motteler/repos/airs_tiling/latB64.mat';
load(gridfile, 'latB2');
% build longitude array
lonB2 = -180:5:180;
% convert from boundaries to centers
tclat = (latB2(2:end) + latB2(1:end-1) ) / 2;
tclon = (lonB2(2:end) + lonB2(1:end-1) ) / 2;

% grid of tile center interpolation points
[cX,cY] = meshgrid(tclon, tclat);

% Smoothed map
addpath /asl/matlib/plotutils
load(llsmap5);
addpath /asl/matlib/maps          % aslmap 
addpath /home/strow/Matlab/Extra/
 [SZ, SS, SFG] = smoothn(dlr);
 
mopts.color = 'k';
mopts.title = 'Anomaly skt rates K/yr 2002:19';
mopts.caxis = [-Inf Inf];
mopts.caxis = [-0.25 0.25];
mopts.cmap  = llsmap5;
mopts.titlesize = 14;

 [SZ, SS, SFG] = smoothn(linfit);
 
fign = 5;
 %fh = equal_area_map(fign, gis.latB2, gis.lonB2, SZ, txtstr);
 fh = aslmap(fign, latB2, lonB2, SZ, [-90 90], [-180 180], mopts);

% ------------- get GISSv4 stemp ---------------------
gopts = mopts;
gopts.title='GSSv4 skt rates K/yr';

fign=4;
fh = aslmap(fign, gis.latB2, gis.lonB2, gis.lnr2, [-90 90], [-180 180], gopts);

% ----- get difference between SM.anom.rates to GISSv4 rates --

dlr = linfit - gis.lnr2;
[SZ, SS, SFG] = smoothn(dlr);

fign = 6;
dopts = mopts;
dopts.title='SKT sm.anom minus GSSv4 rate diff K/yr';
 fh = aslmap(fign, gis.latB2, gis.lonB2, SZ, [-90 90], [-180 180], dopts);
 
 

