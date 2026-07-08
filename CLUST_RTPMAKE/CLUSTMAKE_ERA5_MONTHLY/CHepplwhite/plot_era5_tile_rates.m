% plot_era5_tile_rates

%
%
%
%
%

% for plotting and reference load AIRS tile boundaries
load('/home/motteler/shome/airs_tiling/latB64.mat','latB2');  % latB2 x 65
airs.latB = latB2; clear latB2;
airs.lonB = [-180:5:180];

% get file listing 
d.home = '/home/chepplew/data/rates_anomalies/tiled/era5_mon/fits/';
d.list = dir([d.home 'era5_tile_fits_lonbin*_latbin*.mat']);

% Check all tiles present
lonbin = []; latbin = [];
for i = 1:length(d.list)
 junk = strsplit(d.list(i).name,{'_','.'});
 lonbin(i) = str2num(junk{4}(7:8));
 latbin(i) = str2num(junk{5}(7:8));
end 

% scatter(lonbin,latbin,4)

% load requested field

skt = [];
for ip = 1:length(d.list)
  fn = [d.list(ip).folder '/' d.list(ip).name];

  x = load(fn);
  skt(lonbin(ip), latbin(ip)) = x.fits.cen.skt_bcoef(2);

end

% MAP trend
addpath /asl/matlib/maps          % aslmap
addpath /home/strow/Matlab/Extra/
addpath /asl/matlib/plotutils
load llsmap5

mopts.color = 'k';
mopts.title = 'ERA.5 SKT trend K/yr';
mopts.caxis = [-Inf Inf];
%mopts.caxis = [-0.25 0.25];
mopts.cmap  = llsmap5;
mopts.titlesize = 14;

zdata = skt;

fign = 3;
fh = aslmap(fign, airs.latB, airs.lonB, zdata, [-90 90], [-180 180], mopts);
