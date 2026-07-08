load_tile_daily_data()

% Data are supplied in 16-day collections (22 or 23 collections per year)
% with 64 zonal bands each zone has 72 longitude tiles.
%
cd /home/chepplew/projects/rates_anomalies/tiled

% INPUT: Choose Zone band to load e.g. S90p00
zn_str = 'S90p00';
zn_str = 'N81p75';
zn_str = 'N00p00';

d.home = '/home/chepplew/data/rates_anomalies/tiling/daily_means/';
d.dir  = dir([d.home '*' zn_str '20*_lw40_ssw_daymean.mat']);

robs = [];  rdays = [];  rmn = []; nobs = [];
for ifn = 1:length(d.dir)
   filename = [d.dir(ifn).folder '/' d.dir(ifn).name];
   
   x = load(filename);
   robs   = cat(2, robs, x.robs);
   rmn    = cat(1, rmn, squeeze(nanmean(x.robs,1)) );
   rdays  = [rdays x.days]; 
   nobs   = [nobs x.itot]; 

   fprintf(1,'.')
end      
fairs = x.fairs;
xchs  = x.xchs;
lonB  = x.lonB;

whos r* nobs
