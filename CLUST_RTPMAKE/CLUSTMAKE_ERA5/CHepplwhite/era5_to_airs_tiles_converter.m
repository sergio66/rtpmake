% era5_to_airs_tiles_converter

% combine ERA5 tiles (being 0.25 x 0.25 deg ) to airs.tiles (3 x 5 deg)
% and do some stats for selected airs.tiles
% Eastings are -ve WEST of GM, +ve EAST of GM.
% AIRS.tiles to test: (Eastings first Northings second!)
% 1. seasonal lonbin 26:27, latbin 46:47 (latB:+35.75:+38.50  lon: -50:-55.0 )
% 2. clear    lonbin 54, latbin 24
% 3. cloudy   lonbin:67, latbin 35
airs.test_tiles = [27, 47; 
                   54, 24; 
		   67, 35];   %(lonbin, latbin) (N,:) N=1,2,3

% Choose with test tile (1:seasonal, 2:clear, 3:cloudy)
ibin = 1;

% Choose with ERA to analyse (era5 or erai)
ERA = 'erai';

% AIRS.tile boundaries: (72 x 64 lon x lat tiles)
load('/home/motteler/shome/airs_tiling/latB64.mat','latB2');  % latB2 x 65
airs.latB2 = latB2; clear latB2;
airs.lonB2 = [-180:5:180];

% ERA5.tiles: (1440 x 720 lon x lat tiles)
e5.home = '/asl/models/era5/2016/'; 
e5.ldir = dir([e5.home '*/*_lev_test.nc']);
e5.sdir = dir([e5.home '*/*_sfc.nc']);
e5.latB = [90:-0.25:-89.75];
e5.lonB = [-180:0.25:179.75];       % NB after shifting native ERA5

% ERA.interim tiles: (480 x 241 lon x lat tiles) lat:1st NP (90), last SP (-90).
% Long: [0:0.75:359.25] 
ei.home = '/asl/models/era/2016/'; 
ei.ldir = dir([ei.home '*/*_lev.nc']);
ei.sdir = dir([ei.home '*/*_sfc.nc']);
ei.latB = [90:-0.75:-90];
ei.lonB = [-180:0.75:179.25];

% Choose which ERA and Load some requested fields for a 16-day period
switch ERA
  case 'era5'
    dk = 8;
    gshift = 720;
    er = e5;
    er.skt  = single(zeros(1440,721,16*4));    % lon x lat x time.
    er.clwc = single(zeros(1440,721,60,16*4));
    er.ciwc = single(zeros(1440,721,60,16*4));
    er.t    = single(zeros(1440,721,60,16*4));
    er.q    = single(zeros(1440,721,60,16*4));
  case 'erai'
    dk = 4;
    gshift = 240;
    er = ei;
    er.skt  = single(zeros(480,241,16*4));    % lon x lat x time.
    er.clwc = single(zeros(480,241,60,16*4));
    er.ciwc = single(zeros(480,241,60,16*4));
    er.t    = single(zeros(480,241,60,16*4));
    er.q    = single(zeros(480,241,60,16*4));
end
er.tim = [];

% Match ERA to selected airs tile
iit = find(er.latB >= 35.75 & er.latB <= 38.50);
iin = find(er.lonB >= -50.0 & er.lonB <= -45.0);
iit = find(er.latB >= airs.latB2(46) & er.latB <= airs.latB2(47));
iin = find(er.lonB >= airs.lonB2(27) & er.lonB <= airs.lonB2(28));
% or use all grid points
iit = find(er.latB >= -90  & er.latB <= 90);
iin = find(er.lonB >= -180 & er.lonB <= 180);

k=1;
for i=1:16
  sname =  [er.sdir(i).folder '/' er.sdir(i).name];
  lname =  [er.ldir(i).folder '/' er.ldir(i).name];

  junk  = ncread(sname,'skt');
  er.skt(iix,iiw,k:k+dk-1)   = junk(iix,iiw,:);

  junk  = ncread(sname,'time');
  ntime = datetime(double(junk)*3600, 'ConvertFrom', 'epochtime', ...
         'Epoch', '01-Jan-1900');
  er.tim = [er.tim; ntime]; 

  junk = ncread(lname,'clwc');
  er.clwc(iin,iit,:,k:k+dk-1) = junk(iin,iit,:,:); 

  junk = ncread(lname,'ciwc');
  er.ciwc(iin,iit,:,k:k+dk-1) = junk(iin,iit,:,:); 
  
  junk = ncread(lname,'q');
  er.q(iin,iit,:,k:k+dk-1) = junk(iin,iit,:,:);

  junk = ncread(lname,'t');
  er.t(iin,iit,:,k:k+dk-1) = junk(iin,iit,:,:);
  
  if(k == 1)
    er.lat   = ncread(fname,'latitude');
    er.lon   = ncread(fname,'longitude');
  end

  k = k + dk;
  fprintf(1,'.')
end

% Rotate 180-deg longitude to center on Greenwich Meridian
er.lon = er.lon - 180;
er.skt = circshift(er.skt, gshift,1);

% Get stats for selected tiles from airs.tile.seasonal (+38.5,-50) lat,lon
% 
er.clwc_std = squeeze( nanstd(er.clwc(iin,iit,:,:),[],[1,2,4]) );
er.clwc_ave = squeeze( nanmean(er.clwc(iin,iit,:,:),[1,2,4]) );
for i = 1:60
  ugh = er.clwc(iin,iit,i,:);
  er.clwc_max(i) = max(ugh(:));
  er.clwc_min(i) = min(ugh(:));
end
%
er.t_std = squeeze( nanstd(er.t(iin,iit,:,:),[],[1,2,4]) );
er.t_ave = squeeze( nanmean(er.t(iin,iit,:,:),[1,2,4]) );
for i = 1:60
  ugh = er.t(iin,iit,i,:);
  er.t_max(i) = max(ugh(:));
  er.t_min(i) = min(ugh(:));
end
%
er.q_std = squeeze( nanstd(er.q(iin,iit,:,:),[],[1,2,4]) );
er.q_ave = squeeze( nanmean(er.q(iin,iit,:,:),[1,2,4]) );
for i = 1:60
  ugh = er.q(iin,iit,i,:);
  er.q_max(i) = max(ugh(:));
  er.q_min(i) = min(ugh(:));
end

% Get column sum to find profile with smallest amount  

er.clwc_colsum = squeeze(sum(er.clwc(iin,iit,:,:),3));
er.ciwc_colsum = squeeze(sum(er.ciwc(iin,iit,:,:),3));
in.clwc_min = find(er.clwc_colsum(:) < 1E-19);
in.ciwc_min = find(er.ciwc_colsum(:) < 1E-19);

% and find hottest skt scenes
junk = er.skt(iin,iit,:)
[B, IB] = maxk(junk(:),20);
[I1, I2, I3] =  ind2sub( [480 241 64], IB);
[I1, I2, I3] =  ind2sub( [7 4 64],IB);
     er.skt(iix(5),iiw(3),12)   % I1(1), I2(1), I3(1)
     for J=1:10 disp(num2str(er.skt(iin(I1(J)),iit(I2(J)),I3(J))) ); end
figure;plot(er.clwc_colsum(:),'.'); hold on;plot(er.ciwc_colsum(:),'.')
plot(IB,zeros(30,1),'kd')


% Plot it. Map it
figure; plot(er.q_ave,[1:60],'.-', er.q_max,[1:60],'.-');hold on
  plot(er.q_min,[1:60],'.-', er.q_std,[1:60],'.-');
  set(gca,'YDir','Reverse');grid on;legend('mean','max','min','std');
  ylabel('Level'); title('ERA-I q for era tiles in airs tile')
  
