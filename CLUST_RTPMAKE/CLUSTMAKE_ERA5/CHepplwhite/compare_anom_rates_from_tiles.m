compare_anom_rates_from_tiles()
%
% Compare Anomaly Retrievals to climatology
%
% first:  rtv = load_sergio_anom_retrievals();
% second: gis = load_gisstemp();
% or      era = load_era5_monthly();
% then here: get rates from retrieved stemp

addpath /home/strow/Matlab/Math      % Math_tsfit_lin_robust.m
addpath /home/chepplew/myLib/matlib/math

>> plot(squeeze(tim(32,18,:)),squeeze(stemp(32,18,:)),'.-')

% re-zero the days in sequence.
[an, am, sz] = size(rtv.stemp);        % 64 lats, 72 lons.

% test tile:
in = 32; im = 30;

% set output arrays
linfit = zeros(an,am,4);
fiterr = zeros(an,am);

% record start end end dates of time series. prep to interp to daily.
sdnum = rtv.tim(1);
ednum = rtv.tim(end);

% convert stemp to anomaly
stemp_mn   = nanmean(rtv.stemp,3);
stemp_anom = rtv.stemp - stemp_mn;

%plot(datetime(squeeze(rtv.tim(in,im,:)),'convertfrom','datenum'),... 
%     squeeze(stemp_anom(in,im,:))+stemp>> (in,im),'.-')

for in = 1:an
  for im = 1:am
  
    xdata  = squeeze(rtv.tim(in,im,:) - rtv.tim(in,im,1));
    ydata  = squeeze(stemp_anom(in,im,:));

    try
    [B, Be, stats] = Math_tsfit_lin_robust(xdata,ydata,1);

    %yhat = B(2)*(days(ndt1:ndt2)-days(ndt1))./365 + B(1);
    %yhat = B(2)*(dn_rel./365) + B(1);

    linfit(in, im,:) = B;

% Lag-1 uncertainty estimates
      lxc    = xcorr(stats.resid,1,'coeff');
      lag1   = lxc(1);
      lag1_c = sqrt( (1 + lag1) ./ ( 1 - lag1));
      fiterr(in,im) = stats.se(2).*lag1;
      %l = xcorr(stats_asc.resid,1,'coeff');
      %lag_asc = l(1);
      %lag_asc_c = sqrt( (1 + lag_asc) ./ ( 1 - lag_asc));
      %b_err_desc(lati,loni) = stats_asc.se(2).*lag_asc;

    catch ME
      fprintf(1, '%s\n', ME.message);
      continue; 
    end

    if(~mod(i,100)) fprintf(1,'.'); end
  end
end

% do global map
addpath /home/motteler/shome/obs_stats/source    % equal_area_map
load('/asl/matlib/plotutils/llsmap5.mat');


% Get Lat and Lon boundaries of the tiles
load('/home/motteler/shome/airs_tiling/latB64.mat','latB2');  % latB2 x 65
lonB2 = [-180:5:180];

ydat = squeeze(linfit(:,:,2));
txtstr=['v2 oem stemp rates' sprintf('K/yr')];
fign = 2;
fh = equal_area_map(fign, latB2, lonB2, ydat, txtstr);
colormap(llsmap5); 
caxis([-0.2 0.2])

>>   hold on; plot(dn_rel, yhat,'-')

% plot GISS skt rates

fign = 3;
ydat = gis.lnr2;
txtstr = ['GISSv4 skt rates ' sprintf('K/yr')];
fh = equal_area_map(fign, gis.latB2, gis.lonB2, ydat,txtstr);
colormap(llsmap5); 
caxis([-0.2 0.2])

% Take difference 

dlr = squeeze(linfit(:,:,2)) - gis.lnr2;
fign = 4;
txtstr = ['OEM minus GISSv4 skt rates ' sprintf('K/yr')];
fh = equal_area_map(fign, gis.latB2, gis.lonB2, dlr,txtstr);
colormap(llsmap5); 
caxis([-0.2 0.2])

% Smoothed map
addpath /asl/matlib/maps          % aslmap 
addpath /home/strow/Matlab/Extra/
 [SZ, SS, SFG] = smoothn(dlr);
 
mopts.color = 'k';
mopts.title = 'OEM  skt rates K/yr';
mopts.caxis = [-Inf Inf];
mopts.caxis = [-0.25 0.25];
mopts.cmap  = llsmap5;
mopts.titlesize = 14;

fign = 5;
 %fh = equal_area_map(fign, gis.latB2, gis.lonB2, SZ, txtstr);
 fh = aslmap(fign, gis.latB2, gis.lonB2, SZ, [-90 90], [-180 180], mopts);
 colormap(llsmap5);
 caxis([-0.2 0.2]) 
 
% Get zonal means
ydat = squeeze(linfit(:,:,2));
%ydat = gis.lnr2;

[SZ, SS, SFG] = smoothn(ydat);
sz_mn = nanmean(SZ,2);

clf; plot(gis.lat2,sz_mn, '.-')
>> grid on; xlabel('lat (deg)'); ylabel('skt rate K/yr');title('OEM skt rates 2002:20 Zonal Means')


