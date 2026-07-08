
addpath /asl/matlib/h4tools
addpath /home/chepplew/myLib/matlib/math
addpath /home/chepplew/myLib/matlib/convert_gas_units

% for plotting and reference load AIRS tile boundaries
load('/home/motteler/shome/airs_tiling/latB64.mat','latB2');  % latB2 x 65
airs.latB = latB2; clear latB2;
airs.lonB = [-180:5:180];


d.home  = '/asl/s1/sergio/MakeAvgProfs2002_2020/';
fnpatt  = 'pall_16daytimestep_*.rtp';
d.dir   = dir([d.home fnpatt]);

% ------------------
%  Load in the data
% ------------------
for ifn = 1:length(d.dir)

  if(ifn == 1)
    [hd, ~, pd, ~] = rtpread([d.dir(ifn).folder '/' d.dir(ifn).name]);
    rlat = pd.rlat;
    rlon = pd.rlon;
    plev = pd.plevs(:,1);
  end

  junk = strsplit(d.dir(ifn).name,{'_','.'});
  dt1(ifn)  = datenum(cell2mat(junk(4:6)),'YYYYmmdd');
  dt2(ifn)  = datenum(cell2mat(junk(8:10)),'YYYYmmdd');
  iset(ifn) = str2num(junk{3});

end

[dnum iyt] = sort(dt1);
dnum = (dt1(iyt) + dt2(iyt))/2;

% plot(dnum,'.')

rtim = [];  skt = []; wv = []; o3 = []; o3 = []; rh = [];
rh1km = []; tz  = [];
k = 1;
for ifn = iyt

  [hd, ~, pd, ~] = rtpread([d.dir(ifn).folder '/' d.dir(ifn).name]);
  disp(d.dir(ifn).name)
  
  rtim(k,:)   = pd.rtime;
  skt(k,:)    = pd.stemp;
  tz(k,:,:)   = pd.ptemp;
  wv(k,:,:)   = pd.gas_1;
  co2(k,:,:)  = pd.gas_2;
  o3(k,:,:)   = pd.gas_3;

  % convert wv to relative humidity
  hd.ptype = 2;
  pd.palts(101,:) = pd.palts(100,:) - 240;
  [xh, xh_1km, colwater] = layeramt2RH(hd,pd);
  rh(k,:,:)    = xh;
  rh1km(k,:,:) = xh_1km;

  % convert o3 amt to ppmv
  np = size(pd.gas_3,2);
  [ppmvLAY,ppmvAVG,~,pavgLAY,tavgLAY,~,ppmv75,~] = layers2ppmv(hd,pd,[1:np],3);
  o3pp(k,:,:)  = ppmvLAY;
    
  k = k + 1;
  %fprintf(1,'.')
end

plevs = pd.plevs(:,100);

whos rtim skt tz wv co2 o3 o3pp rh rh1km plevs


% convert wv to fractionals
wv_mn = squeeze(nanmean(wv,1));
for i=1:nt
  for j=1:np
    for k=1:ns
      wv_fc(i,j,k) = wv(i,j,k)/wv_mn(j,k);
    end
  end
end

% select a tile
% Lookup table for job configuration (one less than tile boundaries!)
lonbins = [1:72];
latbins = [1:64];
[X,Y] = meshgrid(lonbins,latbins);

% test tile:cloudy [67,35] lonbin,latbin (N.Dakota [17,49])
ix = find(X == 67 & Y == 35)

% check lat:lon
rlon(ix)
rlat(ix)

% ------------------
% Group processing:
% ------------------
allvars = {'skt','tz','wv','co2','o3','o3pp','rh','rh1km'};
fitflds = {'coef','rerr','anom'};

kn = 6;
varn  = allvars{kn};
cfldn = [allvars{kn} '_' fitflds{1}];
rfldn = [allvars{kn} '_' fitflds{2}];
afldn = [allvars{kn} '_' fitflds{3}];

% QC
yvar  = eval(varn);
inan = find(yvar <= 0);
whos inan
yvar(inan) = NaN;


indt  = find(dnum);
xdat  = dnum(indt)'; %  - dnum(indt(1)) + 1;

% -------------------------
% Fit Profile vars
% -------------------------

if(ndims(yvar) == 3) 
  [nt np ns] = size(yvar);  % time x Z x tile
clear tmp*
warning('off','all')
tic
parfor ilay = 1:np
  warning('off','all')
  for itile = 1:ns

  ydat = squeeze(yvar(indt,ilay,itile)); %  - nanmean(yvar(ilay,itile);
    try
      [Bc, Be, Bstats] = Math_tsfit_lin_robust(xdat-xdat(1)+1,ydat,4);
      tmp1(ilay,itile,:) = Bc;
      %yhat = B(2)*(days(ndt1:ndt2)-days(ndt1))./365 + B(1);
      %linr(i, j,:) = B(2);
% Lag-1 uncertainty estimates
      lxc    = xcorr(Bstats.resid,1,'coeff');
      lag1   = lxc(1);
      lag1_c = sqrt( (1 + lag1) ./ ( 1 - lag1));
      tmp2(ilay,itile)   = Bstats.se(2).*lag1_c;
% to compare with SM anom retrievals add back linear term
      dr = (xdat - xdat(1)+1)./365. *Bc(2);
      tmp3(ilay,itile,:) = Bstats.resid + dr;
    catch ME
      fprintf(1, '%s\n', ME.message);
      %continue;
    end
  end
end
toc
warning('on','all')

fits.(cfldn) = tmp1;
fits.(rfldn) = tmp2;
fits.(afldn) = tmp3;
clear tmp*

end    % end ndims == 3

% -----------------------
% Fit Surface vars
% -----------------------
if(ndims(yvar) == 2) [nt ns] = size(yvar); np = [];   % time x tile

warning('off','all');
tic
for itile = 1:ns

  ydat = squeeze(yvar(indt,itile)); %  - nanmean(yvar(ilay,itile);
    try
      [Bc, Be, Bstats] = Math_tsfit_lin_robust(xdat-xdat(1)+1,ydat,4);
      fits.(cfldn)(itile,:) = Bc;
      %yhat = B(2)*(days(ndt1:ndt2)-days(ndt1))./365 + B(1);
      %linr(i, j,:) = B(2);
% Lag-1 uncertainty estimates
      lxc    = xcorr(Bstats.resid,1,'coeff');
      lag1   = lxc(1);
      lag1_c = sqrt( (1 + lag1) ./ ( 1 - lag1));
      fits.(rfldn)(itile)   = Bstats.se(2).*lag1_c;
% to compare with SM anom retrievals add back linear term
      dr = (xdat - xdat(1)+1)./365. *Bc(2);
      fits.(afldn)(itile,:) = Bstats.resid + dr;
    catch ME
      fprintf(1, '%s\n', ME.message);
      %continue;
    end
end
toc
warning('on','all')


end   % end ndims == 2








% do WV fraction
%  ydat = squeeze(wv_fc(indt,ilay,itile)) - 1; %  - wv_mn(ilay,itile);

% ----------------------------------------------
% pick level to check 800 mb (lev 90)
junk  = squeeze(fits.skt_coef(:,2));
junk  = squeeze(fits.wv_coef(90,:,2));
junk  = squeeze(fits.rh_coef(90,:,2));
junk  = squeeze(fits.tz_coef(90,:,2));
junk  = squeeze(fits.o3pp_coef(20,:,2));
junk  = squeeze(o3pp(30,:,1137));
junk2 = reshape(junk,[64 72]);
junk2 = reshape(junk,[98, 64, 72]);
jnk_zmn = nanmean(junk2,3);


fh = pcolor(airs.lonB(1:72), airs.latB(1:64), junk2);
 set(fh, 'EdgeColor', 'none');colorbar;
 shading 'interp'


% MAP trend
addpath /asl/matlib/maps          % aslmap
addpath /home/strow/Matlab/Extra/
addpath /asl/matlib/plotutils
load llsmap5

mopts.color = 'k';
mopts.title = 'ERA.I RH 800mb trend /yr';
mopts.caxis = [-Inf Inf];
%mopts.caxis = [-0.25 0.25];
mopts.cmap  = llsmap5;
mopts.titlesize = 14;


fh = aslmap(fign, rot90(airs.latB(1:65)), airs.lonB(1:73), junk2, [-90 90], [-180 180], mopts);
>>


% Get zonal average of original gas
o3_zmn = [];
for ilev = 1:101
  junk  = squeeze(o3(30,ilev,:));
  junk2 = reshape(junk,[64 72]);
  for ilat = 1:64
     o3_zmn(ilat,ilev) = nanmean(junk2(ilat,:),2);
  end  
end
fh=pcolor([1:64],[101:-1:1],o3_zmn');

% Get zonal average of linear rate
kn = 6;
g_rte_zmn = [];
[np,nlat,nlon] = size(fits.(cfldn))
junk  = squeeze(fits.(cfldn)(:,:,2));
junk2 = reshape(junk,[np, 64, 72]);

for ilev = 1:np
  g_rte_zmn(ilev,:) = nanmean(junk2(ilev,:,:),3);
end

junk  = squeeze(fits.tz_coef(:,:,2));
junk2 = reshape(junk,[99,64 72]);
tz_zmn_rte = nanmean(junk2,3);
 
fh=pcolor(airs.latB, plev(1:np), g_rte_zmn);
  set(fh, 'EdgeColor', 'none');colorbar;
  shading 'interp'  
  set(gca,'YDir','reverse'); set(gca,'YScale','log')  
  ylim([1 1000])
  xlabel('Latitude');ylabel('Level (hPa)');
  title('ERA.I O3 ppm lin rate /yr')
      
  
