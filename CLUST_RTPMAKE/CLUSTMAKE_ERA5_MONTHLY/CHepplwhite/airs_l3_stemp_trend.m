% AIRS L3 trend maps

%

addpath /home/strow/Matlab/Math
addpath /home/chepplew/myLib/matlib/math

% Get date range if requested
switch nargin
  case 0
    % Define date range to match the retrieval dataset.
    sdate = '2003/01/01'; % '2002/09/08';
    edate = '2020/09/30'; % '2020/10/29';
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
  otherwise
    error('invalid number of arguments [0 or 2]');
end

airs.test_tiles = [27, 47;    % seasonal
                   54, 24;    % clear
		   67, 35;    % cloudy (DCC)
		   64, 34;    % Equator (lonbin, latbin) (N,:) N=1,2,3
		   36,  3;    % Antarctic
		   36, 62];


% output structure with results
fit    = struct;

d.home = '/asl/airs/AIRS3STM/v7/';
d.adir = dir([d.home '*/AIRS.*.L3.*.hdf']);
%d.ddir = dir([d.home '/' 'AIRS_dbase_2*_descending.mat']);

disp('Initializing')
% AIRS.L3 1-deg tile boundaries (  ):
l3.lonB = [-180:1:180.00];
l3.latB = [90:-1:-90];

% ASL AIRS.tile boundaries: (64 x 72 lat x lon tiles)
load('/home/motteler/shome/airs_tiling/latB64.mat','latB2');  % latB2 x 65
airs.latB = latB2; clear latB2;
airs.lonB = [-180:5:180];

% Specified AIRS.tile region to load:
iwb = 3;
lonBin = airs.test_tiles(iwb,1);
latBin = airs.test_tiles(iwb,2);
disp(['lonBin: ' num2str(lonBin) ' latBin: ' num2str(latBin)]);

% get L3 tiles in requested airs.tile
inda = find(l3.latB >= airs.latB(latBin) & l3.latB <= airs.latB(latBin+1) );
indn = find(l3.lonB >= airs.lonB(lonBin) & l3.lonB <= airs.lonB(lonBin+1) );

% assign output variables
airs.dnum   = [];   airs.nod    = [];

airs.skt_a  = [];   airs.skt_d = [];
airs.sct_a  = [];   airs.sct_d = [];
airs.tz_a   = [];  airs.wvmr_a = [];
airs.o3mr_a = [];    airs.rh_a = [];
airs.tz_d   = [];  airs.wvmr_d = [];
airs.o3mr_d = [];    airs.rh_d = [];

fname = [d.adir(1).folder '/' d.adir(1).name];
finfo = hdfinfo(fname);
for i=1:length(finfo.Vgroup(1).Vgroup(1).SDS)
  varnams{i} = finfo.Vgroup(1).Vgroup(1).SDS(i).Name; 
end

k = 1;
for i=1:length(d.adir)
  fname = [d.adir(i).folder '/' d.adir(i).name];
  
  if(i==1);
    airs.lat  = hdfread(fname,'Latitude');
    airs.lon  = hdfread(fname,'Longitude');
    airs.plev = cell2mat(hdfread(fname,'StdPressureLev'));
    airs.hlev = cell2mat(hdfread(fname,'H2OPressureLay'));
  end
  %airs.skt  = [airs.skt hdfread(fname,'SurfSkinTemp_A')];
  airs.skt_a(k,:,:)  = flipud(hdfread(fname,'SurfSkinTemp_A'));
  airs.sct_a(k,:,:)  = flipud(hdfread(fname,'SurfSkinTemp_A_ct'));
  airs.skt_d(k,:,:)  = flipud(hdfread(fname,'SurfSkinTemp_D'));
  airs.sct_d(k,:,:)  = flipud(hdfread(fname,'SurfSkinTemp_D_ct'));
       xyr  = hdfread(fname,'Year');
       xmon = hdfread(fname,'Month');
       xday = hdfread(fname,'Day');
       junk = [double(cell2mat(xyr)) double(cell2mat(xmon)) double(cell2mat(xday))];
  airs.dnum = [airs.dnum; datenum(junk)];
  airs.nod  = [airs.nod; hdfread(fname,'NumOfDays')];

  airs.tz_a(k,:,:,:)    = flipud(hdfread(fname, 'Temperature_A'));
  airs.wvmr_a(k,:,:,:)  = flipud(hdfread(fname, 'H2O_MMR_A'));
  airs.rh_a(k,:,:,:)    = flipud(hdfread(fname, 'RelHum_A'));
  airs.o3mr_a(k,:,:,:)  = flipud(hdfread(fname, 'O3_VMR_A'));

  airs.tz_d(k,:,:,:)    = flipud(hdfread(fname, 'Temperature_D'));
  airs.wvmr_d(k,:,:,:)  = flipud(hdfread(fname, 'H2O_MMR_D'));
  airs.rh_d(k,:,:,:)    = flipud(hdfread(fname, 'RelHum_D'));
  airs.o3mr_d(k,:,:,:)  = flipud(hdfread(fname, 'O3_VMR_D'));
  airs.ch4vr_d(k,:,:,:) = flipud(hdfread(fname, 'CH4_VMR_D'));

  % convert o3 mr to ppmv
  
  %yxxx = toppmv(plevs,txxx, o3xxx, 48,12); 
  airs.o3ppm_d(k,:,:,:) = toppmv(airs.plev,airs.tz_d(k,end:-1:1,:,:), ...
                          airs.o3mr_d(k,end:-1:1,:,:), 48,12); 

  k = k+1;  
  fprintf(1,'.')
end

% --------------------------------------
% Set up fieldnames for group processing
allvars = {'skt_a','skt_d','wvmr_a','tz_d','h2o_d','rh_d','o3mr_d',...
           'o3ppm_d','ch4_d'};
kn = 8;
sfldn = allvars{kn};
cfldn = [allvars{kn} '_coef'];
rfldn = [allvars{kn} '_rerr'];
afldn = [allvars{kn} '_anom'];

% Screen for NaNs (-9999)
inan = find(airs.(sfldn) == -9999 | airs.(sfldn) <= 0);
whos inan
airs.(sfldn)(inan) = NaN;


%{
% -----------------------------
% Group into airs 64 x 72 tiles
iiw=cell(63,71);
ncell = zeros(63,71);
for ij = 1:63
  for ik = 1:71
     iiw{ij,ik} = find(airs.lat >= airs.latB2(ij) & airs.lat <= airs.latB2(ij+1) ...
            &   airs.lon >= airs.lonB2(ik) & airs.lon <= airs.lonB2(ik+1));
	    
     ncell(ij,ik) = length(iiw{ij,ik});

  end
  fprintf(1,'.')
end
%}

% Subset time period if requested:
switch nargin
  case 0
    indt  = find(airs.dnum);
  case 2
    indt  = find(airs.dnum >= D1 & airs.dnum <= D2);
end

% --------------------
% Fit trend - profiles
% --------------------
clear tmp*
if(ndims(airs.(sfldn)) == 4) 
  [nt, nlev, nlat, nlon] = size(airs.(sfldn));
  xdat   = airs.dnum(indt) - airs.dnum(indt(1));
  warning('off','all')
  tic
  parfor ilev = 1:nlev
    for ij = 1:nlat
    for ik = 1:nlon
      ydat = squeeze( airs.(sfldn)(indt, ilev, ij, ik)); % - ...
                     %nanmean(airs.(sfldn)(indt, ilev, ij, ik),1) );
      try
        [Bc, Be, Bstats] = Math_tsfit_lin_robust(xdat,ydat,4);
        tmp1(ilev,ij,ik,:) = Bc;
        %yhat = Bc(2)*(xdat)./365 + Bc(1);
        % Lag-1 uncertainty estimates
        lxc    = xcorr(Bstats.resid,1,'coeff');
        lag1   = lxc(1);
        lag1_c = sqrt( (1 + lag1) ./ ( 1 - lag1));
        tmp2(ilev,ij,ik) = Bstats.se(2).*lag1_c;
        % to compare with SM anom retrievals add back linear term
        dr = (xdat + 1)./365. *Bc(2);
        tmp3(ilev,ij,ik,:) = Bstats.resid + dr;
       catch ME
        fprintf(1, '%s\n', ME.message);
        continue; 
      end

     end  % end ilon
     fprintf(1,'.')
  end  % end ilat
  end  % end ilev
  toc

end   % end if ndims == 4.

warning('on','all')

fit.(cfldn) = tmp1;
fit.(rfldn) = tmp2;
fit.(afldn) = tmp3;
clear tmp*

% --------------------
% Fit trend - surface
% --------------------
if(ndims(airs.(sfldn)) == 3) 
  [nt, nlat, nlon] = size(airs.(sfldn));
  xdat   = airs.dnum(indt) - airs.dnum(indt(1));
  warning('off','all')
  tic
  parfor ij = 1:nlat
    for ik = 1:nlon
      ydat = squeeze( airs.(sfldn)(indt, ij, ik) - ...
                     nanmean(airs.(sfldn)(indt, ij, ik),1) );
      try
        [Bc, Be, Bstats] = Math_tsfit_lin_robust(xdat,ydat,4);
        tmp1(ij,ik,:) = Bc;
        %yhat = Bc(2)*(xdat)./365 + Bc(1);
        % Lag-1 uncertainty estimates
        lxc    = xcorr(Bstats.resid,1,'coeff');
        lag1   = lxc(1);
        lag1_c = sqrt( (1 + lag1) ./ ( 1 - lag1));
        tmp2(ij,ik) = Bstats.se(2).*lag1_c;
        % to compare with SM anom retrievals add back linear term
        dr = (xdat + 1)./365. *Bc(2);
        tmp3(ij,ik,:) = Bstats.resid + dr;
       catch ME
        fprintf(1, '%s\n', ME.message);
        continue; 
      end

     end  % end ilon
     fprintf(1,'.')
  end  % end ilat
  
end   % end if ndims == 3.
toc
warning('on','all')

fit.(cfldn) = tmp1;
fit.(rfldn) = tmp2;
fit.(afldn) = tmp3;
clear tmp*



% ----------------------------------------
% Interpoloate to AIRS.tiles
% convert from boundaries to centers
tclat = (airs.latB(2:end) + airs.latB(1:end-1) ) / 2;
tclon = (airs.lonB(2:end) + airs.lonB(1:end-1) ) / 2;

% grid of tile center interpolation points
[aX,aY] = meshgrid(tclon, tclat);

% L3 geo are centers, to meshgrid.
[sX,sY] = meshgrid(airs.lon(1,:), airs.lat(:,1));

% Interpolate ERA to airs.tiles and add fields to ouput structure 'infit'.
% just do linear rate for now
% 2:skt_d, 4:tz_d, 

for kn=[2 4]
  clear tdata idata
  cfldn = [allvars{kn} '_coef'];
  if(ndims(fit.(cfldn)) == 3)
    tdata = squeeze(fit.(cfldn)(:,:,2));
    idata = interp2(sX, sY, tdata, aX, aY, 'linear');
  end
  if(ndims(fit.(cfldn)) == 4)
    [np, nlat, nlon, nc] = size(fit.(cfldn));
    for i=1:np
      tdata = squeeze(fit.(cfldn)(i,:,:,2));
      idata(i,:,:) = interp2(sX, sY, tdata, aX, aY, 'linear');
    end
  end
  infit.(cfldn) = idata;
  
end

% -----------------------------------------------------------

% -------------------------------------
% Zonal average vs height
% -------------------------------------

tz_d_zmn    = nanmean(airs.tz_d,4);
tz_d_rte    = nanmean(squeeze(fit.tz_d_coef(:,:,:,2)),3);
o3_d_zmn    = nanmean(airs.o3mr_d,4);
o3_d_rte    = nanmean(squeeze(fit.o3mr_d_coef(:,:,:,2)),3);
o3ppm_d_rte = nanmean(squeeze(fit.o3ppm_d_coef(:,:,:,2)),3);

junk     = squeeze(tz_d_zmn(10,[24:-1:1],:));
junk     = (tz_d_rte([24:-1:1],:));

fh=pcolor(airs.lat(:,1), plev, junk );
  set(fh, 'EdgeColor', 'none'); shading interp;
  colorbar
  set(gca,'YDir','reverse')
  set(gca,'Yscale','log');
  ylabel('height (hPa)');xlabel('latitude');
  
  
%{

% sanity mapping checks
airs.latB2(latBin:latBin+1)
airs.lat(inda,1)
airs.lonB2(lonBin:lonBin+1)
airs.lon(1,indn)

addpath /asl/matlib/aslutil
figure; simplemap(airs.skt_a(:,:,1)); caxis([220 320])
xlat = rot90(repmat([89.5:-1:-89.5],360,1));
xlon = repmat([-179.5:1:179.5],180,1);
figure;simplemap(xlat, xlon, (squeeze( airs.skt_a(:,:,1))) );caxis([220 320]);

fh = pcolor(squeeze(airs.skt_d(:,:,1)) );  set(hf, 'EdgeColor', 'none'); colorbar;


size(fit.(cfldn))
junk = squeeze(fit.(cfldn)(3,:,:,2));    % !! beware levels ranges !!
junk2 = [];


fh = pcolor(airs.skt(:,:,1));
 set(fh, 'EdgeColor', 'none')
 shading interp
 colorbar
 caxis([270 320])

figure;plot(linfit,'.')

% MAP trend
addpath /asl/matlib/maps          % aslmap 
addpath /home/strow/Matlab/Extra/
addpath /asl/matlib/plotutils
load llsmap5

mopts.color = 'k';
mopts.title = 'AIRS L3 v7 RH.D 885mb rate 2003:20  %/yr';
mopts.caxis = [-Inf Inf];
%mopts.caxis = [-0.25 0.25];
mopts.cmap  = llsmap5;
mopts.titlesize = 14;

fign = 3;
alon=[-180:1:180]; alat = [-90:1:90];
 %fh = equal_area_map(fign, gis.latB2, gis.lonB2, SZ, txtstr);
 %fh = aslmap(fign, xlat, xlon, linr, [-90 90], [-180 180], mopts);
fh = aslmap(fign, (l3.latB), l3.lonB, junk, [-90 90], [-180 180], mopts);


% ANOMALY time series
figure; hold on;
for i=1:length(inda);
  for j = 1:length(indn)
    plot(datetime(airs.dnum(indt),'convertfrom','datenum'), ...
         squeeze(anom(inda(i),indn(j),:)),'.-');
  end
end
  title(['AIRS.L3 V7 SKT.D all tiles anom bin: ' ...
   sprintf('%2d,%2d',lonBin,latBin)] )
  ylabel('temp (K)');grid on;

junk = squeeze(nanmean(anom(inda,indn,:),[1 2]));
figure;plot(datetime(airs.dnum(indt),'convertfrom','datenum'),junk,'.-')
  grid on; ylabel('temp (K)'); 
  title(['AIRS.L3 skt.d mean anom. bin: ' ...
         sprintf('%2d,%2d',lonBin,latBin)] );
%}
