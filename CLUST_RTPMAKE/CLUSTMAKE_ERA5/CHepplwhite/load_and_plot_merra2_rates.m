


addpath /home/chepplew/myLib/matlib/math

% for plotting and reference load AIRS tile boundaries
load('/home/motteler/shome/airs_tiling/latB64.mat','latB2');  % latB2 x 65
airs.latB = latB2; clear latB2;
airs.lonB = [-180:5:180];

d.home  = '/home/chepplew/data/merra2/';
fpattn1 = 'MERRA2_*.tavg3_3d_asm_Nv.20*01.nc4';
fpattn2 = 'MERRA2_*3d_cld_Np*.nc4';

d(1).dir = dir([d.home fpattn1]);
d(2).dir = dir([d.home fpattn2]);

dnum1 = [];
for ifn = 1:length(d(1).dir)
  junk = strsplit(d(1).dir(ifn).name,{'_','.'});
  dnum1(ifn) = datenum(junk{7},'YYYYmmdd');
end
%
dnum2 = [];
for ifn = 1:length(d(2).dir)
  junk = strsplit(d(2).dir(ifn).name,{'_','.'});
  dnum2(ifn) = datenum(junk{7},'YYYYmm');
end

[mer.dnum1 iyt1] = sort(dnum1);
%
[mer.dnum2 iyt2] = sort(dnum2);
%% plot(mer.dnum1,'.')

tza_sub = [];
k = 1;
for ifn = iyt1
  fname = [d(1).dir(ifn).folder '/' d(1).dir(ifn).name];
  fninfo = ncinfo(fname);
  if(k == 1)
    lon  = ncread(fname,'lon');
    lat  = ncread(fname,'lat');
    lev  = ncread(fname,'lev');      % 1:TOA, 72:SFC
  end
  xtim  = ncread(fname,'time');
  
  o3   = ncread(fname,'O3');
  %rh   = ncread(fname,'rh');

  tz   = ncread(fname,'T');        % lon x lat x lev
  tza  = nanmean(tz,4);

  for i = 1:length(airs.latB)-1
    for j = 1:length(airs.lonB)-1
      iix = find(lat > airs.latB(i) & lat <= airs.latB(i+1));
      iiy = find(lon > airs.lonB(j) & lon <= airs.lonB(j+1));

      tza_sub(k,j,i,:) = squeeze(nanmean(tza(iiy, iix, :),[1,2]));

    end
  end

  k = k+1;
  fprintf(1,'.')
end

% -----------
% RH from CLD
% -----------
rha_sub = [];
k = 1;
for ifn = iyt2
  fname = [d(2).dir(ifn).folder '/' d(2).dir(ifn).name];
  fninfo = ncinfo(fname);
  if(k == 1)
    lon  = ncread(fname,'lon');
    lat  = ncread(fname,'lat');
    lev  = ncread(fname,'lev');      % 1:TOA, 72:SFC
  end
  xtim  = ncread(fname,'time');
  rha   = ncread(fname,'RH');

  for i = 1:length(airs.latB)-1
    for j = 1:length(airs.lonB)-1
      iix = find(lat > airs.latB(i) & lat <= airs.latB(i+1));
      iiy = find(lon > airs.lonB(j) & lon <= airs.lonB(j+1));

      rha_sub(k,j,i,:) = squeeze(nanmean(rha(iiy, iix, :),[1,2]));

    end
  end

  k = k+1;
  fprintf(1,'.')
end



  
% preselect field names
allvars = {'tza_sub','rha_sub'};

kn = 1;
sfldn = allvars{kn};
cfldn = [sfldn '_coef'];
rfldn = [sfldn '_err'];
afldn = [sfldn '_anom'];
yvar  = mer.(sfldn);

% QC
inan = find(yvar <=0 | yvar == NaN);


indt = find(mer.dnum);
% --------------------
% Fit trend - profiles
% --------------------
clear tmp*
if(ndims(yvar) == 4) 
  [nt, nlon, nlat, nlev] = size(yvar);
  xdat   = mer.dnum2(indt) - mer.dnum2(indt(1));
  warning('off','all')
  tic
  parfor ilev = 1:nlev
    for ij = 1:nlat
    for ik = 1:nlon
      ydat = squeeze( yvar(indt, ik, ij, ilev)); % - ...
                     %nanmean(airs.(sfldn)(indt, ilev, ij, ik),1) );
      try
        [Bc, Be, Bstats] = Math_tsfit_lin_robust(xdat,ydat,4);
        tmp1(ik,ij,ilev,:) = Bc;
        %yhat = Bc(2)*(xdat)./365 + Bc(1);
        % Lag-1 uncertainty estimates
        lxc    = xcorr(Bstats.resid,1,'coeff');
        lag1   = lxc(1);
        lag1_c = sqrt( (1 + lag1) ./ ( 1 - lag1));
        tmp2(ik,ij,ilev) = Bstats.se(2).*lag1_c;
        % to compare with SM anom retrievals add back linear term
        dr = (xdat + 1)./365. *Bc(2);
        tmp3(ik,ij,ilev,:) = Bstats.resid + dr';
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







% MAP trend
plat = 0.5*(airs.latB(1:end-1) + airs.latB(2:end));
plon = 0.5*(airs.lonB(1:end-1) + airs.lonB(2:end));


addpath /asl/matlib/maps          % aslmap 
addpath /home/strow/Matlab/Extra/
addpath /asl/matlib/plotutils
load llsmap5

mopts.color = 'k';
mopts.title = 'MERRA2';
mopts.caxis = [-Inf Inf];
%mopts.caxis = [-0.25 0.25];
mopts.cmap  = llsmap5;
mopts.titlesize = 14;

junk = squeeze(tza_sub(:,:,72));
fign = 2;
fh = aslmap(fign, (airs.latB(1:end-1)), airs.lonB(1:end-1), junk, [-90 90], [-180 180], mopts);
fh = aslmap(fign, (airs.latB)', airs.lonB', (junk)', [-90 90], [-180 180], mopts);

plev = 6;
junk = squeeze(fit.(cfldn)(:,:,plev,2));
