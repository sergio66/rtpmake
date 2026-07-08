%
% first: [geo, o3, prodn] =  load_all_mls_lev3_o3(ptype)
%
%
%
%

addpath /home/chepplew/myLib/matlib/math

ptype = 'PZMN'
if(ndims(o3) == 3)
  [nlat, nlev, ntim] = size(o3);
end

sfldn = 'o3pp'
cfldn = [sfldn '_coef'];
rfldn = [sfldn '_rerr'];
afldn = [sfldn '_anom'];


% QC
yvar  = o3;
inan = find(yvar <= 0);
whos inan
yvar(inan) = NaN;

dnum  = datenum(geo.tim);
indt  = find(dnum);
xdat  = dnum(indt); %  - dnum(indt(1)) + 1;

% -------------------------
% Fit Profile vars
% -------------------------

clear tmp*
warning('off','all')
tic
for ilev = 1:nlev
  for ilat = 1:nlat

  ydat = squeeze(yvar(ilat, ilev, indt)); %  - nanmean(yvar(ilay,itile);
    try
      [Bc, Be, Bstats] = Math_tsfit_lin_robust(xdat-xdat(1)+1,ydat,4);
      tmp1(ilat,ilev,:) = Bc;
      %yhat = B(2)*(days(ndt1:ndt2)-days(ndt1))./365 + B(1);
      %linr(i, j,:) = B(2);
% Lag-1 uncertainty estimates
      lxc    = xcorr(Bstats.resid,1,'coeff');
      lag1   = lxc(1);
      lag1_c = sqrt( (1 + lag1) ./ ( 1 - lag1));
      tmp2(ilat,ilev)   = Bstats.se(2).*lag1_c;
% to compare with SM anom retrievals add back linear term
      dr = (xdat - xdat(1)+1)./365. *Bc(2);
      tmp3(ilat,ilev,:) = Bstats.resid + dr';
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


%{
fh=pcolor(geo.lat, geo.lev, squeeze(o3(:,:,10))' )
 set(fh, 'EdgeColor', 'none');colorbar;
 shading 'interp'
 set(gca,'YDir','reverse'); set(gca,'YScale','log');
 ylim([0.05 200]);
 
junk = squeeze(fits.(cfldn)(:,:,2));
xlat = geo.lat(2:end-1);
xlev = geo.lev(5:end-6);
fh = pcolor(xlat,xlev, junk)





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


fh = aslmap(fign, rot90(geo.lat, geo.lon, junk2, [-90 90], [-180 180], mopts);
%}
