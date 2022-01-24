function [h,ha,p,pa] = airsL3_climatology(yIN,mIN,dIN,gIN);

%yIN = 2015;
%mIN = 04;
%dIN = 24;
%gIN = 177;

addpath /asl/matlib/h4tools

addpath /home/sergio/MATLABCODE/JPL_DUST_Nov2014/DUSTFLAG
addpath /home/sergio/MATLABCODE/JPL_DUST_Nov2014/MAIN
addpath /home/sergio/MATLABCODE/JPL_DUST_Nov2014/PLOTTER

for ii = 1 : 4
  figure(ii); colormap(jet);
end

[zdust,hVolc,p,zz35,filename] = is_there_dust_hdf(yIN,mIN,dIN,gIN,'F');
p.rlon = wrapTo180(p.rlon);
figure(1); scatter_coast(p.rlon,p.rlat,50,zz35);          title('BT1231 - BT961');
figure(2); scatter_coast(p.rlon,p.rlat,50,p.dust_flag);   title('DustFLag');
figure(3); scatter_coast(p.rlon,p.rlat,50,p.dust_score);  title('Dust Score')
figure(4); scatter_coast(p.rlon,p.rlat,50,p.BT_diff_SO2); title('BT Diff SO2')

addpath /asl/matlib/aslutil/
if ~exist('Airs_Temp_A')
  disp('loading in climatolgy from /asl/s1/sergio/AIRS_L3/airs_L3v6_*_Sept2014.mat');
  load /asl/s1/sergio/AIRS_L3/airs_L3v6_Sept2014.mat
  load /asl/s1/sergio/AIRS_L3/airs_L3v6_extra_Sept2014.mat
end

Airs_PQ = [1000, 925, 850, 700, 600, 500, 400, 300, 250, 200, 150, 100];
Airs_PT = [1000, 925, 850, 700, 600, 500, 400, 300, 250, 200, 150, 100, ...
           70, 50, 30, 20, 15, 10, 7, 5, 3, 2, 1.5,1];
AIRS_PO3 = Airs_PT;

% Airs_Date_Start = datenum(2007,01,01);
% Airs_Date_End   = datenum(2007,12,31);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[zlon,zlat] = meshgrid(-179.5:+179.5,-89.5:+89.5);

nine_year = 5 : 108+5-1;
%junk = [yy(nine_year(1)) mm(nine_year(1)) yy(nine_year(end)) mm(nine_year(end))];
%fprintf(1,'yy/mm start = %4i/%2i   end = %4i/%2i \n',junk);

ten_year = 5 : 120+5-1;
ten_year = 1 : 120;

iimm = find(mm == mIN);
dn = nanmean(p.solzen);
if dn < 90
  stempXY = squeeze(Airs_STemp_A(iimm,:,:));
  spresXY = squeeze(Airs_SPres_A(iimm,:,:));
  wvXY    = squeeze(Airs_H2OVap_A(iimm,:,:,:));
  o3XY    = squeeze(Airs_Ozone_A(iimm,:,:,:));  
  tXY     = squeeze(Airs_Temp_A(iimm,:,:,:));
else
  stempXY = squeeze(Airs_STemp_D(iimm,:,:));
  spresXY = squeeze(Airs_SPres_D(iimm,:,:));  
  wvXY    = squeeze(Airs_H2OVap_D(iimm,:,:,:));
  o3XY    = squeeze(Airs_Ozone_D(iimm,:,:,:));  
  tXY     = squeeze(Airs_Temp_D(iimm,:,:,:));
end

woo = find(stempXY < 0); stempXY(woo) = NaN;
woo = find(spresXY < 0); spresXY(woo) = NaN;
woo = find(wvXY < 0); wvXY(woo) = NaN;
woo = find(o3XY < 0); o3XY(woo) = NaN;
woo = find(tXY < 0); tXY(woo) = NaN;

stempXY = squeeze(nanmean(stempXY,1));
spresXY = squeeze(nanmean(spresXY,1));
tXY     = squeeze(nanmean(tXY,1));
wvXY    = squeeze(nanmean(wvXY,1));
o3XY    = squeeze(nanmean(o3XY,1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p.rlon = wrapTo180(p.rlon);

rlon = p.rlon; rlon(rlon < min(zlon(:))) = min(zlon(:)); rlon(rlon > max(zlon(:))) = max(zlon(:));
rlat = p.rlat; rlat(rlat < min(zlat(:))) = min(zlat(:)); rlat(rlat > max(zlat(:))) = max(zlat(:));

for iL = 1 : 24
  woo = squeeze(tXY(iL,:,:));
  woo = flipud(woo);
  bad = find(woo < 0); 
  if length(bad) > 0
    good = find(woo > 0);
    %boo = interp2(zlon(good),zlat(good),woo(good),zlon(bad),zlat(bad));
    %woo(bad) = boo;
    woo(bad) = Inf;
  end
  tZ(iL,:) = interp2(zlon,zlat,woo,rlon,rlat,'linear',NaN);
  boo = find(tZ(iL,:) == Inf | tZ(iL,:) == NaN); tZ(iL,boo) = NaN;
  figure(1); scatter_coast(p.rlon,p.rlat,50,tZ(iL,:)); ax = axis; cx = caxis;
    title(['T lev=' num2str(iL) ' p=' num2str(Airs_PT(iL))]);  
  figure(2); simplemap(woo); axis(ax); caxis(cx); title(num2str(Airs_PT(iL)));
    title(['T lev=' num2str(iL) ' p=' num2str(Airs_PT(iL))]);
  pause(0.1);  
end
%disp('ret to continue : '); pause

for iL = 1 : 24
  woo = squeeze(o3XY(iL,:,:));
  woo = flipud(woo);
  bad = find(woo < 0); 
  if length(bad) > 0
    good = find(woo > 0);
    %boo = interp2(zlon(good),zlat(good),woo(good),zlon(bad),zlat(bad));
    %woo(bad) = boo;
    woo(bad) = Inf;
  end
  o3Z(iL,:) = interp2(zlon,zlat,woo,rlon,rlat,'linear',NaN);
  boo = find(o3Z(iL,:) == Inf | o3Z(iL,:) == NaN); o3Z(iL,boo) = NaN;
  figure(1); scatter_coast(p.rlon,p.rlat,50,o3Z(iL,:)); ax = axis; cx = caxis;
    title(['O3 lev=' num2str(iL) ' p=' num2str(Airs_PT(iL))]);  
  figure(2); simplemap(woo); axis(ax); caxis(cx); title(num2str(Airs_PT(iL)));
    title(['O3 lev=' num2str(iL) ' p=' num2str(Airs_PT(iL))]);
  pause(0.1);
end
%disp('ret to continue : '); pause

for iL = 1 : 12
  woo = squeeze(wvXY(iL,:,:));
  woo = flipud(woo);
  bad = find(woo < 0); 
  if length(bad) > 0
    good = find(woo > 0);
    %boo = interp2(zlon(good),zlat(good),woo(good),zlon(bad),zlat(bad));
    %woo(bad) = boo;
    woo(bad) = NaN;
    woo(bad) = Inf;
  end
  %simplemap(woo); %caxis([220 320]); colorbar; colormap(jet)
  wvZ(iL,:) = interp2(zlon,zlat,woo,rlon,rlat,'linear',NaN);
  boo = find(wvZ(iL,:) == Inf | wvZ(iL,:) == NaN); wvZ(iL,boo) = NaN;
  figure(1); scatter_coast(p.rlon,p.rlat,50,wvZ(iL,:)); ax = axis; cx = caxis;
    title(['WV lev=' num2str(iL) ' p=' num2str(Airs_PQ(iL))]);
  figure(2); simplemap(woo); axis(ax); caxis(cx);
    title(['WV lev=' num2str(iL) ' p=' num2str(Airs_PQ(iL))]);  
  pause(0.1);
end

woo = flipud(stempXY);
bad = find(woo < 0); 
if length(bad) > 0
  good = find(woo > 0);
  woo(bad) = Inf;
end
stemp = interp2(zlon,zlat,woo,rlon,rlat); 
boo = find(stemp == Inf | stemp == NaN); stemp(boo) = NaN;
  figure(1); scatter_coast(p.rlon,p.rlat,50,stemp); ax = axis; cx = caxis;
  figure(2); simplemap(woo); axis(ax); caxis(cx);

woo = flipud(spresXY);
bad = find(woo < 0); 
if length(bad) > 0
  good = find(woo > 0);
  woo(bad) = Inf;
end
spres = interp2(zlon,zlat,woo,rlon,rlat); 
boo = find(spres == Inf | spres == NaN); spres(boo) = NaN;
  figure(3); scatter_coast(p.rlon,p.rlat,50,spres); ax = axis; cx = caxis;
  figure(4); simplemap(woo); axis(ax); caxis(cx);
%disp('ret to continue : '); pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

extraWV = ones(length(Airs_PT)-length(Airs_PQ),length(p.rlat)) * -9999;
extraWV = ones(length(Airs_PT)-length(Airs_PQ),1) * wvZ(length(Airs_PQ),:);
new_wvZ = [wvZ; extraWV];

p.stemp = stemp;
p.spres = spres;
p.plevs = flipud(Airs_PT' * ones(1,length(p.rlat)));
p.ptemp = flipud(tZ);
p.gas_1 = flipud(new_wvZ);
p.gas_3 = flipud(o3Z);
Airs_PTx = fliplr(Airs_PT);
for ii = 1 : length(p.rlat)
  dat = p.ptemp(:,ii);
  woo = find(Airs_PTx < p.spres(ii) & isfinite(dat'));
  woo = length(woo);
  orig_nlevs(ii) = woo;
  orig_tgnd(ii)  = dat(woo);
  final_tgnd(ii) = dat(woo);  
  if woo < 24 & dat(woo+1) > 0
    woo = woo + 1;
    final_tgnd(ii) = dat(woo);      
  end
  p.nlevs(ii) = woo;
end
p.plat = p.rlat;
p.plon = p.rlon;
p.wspeed = zeros(size(p.rlat));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%pa = {{'profiles','rtime','seconds since 1993'}};
pa = {{'profiles','rtime','seconds since 1958'}};
ha = {{'header','hdf file',filename}};

addpath /asl/packages/time
addpath /home/sbuczko1/git/rtp_prod2/emis/
addpath /home/sbuczko1/git/rtp_prod2/util
[p,pa] = rtp_add_emis(p,pa);

h.ptype   = 0;  %% levels
h.pfields = 5;  %% 1 (prof) + 4 (obs)
h.nchan = 2378;
h.ichan = (1:2378)';
h.ngas = 2;
h.gasid = [1 3]';
h.gunit = [20 12]';

klayers = '/asl/packages/klayers/Bin/klayers_airs';
sarta   = '/asl/packages/sartaV108/Bin/sarta_apr08_m140_wcon_nte';
fip = mktemp('junk.ip.rtp');
fop = mktemp('junk.op.rtp');
frp = mktemp('junk.rp.rtp');

rtpwrite(fip,h,ha,p,pa);

disp('klayers ....')
klayerser = ['!' klayers ' fin=' fip ' fout=' fop ' >& ugh']; eval(klayerser)
disp('sarta ....')
sartaer   = ['!' sarta '   fin=' fop ' fout=' frp]; eval(sartaer)

[hnew,hanew,pnew,panew] = rtpread(frp);
p.rcalc = pnew.rcalc;
h.pfields = 7;  %% 1 (prof) + 4 (obs) + 2 (calc)

if isfield(hnew,'vchan') & ~isfield(h,'vchan')
  h.vchan = hnew.vchan;
end

rmer = ['!/bin/rm ' fip ' ' fop ' ' frp];
eval(rmer);

disp('done!!!!!')
