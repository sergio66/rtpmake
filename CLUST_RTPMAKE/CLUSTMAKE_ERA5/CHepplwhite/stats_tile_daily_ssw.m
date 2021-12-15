function [] = stats_tiled_daily_ssw()

% Statistical analysis of daily averaged tiled radiances for SSW

addpath /asl/matlib/aslutil
addpath /home/strow/Matlab/Math

% choose date window
xt1 = datetime(2003,01,01);
xt2 = datetime(2014,12,31);

xd1 = datenum(xt1);
xd2 = datenum(xt2);

% subset to date range:
iid = find(rdays);
iid = find(xd1 < rdays & rdays < xd2);

% Convert to BT
clear bto;     % [nlon x nchan x ndays]
for lz=1:72 
  bto(lz,:,:) = rad2bt(fairs(xchs), squeeze(robs(lz,iid,:))'); 
end
btm = rad2bt(fairs(xchs), rmn');

% plot(datetime(rdays(iid),'convertfrom','datenum'), squeeze(bto(lz,ich,:)),'.-')
% hold on; plot(datetime(rdays(iid),'convertfrom','datenum'), btm(ich,:),'-')
 
clear iie;
for lz = 1:72 
  iie(lz) = find(bto(lz,j,:) > 230,1);
end



%{
addpath /home/strow/Matlab/Math                 % Math_tsfit_lin_robust.m

% choose tile lonitude and channel
lz = 1; ich = 7;

% May need to check for NaN
inan=find(isnan(bto(lz,ich,:)));
if(~isempty(inan))
  rdays(inan) = [];
  bto(:,:,inan) = [];
  btm(:,inan) = [];
end
whos rdays bt*
 
X     = rdays-rdays(1);
Y     = squeeze(bto(lz,ich,:))';
bt_mn = nanmean(squeeze(bto(lz,ich,:)))';
Ysp = squeeze(bto(lz,ich,:))-bt_mn;
[Bb, Bstats] = Math_tsfit_lin_robust(X,Ysp,2);
YS_bar = Math_timeseries_2(X,Bb);

% Empirical mode decomposition and hilbert spectrum
[imf,resid,info] = emd(Ysp,'Interpolation','pchip','MaxNumIMF',10);

% Polar band mean
X      = rdays-rdays(1);
btm_mn = nanmean(btm,2);
Ysp    = btm - btm_mn;
[Bb, Bstats] = Math_tsfit_lin_robust(X,Ysp(ich,:),2);
YS_bar = Math_timeseries_2(X,Bb);
[imf,resid,info] = emd(Ysp(ich,:),'Interpolation','pchip','MaxNumIMF',10);
[vmf,vesid]      = vmd(Ysp(ich,:),'NumIMF',10);

Yrec = sum(imf(:,7:10),2);

hold on; 
 plot(datetime(rdays(iid),'convertfrom','datenum'), YS_bar+btm_mn(ich,:),'-')


%}
