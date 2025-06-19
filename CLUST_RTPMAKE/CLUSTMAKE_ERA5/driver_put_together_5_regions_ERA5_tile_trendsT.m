addpath /home/sergio/MATLABCODE/COLORMAP/LLS

%% this cannot do ALL profiles so just do a set at a time ef 0001-1536, 1537-3072, 3073-4608
load latB64.mat
rlat65 = latB2; rlon73 = -180 : 5 : +180;
rlon = -180 : 5 : +180;  rlat = latB2;
rlon = 0.5*(rlon(1:end-1)+rlon(2:end));
rlat = 0.5*(rlat(1:end-1)+rlat(2:end));

[Y,X] = meshgrid(rlat,rlon);
Y0 = Y;
X0 = X;
Y = Y(:); X = X(:);

%%%%%%%%%%%%%%%%%%%%%%%%%

plays = load('/home/sergio/MATLABCODE/airslevels.dat');
plays = [1100*1.001; plays];
plays = flipud(meanvaluebin(plays));

load llsmap5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tm2 = load('trends_summary_10percent_from_ERA5clearcalc_T_region-2.mat');
tm1 = load('trends_summary_10percent_from_ERA5clearcalc_T_region-1.mat');
t00 = load('trends_summary_10percent_from_ERA5clearcalc_T_region00.mat');
tp1 = load('trends_summary_10percent_from_ERA5clearcalc_T_region01.mat');
tp2 = load('trends_summary_10percent_from_ERA5clearcalc_T_region02.mat');

%put_together_skt_trends

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tz1  = [tm2.ptemp_trend_Q90_cld; tm1.ptemp_trend_Q90_cld; t00.ptemp_trend_Q90_cld; tp1.ptemp_trend_Q90_cld; tp2.ptemp_trend_Q90_cld]';
  figure(1); clf; tzz1 = reshape(tz1,101,72,64); tzz1 = squeeze(nanmean(tzz1,2)); pcolor(rlat,plays,tzz1); title('Q90 BT1231 cld dT/dt [K/yr]')
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([10 1000]); colorbar; caxis([-1 +1]*0.15); shading interp; xlabel('Lat'); ylabel('Pressure [mb]');
tz2  = [tm2.ptemp_trend_Q90_clr; tm1.ptemp_trend_Q90_clr; t00.ptemp_trend_Q90_clr; tp1.ptemp_trend_Q90_clr; tp2.ptemp_trend_Q90_clr]';
  figure(2); clf; tzz2 = reshape(tz2,101,72,64); tzz2 = squeeze(nanmean(tzz2,2)); pcolor(rlat,plays,tzz2); title('Q90 BT1231 clr dT/dt [K/yr]')
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([10 1000]); colorbar; caxis([-1 +1]*0.15); shading interp; xlabel('Lat'); ylabel('Pressure [mb]');
tz3  = [tm2.ptemp_trend_Q90    ; tm1.ptemp_trend_Q90    ; t00.ptemp_trend_Q90    ; tp1.ptemp_trend_Q90    ; tp2.ptemp_trend_Q90    ]';
  figure(3); clf; tzz3 = reshape(tz3,101,72,64); tzz3 = squeeze(nanmean(tzz3,2)); pcolor(rlat,plays,tzz3); title('Q90 dT/dt [K/yr]')
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([10 1000]); colorbar; caxis([-1 +1]*0.15); shading interp; xlabel('Lat'); ylabel('Pressure [mb]');
tz4  = [tm2.ptemp_trend_mean   ; tm1.ptemp_trend_mean   ; t00.ptemp_trend_mean   ; tp1.ptemp_trend_mean   ; tp2.ptemp_trend_mean   ]';
  figure(4); clf; tzz4 = reshape(tz4,101,72,64); tzz4 = squeeze(nanmean(tzz4,2)); pcolor(rlat,plays,tzz4); title('mean dT/dt [K/yr]')
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([10 1000]); colorbar; caxis([-1 +1]*0.15); shading interp; xlabel('Lat'); ylabel('Pressure [mb]');
tz5  = [tm2.ptemp_trend_cntr   ; tm1.ptemp_trend_cntr   ; t00.ptemp_trend_cntr   ; tp1.ptemp_trend_cntr   ; tp2.ptemp_trend_cntr   ]';
  figure(5); clf; tzz5 = reshape(tz5,101,72,64); tzz5 = squeeze(nanmean(tzz5,2)); pcolor(rlat,plays,tzz5); title('cntr dT/dt [K/yr]')
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([10 1000]); colorbar; caxis([-1 +1]*0.15); shading interp; xlabel('Lat'); ylabel('Pressure [mb]');
for ii = 1 : 5; figure(ii); set(gca,'fontsize',10); end
disp('ret to continue'); pause;

figure(1); clf; tzz1 = reshape(tz1-0*tz1,101,72,64); tzz1 = squeeze(nanmean(tzz1,2)); pcolor(rlat,plays,tzz1); title('Q90 BT1231 cld dT/dt [K/yr]')
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([10 1000]); colorbar; caxis([-1 +1]*0.15); shading interp; xlabel('Lat'); ylabel('Pressure [mb]');
figure(2); clf; tzz2 = reshape(tz1-tz2,101,72,64); tzz2 = squeeze(nanmean(tzz2,2)); pcolor(rlat,plays,tzz2); title('Abs diff Q90 BT1231 clr dT/dt [K/yr]')
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([10 1000]); colorbar; caxis([-1 +1]*0.15/10); shading interp; xlabel('Lat'); ylabel('Pressure [mb]');
figure(3); clf; tzz3 = reshape(tz1-tz3,101,72,64); tzz3 = squeeze(nanmean(tzz3,2)); pcolor(rlat,plays,tzz3); title('Abs diff Q90 dT/dt [K/yr]')
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([10 1000]); colorbar; caxis([-1 +1]*0.15/10); shading interp; xlabel('Lat'); ylabel('Pressure [mb]');
figure(4); clf; tzz4 = reshape(tz1-tz4,101,72,64); tzz4 = squeeze(nanmean(tzz4,2)); pcolor(rlat,plays,tzz4); title('Abs diff mean dT/dt [K/yr]')
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([10 1000]); colorbar; caxis([-1 +1]*0.15/10); shading interp; xlabel('Lat'); ylabel('Pressure [mb]');
figure(5); clf; tzz5 = reshape(tz1-tz5,101,72,64); tzz5 = squeeze(nanmean(tzz5,2)); pcolor(rlat,plays,tzz5); title('Abs diff cntr dT/dt [K/yr]')
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([10 1000]); colorbar; caxis([-1 +1]*0.15/10); shading interp; xlabel('Lat'); ylabel('Pressure [mb]');
for ii = 1 : 5; figure(ii); set(gca,'fontsize',10); end
disp('ret to continue'); pause;

figure(1); clf; tzz1 = reshape(tz1-0*tz1,101,72,64); tzz1 = squeeze(nanmean(tzz1,2)); pcolor(rlat,plays,tzz1); title('Q90 BT1231 cld dT/dt [K/yr]')
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([10 1000]); colorbar; caxis([-1 +1]*0.15); shading interp; xlabel('Lat'); ylabel('Pressure [mb]');
figure(2); clf; tzz2 = reshape(tz1-tz2,101,72,64); tzz2 = squeeze(nanmean(tzz2,2)); tzz2 = tzz2./tzz1 * 100; pcolor(rlat,plays,tzz2); title('Percent diff Q90 BT1231 clr dT/dt [K/yr]')
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([10 1000]); colorbar; caxis([-1 +1]*10); shading interp; xlabel('Lat'); ylabel('Pressure [mb]');
figure(3); clf; tzz3 = reshape(tz1-tz3,101,72,64); tzz3 = squeeze(nanmean(tzz3,2)); tzz3 = tzz3./tzz1 * 100; pcolor(rlat,plays,tzz3); title('Percent diff Q90 dT/dt [K/yr]')
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([10 1000]); colorbar; caxis([-1 +1]*10); shading interp; xlabel('Lat'); ylabel('Pressure [mb]');
figure(4); clf; tzz4 = reshape(tz1-tz4,101,72,64); tzz4 = squeeze(nanmean(tzz4,2)); tzz4 = tzz4./tzz1 * 100; pcolor(rlat,plays,tzz4); title('Percent diff mean dT/dt [K/yr]')
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([10 1000]); colorbar; caxis([-1 +1]*10); shading interp; xlabel('Lat'); ylabel('Pressure [mb]');
figure(5); clf; tzz5 = reshape(tz1-tz5,101,72,64); tzz5 = squeeze(nanmean(tzz5,2)); tzz5 = tzz5./tzz1 * 100; pcolor(rlat,plays,tzz5); title('Percent diff cntr dT/dt [K/yr]')
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([10 1000]); colorbar; caxis([-1 +1]*10); shading interp; xlabel('Lat'); ylabel('Pressure [mb]');
for ii = 1 : 5; figure(ii); set(gca,'fontsize',10); end
