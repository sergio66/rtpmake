addpath /home/sergio/KCARTA/MATLAB
addpath /home/sergio/MATLABCODE/
addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS

%{
[hjunk,hajunk,pjunk,pajunk] = rtpread('/asl/rtp/airs_l1c_v6/allfov/2020/003/allfov_ecmwf_airicrad_d2020003_061.rtp'); %% check these are 91 levels
junker = ['!/asl/packages/klayersV205/BinV201/klayers_airs fin=/asl/rtp/airs_l1c_v6/allfov/2020/003/allfov_ecmwf_airicrad_d2020003_061.rtp fout=junkecm.op.rtp >& ugh'];
junker = ['!/asl/packages/klayersV205/BinV201/klayers_airs fin=/asl/rtp/airs_l1c_v6/allfov/2020/003/allfov_ecmwf_airicrad_d2020003_072.rtp fout=junkecm.op.rtp >& ugh'];
junker = ['!/asl/packages/klayersV205/BinV201/klayers_airs fin=/asl/rtp/airs_l1c_v6/allfov/2020/003/allfov_ecmwf_airicrad_d2020003_060.rtp fout=junkecm.op.rtp >& ugh'];
junker = ['!/asl/packages/klayersV205/BinV201/klayers_airs fin=/asl/rtp/airs_l1c_v6/allfov/2020/003/allfov_ecmwf_airicrad_d2020003_216.rtp fout=junkecm.op.rtp >& ugh'];

junker = ['!/asl/packages/klayersV205/BinV201/klayers_airs fin=ERAorECM/allfov_ecmwf_airicrad_day_2016018_180.rtp fout=junkecm.op.rtp >& ugh'];
junker = ['!/asl/packages/klayersV205/BinV201/klayers_airs fin=ERAorECM/allfov_era_airicrad_day_2016018_180.rtp fout=junkecm.op.rtp >& ugh'];

eval(junker)
[hjunk,hajunk,pjunk,pajunk] = rtpread('junkecm.op.rtp');
mmwECM = mmwater_rtp(hjunk,pjunk);
scatter_coast(pjunk.rlon,pjunk.rlat,25,mmwECM); title('mmw ECM')
%}


if ~exist('era5')
  era5 = load('/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Center/ASC/era5_tile_center_16day_timestep_261.mat');
end
if ~exist('era5month')
  era5month = load('/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Center/ASC/era5_tile_center_monthly_137.mat');
end

era1  = load('/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA/Tile_Center/ASC/era_tile_center_timestep_261.mat'); 
era0  = load('/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA/Tile_Center/ASC/era_tile_center_timestep_261_closestINtime.mat'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; disp('comparing ERA5 : monthly vs daily')
clear ppmv* tz_* mmw* t*clr

figure(3); plot(era5.pnew_op.rlat,(era5.pnew_op.rtime-era5month.pnew_op.rtime)/86400); title('ERA5 : daily-monthly rtime diff in days')
%figure(3); scatter_coast(era5.pnew_op.rlon,era5.pnew_op.rlat,50,(era5.pnew_op.rtime-era5month.pnew_op.rtime)/86400); title('ERA5 : daily-monthly rtime diff in days')
%  caxis([-20 +20]); colormap(jet);
rett

figure(3); plot(era5.pnew_op.rlat,era5.pnew_op.stemp-era5month.pnew_op.stemp); title('ERA5 : daily-monthly stemp')
figure(3); scatter_coast(era5.pnew_op.rlon,era5.pnew_op.rlat,50,era5.pnew_op.stemp-era5month.pnew_op.stemp); title('ERA5 : daily-monthly stemp')
  caxis([-20 +20]); colormap(usa2);

mmw5 = mmwater_rtp(era5.hnew_op,era5.pnew_op);
mmw  = mmwater_rtp(era5month.hnew_op,era5month.pnew_op);
figure(6); scatter_coast(era5.pnew_op.rlon,era5.pnew_op.rlat,50,mmw5-mmw); title('ERA5 : daily-monthly mmw');  colormap(usa2);
figure(4); scatter_coast(era5.pnew_op.rlon,era5.pnew_op.rlat,50,mmw5); title('ERA5 daily mmw');   colormap(jet); caxis([0 70])
figure(5); scatter_coast(era5.pnew_op.rlon,era5.pnew_op.rlat,50,mmw);  title('ERA5 monthly mmw'); colormap(jet); caxis([0 70])
rett

t5clr = rad2bt(era5.hnew_op.vchan,era5.pnew_op.sarta_rclearcalc);
tclr  = rad2bt(era5.hnew_op.vchan,era5month.pnew_op.sarta_rclearcalc);
figure(5); scatter_coast(era5.pnew_op.rlon,era5.pnew_op.rlat,50,t5clr(1520,:)-tclr(1520,:)); title('ERA5 : daily-monthly BT1231')
  caxis([-20 +20]); colormap(usa2);

ocean = find(era5.pnew_op.landfrac == 0 & era5month.pnew_op.landfrac == 0);
land = find(era5.pnew_op.landfrac == 1 & era5month.pnew_op.landfrac == 1);
figure(6); plot(era5.hnew_op.vchan,nanmean(t5clr'-tclr'),'b',era5.hnew_op.vchan,nanstd(t5clr'-tclr'),'c'); title('Spectra Land/Ocean ERA5 : daily-monthly')
figure(6); plot(era5.hnew_op.vchan,nanmean(t5clr(:,ocean)'-tclr(:,ocean)'),'b',era5.hnew_op.vchan,nanstd(t5clr(:,ocean)'-tclr(:,ocean)'),'c'); title('Spectra Ocean ERA5 : daily-monthly')
rett

ppmv1_era5 = layers2ppmv(era5.hnew_op,era5.pnew_op,1:4608,1);
ppmv1_era  = layers2ppmv(era5month.hnew_op,era5month.pnew_op,1:4608,1);
figure(1); plot(nanmean(ppmv1_era5./ppmv1_era,2),1:98,1+nanstd(ppmv1_era5./ppmv1_era,[],2),1:98,'b--'); set(gca,'ydir','reverse');
  hl = legend('mean','std','location','best','fontsize',8); title('WV ppmv ratio ERA5 daily/monthly'); plotaxis2;

ppmv2_era5 = layers2ppmv(era5.hnew_op,era5.pnew_op,1:4608,2);
ppmv2_era  = layers2ppmv(era5month.hnew_op,era5month.pnew_op,1:4608,2);
figure(2); plot(nanmean(ppmv2_era5./ppmv2_era,2),1:98,1+nanstd(ppmv2_era5./ppmv2_era,[],2),1:98,'b--'); set(gca,'ydir','reverse'); xlim([0.995 1.005])
  hl = legend('mean','std','location','best','fontsize',8); title('CO2 ppmv ratio ERA5 daily/monthly'); plotaxis2;

ppmv3_era5 = layers2ppmv(era5.hnew_op,era5.pnew_op,1:4608,3);
ppmv3_era  = layers2ppmv(era5month.hnew_op,era5month.pnew_op,1:4608,3);
figure(3); plot(nanmean(ppmv3_era5./ppmv3_era,2),1:98,1+nanstd(ppmv3_era5./ppmv3_era,[],2),1:98,'b--'); set(gca,'ydir','reverse'); 
  hl = legend('mean','std','location','best','fontsize',8); title('O3 ppmv ratio ERA5 daily/monthly'); plotaxis2;

tz_era5 = era5.pnew_op.ptemp;
tz_era  = era5month.pnew_op.ptemp;
figure(4); plot(nanmean(tz_era5-tz_era,2),1:101,nanstd(tz_era5-tz_era,[],2),1:101,'b--'); set(gca,'ydir','reverse'); xlim([-10 +10])
  hl = legend('mean','std','location','best','fontsize',8); title('Tz diff ERA5 : daily-monthly'); plotaxis2;
rett

clear mmw5 mmw

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; disp('comparing ERA : interpolate in time vs closest in time')
clear ppmv* tz_* mmw* t*clr

figure(3); plot(era1.pnew_op.rlat,(era1.pnew_op.rtime-era0.pnew_op.rtime)/86400); title('ERA1-ERA rtime diff in days')
%figure(3); scatter_coast(era1.pnew_op.rlon,era1.pnew_op.rlat,50,(era1.pnew_op.rtime-era0.pnew_op.rtime)/86400); title('ERA1-ERA rtime diff in days')
%  caxis([-20 +20]); colormap(jet);
rett

figure(3); plot(era1.pnew_op.rlat,era1.pnew_op.stemp-era0.pnew_op.stemp); title('ERA1-ERA stemp')
figure(3); scatter_coast(era1.pnew_op.rlon,era1.pnew_op.rlat,50,era1.pnew_op.stemp-era0.pnew_op.stemp); title('ERA1-ERA stemp')
  caxis([-20 +20]); colormap(usa2);

mmw5 = mmwater_rtp(era1.hnew_op,era1.pnew_op);
mmw  = mmwater_rtp(era0.hnew_op,era0.pnew_op);
figure(6); scatter_coast(era1.pnew_op.rlon,era1.pnew_op.rlat,50,mmw5-mmw); title('ERA1-ERA mmw');  colormap(usa2);
figure(4); scatter_coast(era1.pnew_op.rlon,era1.pnew_op.rlat,50,mmw5); title('ERA1 mmw');  colormap(jet);  caxis([0 70])
figure(5); scatter_coast(era1.pnew_op.rlon,era1.pnew_op.rlat,50,mmw); title('ERA mmw');    colormap(jet);  caxis([0 70])
rett

t5clr = rad2bt(era1.hnew_op.vchan,era1.pnew_op.sarta_rclearcalc);
tclr  = rad2bt(era1.hnew_op.vchan,era0.pnew_op.sarta_rclearcalc);
figure(5); scatter_coast(era1.pnew_op.rlon,era1.pnew_op.rlat,50,t5clr(1520,:)-tclr(1520,:)); title('ERA1-ERA0 BT1231')
  caxis([-20 +20]); colormap(usa2);

ocean = find(era1.pnew_op.landfrac == 0 & era0.pnew_op.landfrac == 0);
land = find(era1.pnew_op.landfrac == 1 & era0.pnew_op.landfrac == 1);
figure(6); plot(era1.hnew_op.vchan,nanmean(t5clr'-tclr'),'b',era1.hnew_op.vchan,nanstd(t5clr'-tclr'),'c'); title('Spectra Land/Ocean ERA1-ERA0')
figure(6); plot(era1.hnew_op.vchan,nanmean(t5clr(:,ocean)'-tclr(:,ocean)'),'b',era1.hnew_op.vchan,nanstd(t5clr(:,ocean)'-tclr(:,ocean)'),'c'); title('Spectra Ocean ERA1-ERA0')
rett

ppmv1_era5 = layers2ppmv(era1.hnew_op,era1.pnew_op,1:4608,1);
ppmv1_era  = layers2ppmv(era0.hnew_op,era0.pnew_op,1:4608,1);
figure(1); plot(nanmean(ppmv1_era5./ppmv1_era,2),1:98,1+nanstd(ppmv1_era5./ppmv1_era,[],2),1:98,'b--'); set(gca,'ydir','reverse');
  hl = legend('mean','std','location','best','fontsize',8); title('WV ppmv ratio ERA1/ERA0'); plotaxis2;

ppmv2_era5 = layers2ppmv(era1.hnew_op,era1.pnew_op,1:4608,2);
ppmv2_era  = layers2ppmv(era0.hnew_op,era0.pnew_op,1:4608,2);
figure(2); plot(nanmean(ppmv2_era5./ppmv2_era,2),1:98,1+nanstd(ppmv2_era5./ppmv2_era,[],2),1:98,'b--'); set(gca,'ydir','reverse'); xlim([0.995 1.005])
  hl = legend('mean','std','location','best','fontsize',8); title('CO2 ppmv ratio ERA1/ERA0'); plotaxis2;

ppmv3_era5 = layers2ppmv(era1.hnew_op,era1.pnew_op,1:4608,3);
ppmv3_era  = layers2ppmv(era0.hnew_op,era0.pnew_op,1:4608,3);
figure(3); plot(nanmean(ppmv3_era5./ppmv3_era,2),1:98,1+nanstd(ppmv3_era5./ppmv3_era,[],2),1:98,'b--'); set(gca,'ydir','reverse'); 
  hl = legend('mean','std','location','best','fontsize',8); title('O3 ppmv ratio ERA1/ERA0'); plotaxis2;

tz_era5 = era1.pnew_op.ptemp;
tz_era  = era0.pnew_op.ptemp;
figure(4); plot(nanmean(tz_era5-tz_era,2),1:101,nanstd(tz_era5-tz_era,[],2),1:101,'b--'); set(gca,'ydir','reverse'); xlim([-10 +10])
  hl = legend('mean','std','location','best','fontsize',8); title('Tz diff ERA1-ERA0'); plotaxis2;
rett

clear mmw5 mmw

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; disp('comparing ERA5 daily vs ERA daily')
clear ppmv* tz_* mmw* t*clr

era = era1;

figure(1); plot(era5.pnew_op.rlat,era5.pnew_op.stemp,'b',era5.pnew_op.rlat,rad2bt(1231,era5.pnew_op.sarta_rclearcalc(1520,:)),'g',era5.pnew_op.rlat,rad2bt(1231,era5.pnew_op.rcalc(1520,:)),'r','linewidth',2)
    hl = legend('stemp','clear BT1231','allsky BT1231','location','best'); title('ERA5')

%figure(2); plot(era.pnew_op.rlat,era.pnew_op.stemp,'b',era.pnew_op.rlat,rad2bt(1231,era.pnew_op.sarta_rclearcalc(1520,:)),'g',era.pnew_op.rlat,rad2bt(1231,era.pnew_op.rcalc(1520,:)),'r','linewidth',2)
figure(2); plot(era.pnew_op.rlat,era.pnew_op.stemp,'b',era.pnew_op.rlat,rad2bt(1231,era.pnew_op.sarta_rclearcalc(1520,:)),'g',era.pnew_op.rlat,rad2bt(1231,era.pnew_op.sarta_rclearcalc(1520,:)),'r','linewidth',2)
    hl = legend('stemp','clear BT1231','allsky BT1231','location','best'); title('ERA')

figure(3); plot(era5.pnew_op.rlat,(era5.pnew_op.rtime-era.pnew_op.rtime)/86400); title('ERA5-ERA rtime diff in days')
%figure(3); scatter_coast(era5.pnew_op.rlon,era5.pnew_op.rlat,50,(era5.pnew_op.rtime-era.pnew_op.rtime)/86400); title('ERA5-ERA rtime diff in days')
%  caxis([-20 +20]); colormap(jet);
rett

figure(3); plot(era5.pnew_op.rlat,era5.pnew_op.stemp-era.pnew_op.stemp); title('ERA5-ERA stemp')
figure(3); scatter_coast(era5.pnew_op.rlon,era5.pnew_op.rlat,50,era5.pnew_op.stemp-era.pnew_op.stemp); title('ERA5-ERA stemp')
  caxis([-20 +20]); colormap(usa2);

mmw5 = mmwater_rtp(era5.hnew_op,era5.pnew_op);
mmw  = mmwater_rtp(era.hnew_op,era.pnew_op);
figure(6); scatter_coast(era5.pnew_op.rlon,era5.pnew_op.rlat,50,mmw5-mmw); title('ERA5-ERA mmw');  colormap(usa2);
figure(4); scatter_coast(era5.pnew_op.rlon,era5.pnew_op.rlat,50,mmw5); title('ERA5 mmw');  colormap(jet);  caxis([0 70])
figure(5); scatter_coast(era5.pnew_op.rlon,era5.pnew_op.rlat,50,mmw); title('ERA mmw');    colormap(jet);  caxis([0 70])

rett
ix = find(era5.pnew_op.rlon >= 150 & era5.pnew_op.rlat >= 0,1); [mmw(ix) mmw5(ix)]
semilogy(era.pnew_ip.ptemp(:,ix),era.pnew_ip.plevs(:,ix),'b.-',era5.pnew_ip.ptemp(:,ix),era5.pnew_ip.plevs(:,ix),'rx-'); set(gca,'ydir','reverse');
loglog(era.pnew_ip.gas_1(:,ix),era.pnew_ip.plevs(:,ix),'b.-',era5.pnew_ip.gas_1(:,ix),era5.pnew_ip.plevs(:,ix),'rx-'); set(gca,'ydir','reverse');
loglog(era.pnew_ip.gas_3(:,ix),era.pnew_ip.plevs(:,ix),'b.-',era5.pnew_ip.gas_3(:,ix),era5.pnew_ip.plevs(:,ix),'rx-'); set(gca,'ydir','reverse');
nlays = era.pnew_op.nlevs(ix)-1; nlays5 = era5.pnew_op.nlevs(ix)-1;
semilogy(era.pnew_op.ptemp(1:nlays,ix),era.pnew_op.plevs(1:nlays,ix),'b.-',era5.pnew_op.ptemp(1:nlays5,ix),era5.pnew_op.plevs(1:nlays5,ix),'rx-'); set(gca,'ydir','reverse');
loglog(era.pnew_op.gas_1(1:nlays,ix),era.pnew_op.plevs(1:nlays,ix),'b.-',era5.pnew_op.gas_1(1:nlays5,ix),era5.pnew_op.plevs(1:nlays5,ix),'rx-'); set(gca,'ydir','reverse');
loglog(era.pnew_op.gas_3(1:nlays,ix),era.pnew_op.plevs(1:nlays,ix),'b.-',era5.pnew_op.gas_3(1:nlays5,ix),era5.pnew_op.plevs(1:nlays5,ix),'rx-'); set(gca,'ydir','reverse');
rett

t5clr = rad2bt(era5.hnew_op.vchan,era5.pnew_op.sarta_rclearcalc);
tclr  = rad2bt(era5.hnew_op.vchan,era.pnew_op.sarta_rclearcalc);
figure(5); scatter_coast(era5.pnew_op.rlon,era5.pnew_op.rlat,50,t5clr(1520,:)-tclr(1520,:)); title('ERA5-ERA BT1231')
  caxis([-20 +20]); colormap(usa2);

ocean = find(era5.pnew_op.landfrac == 0 & era.pnew_op.landfrac == 0);
land = find(era5.pnew_op.landfrac == 1 & era.pnew_op.landfrac == 1);
figure(6); plot(era5.hnew_op.vchan,nanmean(t5clr'-tclr'),'b',era5.hnew_op.vchan,nanstd(t5clr'-tclr'),'c'); title('Spectra Land/Ocean ERA5-ERA')
  plotaxis2; hl = legend('bias','std','location','best','fontsize',10);
figure(6); plot(era5.hnew_op.vchan,nanmean(t5clr(:,ocean)'-tclr(:,ocean)'),'b',era5.hnew_op.vchan,nanstd(t5clr(:,ocean)'-tclr(:,ocean)'),'c'); title('Spectra Ocean ERA5-ERA')
  plotaxis2;   hl = legend('bias','std','location','best','fontsize',10);
rett

ppmv1_era5 = layers2ppmv(era5.hnew_op,era5.pnew_op,1:4608,1);
ppmv1_era  = layers2ppmv(era.hnew_op,era.pnew_op,1:4608,1);
figure(1); plot(nanmean(ppmv1_era5./ppmv1_era,2),1:98,1+nanstd(ppmv1_era5./ppmv1_era,[],2),1:98,'b--'); set(gca,'ydir','reverse');
  hl = legend('mean','std','location','best','fontsize',8); title('WV ppmv ratio ERA5/ERA'); plotaxis2;

ppmv2_era5 = layers2ppmv(era5.hnew_op,era5.pnew_op,1:4608,2);
ppmv2_era  = layers2ppmv(era.hnew_op,era.pnew_op,1:4608,2);
figure(2); plot(nanmean(ppmv2_era5./ppmv2_era,2),1:98,1+nanstd(ppmv2_era5./ppmv2_era,[],2),1:98,'b--'); set(gca,'ydir','reverse'); xlim([0.995 1.005])
  hl = legend('mean','std','location','best','fontsize',8); title('CO2 ppmv ratio ERA5/ERA'); plotaxis2;

ppmv3_era5 = layers2ppmv(era5.hnew_op,era5.pnew_op,1:4608,3);
ppmv3_era  = layers2ppmv(era.hnew_op,era.pnew_op,1:4608,3);
figure(3); plot(nanmean(ppmv3_era5./ppmv3_era,2),1:98,1+nanstd(ppmv3_era5./ppmv3_era,[],2),1:98,'b--'); set(gca,'ydir','reverse'); 
  hl = legend('mean','std','location','best','fontsize',8); title('O3 ppmv ratio ERA5/ERA'); plotaxis2;
xlim([0 3])
tz_era5 = era5.pnew_op.ptemp;
tz_era  = era.pnew_op.ptemp;
figure(4); plot(nanmean(tz_era5-tz_era,2),1:101,nanstd(tz_era5-tz_era,[],2),1:101,'b--'); set(gca,'ydir','reverse'); xlim([-10 +10])
  hl = legend('mean','std','location','best','fontsize',8); title('Tz diff ERA5-ERA'); plotaxis2;
rett


ix = 2500;
ix = 50 + 10*64;
semilogy(era.pnew_ip.ptemp(:,ix),era.pnew_ip.plevs(:,ix),'b.-',era5.pnew_ip.ptemp(:,ix),era5.pnew_ip.plevs(:,ix),'rx-'); set(gca,'ydir','reverse');
loglog(era.pnew_ip.gas_1(:,ix),era.pnew_ip.plevs(:,ix),'b.-',era5.pnew_ip.gas_1(:,ix),era5.pnew_ip.plevs(:,ix),'rx-'); set(gca,'ydir','reverse');
loglog(era.pnew_ip.gas_3(:,ix),era.pnew_ip.plevs(:,ix),'b.-',era5.pnew_ip.gas_3(:,ix),era5.pnew_ip.plevs(:,ix),'rx-'); set(gca,'ydir','reverse');

