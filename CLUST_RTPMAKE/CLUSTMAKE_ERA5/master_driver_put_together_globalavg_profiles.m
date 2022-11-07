addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS/
addpath /asl/matlib/h4tools
addpath /asl/matlib/rtptools

simulateYear = 2012;

if ~exist('allprof11')
  for ii = 1 : 64
    if mod(ii,10) == 0
      fprintf(1,'+')
    else
      fprintf(1,'.')
    end
    fname = ['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Center12months/DESC/' num2str(simulateYear,'%04d') '/era5_full12months_latbin_' num2str(ii,'%02d') '_tile_center_profilesQcumulative_1_11.mat'];
    a = load(fname);
    if ii == 1
      hnew_op = a.hnew_op;
      allprof01 = a.prof01;
      allprof02 = a.prof02;
      allprof03 = a.prof03;
      allprof04 = a.prof04;
      allprof05 = a.prof05;
      allprof06 = a.prof06;
      allprof07 = a.prof07;
      allprof08 = a.prof08;
      allprof09 = a.prof09;
      allprof10 = a.prof10;
      allprof11 = a.prof11;
    else
      [~,allprof01] = cat_rtp(hnew_op,allprof01,hnew_op,a.prof01);
      [~,allprof02] = cat_rtp(hnew_op,allprof02,hnew_op,a.prof02);
      [~,allprof03] = cat_rtp(hnew_op,allprof03,hnew_op,a.prof03);
      [~,allprof04] = cat_rtp(hnew_op,allprof04,hnew_op,a.prof04);
      [~,allprof05] = cat_rtp(hnew_op,allprof05,hnew_op,a.prof05);
      [~,allprof06] = cat_rtp(hnew_op,allprof06,hnew_op,a.prof06);
      [~,allprof07] = cat_rtp(hnew_op,allprof07,hnew_op,a.prof07);
      [~,allprof08] = cat_rtp(hnew_op,allprof08,hnew_op,a.prof08);
      [~,allprof09] = cat_rtp(hnew_op,allprof09,hnew_op,a.prof09);
      [~,allprof10] = cat_rtp(hnew_op,allprof10,hnew_op,a.prof10);
      [~,allprof11] = cat_rtp(hnew_op,allprof11,hnew_op,a.prof11);
    end
  end
  clear a
  fprintf(1,'\n');
  
  rtpwrite(['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Center12months/DESC/' num2str(simulateYear,'%04d') '/all4608_era5_full12months_Qcumulative01.rtp'],hnew_op,[],allprof01,[]);
  rtpwrite(['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Center12months/DESC/' num2str(simulateYear,'%04d') '/all4608_era5_full12months_Qcumulative02.rtp'],hnew_op,[],allprof02,[]);
  rtpwrite(['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Center12months/DESC/' num2str(simulateYear,'%04d') '/all4608_era5_full12months_Qcumulative03.rtp'],hnew_op,[],allprof03,[]);
  rtpwrite(['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Center12months/DESC/' num2str(simulateYear,'%04d') '/all4608_era5_full12months_Qcumulative04.rtp'],hnew_op,[],allprof04,[]);
  rtpwrite(['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Center12months/DESC/' num2str(simulateYear,'%04d') '/all4608_era5_full12months_Qcumulative05.rtp'],hnew_op,[],allprof05,[]);
  rtpwrite(['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Center12months/DESC/' num2str(simulateYear,'%04d') '/all4608_era5_full12months_Qcumulative06.rtp'],hnew_op,[],allprof06,[]);
  rtpwrite(['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Center12months/DESC/' num2str(simulateYear,'%04d') '/all4608_era5_full12months_Qcumulative07.rtp'],hnew_op,[],allprof07,[]);
  rtpwrite(['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Center12months/DESC/' num2str(simulateYear,'%04d') '/all4608_era5_full12months_Qcumulative08.rtp'],hnew_op,[],allprof08,[]);
  rtpwrite(['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Center12months/DESC/' num2str(simulateYear,'%04d') '/all4608_era5_full12months_Qcumulative09.rtp'],hnew_op,[],allprof09,[]);
  rtpwrite(['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Center12months/DESC/' num2str(simulateYear,'%04d') '/all4608_era5_full12months_Qcumulative10.rtp'],hnew_op,[],allprof10,[]);
  rtpwrite(['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Center12months/DESC/' num2str(simulateYear,'%04d') '/all4608_era5_full12months_Qcumulative11.rtp'],hnew_op,[],allprof11,[]);
end
  
RH01 = layeramt2RH(hnew_op,allprof01); co2ppmv01 = layers2ppmv(hnew_op,allprof01,1:length(allprof01.stemp),2);
RH02 = layeramt2RH(hnew_op,allprof02); co2ppmv02 = layers2ppmv(hnew_op,allprof02,1:length(allprof02.stemp),2);
RH03 = layeramt2RH(hnew_op,allprof03); co2ppmv03 = layers2ppmv(hnew_op,allprof03,1:length(allprof03.stemp),2);
RH04 = layeramt2RH(hnew_op,allprof04); co2ppmv04 = layers2ppmv(hnew_op,allprof04,1:length(allprof04.stemp),2);
RH05 = layeramt2RH(hnew_op,allprof05); co2ppmv05 = layers2ppmv(hnew_op,allprof05,1:length(allprof05.stemp),2);
RH06 = layeramt2RH(hnew_op,allprof06); co2ppmv06 = layers2ppmv(hnew_op,allprof06,1:length(allprof06.stemp),2);
RH07 = layeramt2RH(hnew_op,allprof07); co2ppmv07 = layers2ppmv(hnew_op,allprof07,1:length(allprof07.stemp),2);
RH08 = layeramt2RH(hnew_op,allprof08); co2ppmv08 = layers2ppmv(hnew_op,allprof08,1:length(allprof08.stemp),2);
RH09 = layeramt2RH(hnew_op,allprof09); co2ppmv09 = layers2ppmv(hnew_op,allprof09,1:length(allprof09.stemp),2);
RH10 = layeramt2RH(hnew_op,allprof10); co2ppmv10 = layers2ppmv(hnew_op,allprof10,1:length(allprof10.stemp),2);
RH11 = layeramt2RH(hnew_op,allprof11); co2ppmv11 = layers2ppmv(hnew_op,allprof11,1:length(allprof11.stemp),2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[hold,haold,pold,paold] = rtpread('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/summary_19years_all_lat_all_lon_2002_2021_monthlyERA5.rp.rtp');
figure(1); scatter_coast(pold.rlon,pold.rlat,50,pold.stemp);      colormap(jet); title('Orig Stemp, used in STM'); caxis([220 330])
figure(2); scatter_coast(pold.rlon,pold.rlat,50,allprof11.stemp); colormap(jet); title('New Q11 Stemp'); caxis([220 330])
figure(3); scatter_coast(pold.rlon,pold.rlat,50,pold.stemp-allprof11.stemp); caxis([-1 +1]*20); title('OrigStemp - NewQ11 stemp')
figure(4); scatter_coast(pold.rlon,pold.rlat,50,rad2bt(1231,allprof11.rcalc(1520,:))); colormap(jet); title('New Q11 BT1231 cld cal'); caxis([220 330])
figure(5); scatter_coast(pold.rlon,pold.rlat,50,rad2bt(1231,allprof01.rcalc(1520,:))); colormap(jet); title('New Q01 BT1231 cld cal'); caxis([220 330])
figure(6); scatter_coast(pold.rlon,pold.rlat,50,rad2bt(1231,allprof11.sarta_rclearcalc(1520,:))); colormap(jet); title('New Q11 BT1231 clr cal'); caxis([220 330])

figure(1); caxis([220 330]); figure(2); caxis([220 330]); figure(4); caxis([220 330]); figure(5); caxis([220 330]); figure(6); caxis([220 330]); 
figure(1); caxis([200 300]); figure(2); caxis([200 300]); figure(4); caxis([200 300]); figure(5); caxis([200 300]); figure(6); caxis([200 300]); 
disp('ret to continue'); pause;

clr1231woo00 = reshape(pold.rcalc(1520,:),72,64); clr1231woo00(clr1231woo00 < 0) = NaN; clr1231woo00 = rad2bt(1231,squeeze(nanmean(clr1231woo00,1))); %%% <<<< this is clr N cld

stwoo11 = reshape(allprof11.stemp,72,64); stwoo11(stwoo11 < 0) = NaN; stwoo11 = squeeze(nanmean(stwoo11,1));
stwoo10 = reshape(allprof10.stemp,72,64); stwoo10(stwoo10 < 0) = NaN; stwoo10 = squeeze(nanmean(stwoo10,1));
stwoo09 = reshape(allprof09.stemp,72,64); stwoo09(stwoo09 < 0) = NaN; stwoo09 = squeeze(nanmean(stwoo09,1));
stwoo08 = reshape(allprof08.stemp,72,64); stwoo08(stwoo08 < 0) = NaN; stwoo08 = squeeze(nanmean(stwoo08,1));
stwoo07 = reshape(allprof07.stemp,72,64); stwoo07(stwoo07 < 0) = NaN; stwoo07 = squeeze(nanmean(stwoo07,1));
stwoo01 = reshape(allprof01.stemp,72,64); stwoo01(stwoo01 < 0) = NaN; stwoo01 = squeeze(nanmean(stwoo01,1));

clr1231woo11 = reshape(allprof11.sarta_rclearcalc(1520,:),72,64); clr1231woo11(clr1231woo11 < 0) = NaN; clr1231woo11 = rad2bt(1231,squeeze(nanmean(clr1231woo11,1)));
clr1231woo10 = reshape(allprof10.sarta_rclearcalc(1520,:),72,64); clr1231woo10(clr1231woo10 < 0) = NaN; clr1231woo10 = rad2bt(1231,squeeze(nanmean(clr1231woo10,1)));
clr1231woo09 = reshape(allprof09.sarta_rclearcalc(1520,:),72,64); clr1231woo09(clr1231woo09 < 0) = NaN; clr1231woo09 = rad2bt(1231,squeeze(nanmean(clr1231woo09,1)));
clr1231woo08 = reshape(allprof08.sarta_rclearcalc(1520,:),72,64); clr1231woo08(clr1231woo08 < 0) = NaN; clr1231woo08 = rad2bt(1231,squeeze(nanmean(clr1231woo08,1)));
clr1231woo07 = reshape(allprof07.sarta_rclearcalc(1520,:),72,64); clr1231woo07(clr1231woo07 < 0) = NaN; clr1231woo07 = rad2bt(1231,squeeze(nanmean(clr1231woo07,1)));
clr1231woo01 = reshape(allprof01.sarta_rclearcalc(1520,:),72,64); clr1231woo01(clr1231woo01 < 0) = NaN; clr1231woo01 = rad2bt(1231,squeeze(nanmean(clr1231woo01,1)));
figure(1); plot(1:64,clr1231woo00,'k+-',1:64,clr1231woo11,1:64,clr1231woo10,1:64,clr1231woo09,1:64,clr1231woo08,1:64,clr1231woo07,1:64,clr1231woo01); title('BT1231 clr calcs')

cld1231woo11 = reshape(allprof11.rcalc(1520,:),72,64); cld1231woo11(cld1231woo11 < 0) = NaN; cld1231woo11 = rad2bt(1231,squeeze(nanmean(cld1231woo11,1)));
cld1231woo10 = reshape(allprof10.rcalc(1520,:),72,64); cld1231woo10(cld1231woo10 < 0) = NaN; cld1231woo10 = rad2bt(1231,squeeze(nanmean(cld1231woo10,1)));
cld1231woo09 = reshape(allprof09.rcalc(1520,:),72,64); cld1231woo09(cld1231woo09 < 0) = NaN; cld1231woo09 = rad2bt(1231,squeeze(nanmean(cld1231woo09,1)));
cld1231woo08 = reshape(allprof08.rcalc(1520,:),72,64); cld1231woo08(cld1231woo08 < 0) = NaN; cld1231woo08 = rad2bt(1231,squeeze(nanmean(cld1231woo08,1)));
cld1231woo07 = reshape(allprof07.rcalc(1520,:),72,64); cld1231woo07(cld1231woo07 < 0) = NaN; cld1231woo07 = rad2bt(1231,squeeze(nanmean(cld1231woo07,1)));
cld1231woo01 = reshape(allprof01.rcalc(1520,:),72,64); cld1231woo01(cld1231woo01 < 0) = NaN; cld1231woo01 = rad2bt(1231,squeeze(nanmean(cld1231woo01,1)));
figure(2); plot(1:64,clr1231woo00,'k+-',1:64,cld1231woo11,1:64,cld1231woo10,1:64,cld1231woo09,1:64,cld1231woo08,1:64,cld1231woo07,1:64,cld1231woo01); title('BT1231 cld calcs')

figure(3); plot(1:64,clr1231woo11,'b',1:64,cld1231woo11,'r',1:64,stwoo11,'k','linewidth',2); hl = legend('Q11 clr','Q11 cld','Q11 stemp'); title('BT1231'); ylim([240 310])
figure(4); plot(1:64,clr1231woo01,'b',1:64,cld1231woo01,'r',1:64,stwoo01,'k','linewidth',2); hl = legend('Q01 clr','Q01 cld','Q01 stemp'); title('BT1231'); ylim([240 310])
disp('ret to continue'); pause;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% remember, we are looking at populating quantiles then averaging. so prof11.rtime is NOT same as prof10.rtime etc
[yy11,mm11,dd11,hh11] = tai2utcSergio(allprof11.rtime);
[yy10,mm10,dd10,hh10] = tai2utcSergio(allprof10.rtime);
[yy09,mm09,dd09,hh09] = tai2utcSergio(allprof09.rtime);
[yy08,mm08,dd08,hh08] = tai2utcSergio(allprof08.rtime);
[yy07,mm07,dd07,hh07] = tai2utcSergio(allprof07.rtime);
[yy01,mm01,dd01,hh01] = tai2utcSergio(allprof01.rtime);

figure(1); pcolor_coast(allprof11.rlon,allprof11.rlat,allprof11.rtime); colormap jet; title('rtime')
figure(2); pcolor_coast(allprof11.rlon,allprof11.rlat,allprof11.rtime-allprof10.rtime); colormap usa2; caxis([-1 +1]*200)
figure(3); pcolor_coast(allprof11.rlon,allprof11.rlat,allprof11.rtime-allprof09.rtime); colormap usa2; caxis([-1 +1]*200)
figure(4); pcolor_coast(allprof11.rlon,allprof11.rlat,allprof11.rtime-allprof06.rtime); colormap usa2; caxis([-1 +1]*200)
figure(5); pcolor_coast(allprof11.rlon,allprof11.rlat,allprof11.rtime-allprof07.rtime); colormap usa2; caxis([-1 +1]*200)
disp('showing diffs in rtime (seconds) between the quantiles .. no need to worry since have to mix and match different simulations to create quantiles');
disp('ret to continue'); pause;

figure(1); scatter_coast(allprof11.rlon,allprof11.rlat,50,mm11+dd11/30); colormap jet; title('Q11 hottest month')
figure(2); scatter_coast(allprof11.rlon,allprof11.rlat,50,mm10+dd10/30); colormap jet; title('Q10 hottest month')
figure(3); scatter_coast(allprof11.rlon,allprof11.rlat,50,mm08+dd09/30); colormap jet; title('Q09 hottest month')
figure(4); scatter_coast(allprof11.rlon,allprof11.rlat,50,mm07+dd08/30); colormap jet; title('Q08 hottest month')
figure(5); scatter_coast(allprof11.rlon,allprof11.rlat,50,mm01+dd01/30); colormap jet; title('Q01 hottest month')
disp('ret to continue'); pause;

daysSinceX = change2days(yy11,mm11,dd11,simulateYear); figure(1); scatter_coast(allprof11.rlon,allprof11.rlat,50,daysSinceX); colormap jet; title(['Q11 : days since Jan 01 ' num2str(simulateYear)])
daysSinceX = change2days(yy10,mm10,dd10,simulateYear); figure(2); scatter_coast(allprof10.rlon,allprof10.rlat,50,daysSinceX); colormap jet; title(['Q10 : days since Jan 01 ' num2str(simulateYear)])
daysSinceX = change2days(yy09,mm09,dd09,simulateYear); figure(3); scatter_coast(allprof09.rlon,allprof09.rlat,50,daysSinceX); colormap jet; title(['Q09 : days since Jan 01 ' num2str(simulateYear)])
daysSinceX = change2days(yy08,mm08,dd08,simulateYear); figure(4); scatter_coast(allprof08.rlon,allprof08.rlat,50,daysSinceX); colormap jet; title(['Q08 : days since Jan 01 ' num2str(simulateYear)])
daysSinceX = change2days(yy07,mm07,dd07,simulateYear); figure(5); scatter_coast(allprof07.rlon,allprof07.rlat,50,daysSinceX); colormap jet; title(['Q07 : days since Jan 01 ' num2str(simulateYear)])
daysSinceX = change2days(yy01,mm01,dd01,simulateYear); figure(6); scatter_coast(allprof01.rlon,allprof01.rlat,50,daysSinceX); colormap jet; title(['Q01 : days since Jan 01 ' num2str(simulateYear)])
disp('ret to continue'); pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1); scatter_coast(allprof11.rlon,allprof11.rlat,50,allprof11.stemp); colormap jet
figure(2); scatter_coast(allprof11.rlon,allprof11.rlat,50,allprof11.stemp-allprof10.stemp); colormap usa2; caxis([-1 +1]*5); title('stemp : Q(1.00)-Q(0.97)')
figure(3); scatter_coast(allprof11.rlon,allprof11.rlat,50,allprof11.stemp-allprof09.stemp); colormap usa2; caxis([-1 +1]*5); title('stemp : Q(1.00)-Q(0.95)')
figure(4); scatter_coast(allprof11.rlon,allprof11.rlat,50,allprof11.stemp-allprof06.stemp); colormap usa2; caxis([-1 +1]*5); title('stemp : Q(1.00)-Q(0.90)')
figure(5); scatter_coast(allprof11.rlon,allprof11.rlat,50,allprof11.stemp-allprof07.stemp); colormap usa2; caxis([-1 +1]*5); title('stemp : Q(1.00)-Q(0.80)')

figure(1); pcolor_coast(allprof11.rlon,allprof11.rlat,allprof11.stemp); colormap jet
figure(2); pcolor_coast(allprof11.rlon,allprof11.rlat,allprof11.stemp-allprof10.stemp); colormap usa2; caxis([-1 +1]*5); title('stemp : Q(1.00)-Q(0.97)')
figure(3); pcolor_coast(allprof11.rlon,allprof11.rlat,allprof11.stemp-allprof09.stemp); colormap usa2; caxis([-1 +1]*5); title('stemp : Q(1.00)-Q(0.95)')
figure(4); pcolor_coast(allprof11.rlon,allprof11.rlat,allprof11.stemp-allprof06.stemp); colormap usa2; caxis([-1 +1]*5); title('stemp : Q(1.00)-Q(0.90)')
figure(5); pcolor_coast(allprof11.rlon,allprof11.rlat,allprof11.stemp-allprof07.stemp); colormap usa2; caxis([-1 +1]*5); title('stemp : Q(1.00)-Q(0.80)')
disp('ret to continue'); pause;

woo11 = reshape(allprof11.ptemp,101,72,64); woo11(woo11 < 150) = NaN; woo11 = squeeze(nanmean(woo11,2));
woo10 = reshape(allprof10.ptemp,101,72,64); woo10(woo10 < 150) = NaN; woo10 = squeeze(nanmean(woo10,2));
woo09 = reshape(allprof09.ptemp,101,72,64); woo09(woo09 < 150) = NaN; woo09 = squeeze(nanmean(woo09,2));
woo08 = reshape(allprof08.ptemp,101,72,64); woo08(woo08 < 150) = NaN; woo08 = squeeze(nanmean(woo08,2));
woo07 = reshape(allprof07.ptemp,101,72,64); woo07(woo07 < 150) = NaN; woo07 = squeeze(nanmean(woo07,2));
figure(1); pcolor(woo11); colorbar; shading interp; set(gca,'ydir','reverse');
figure(2); pcolor(woo11-woo10); colorbar; shading interp; set(gca,'ydir','reverse'); caxis([-1 +1]);  title('T(z,lat) : Q(1.00)-Q(0.97)')
figure(3); pcolor(woo11-woo09); colorbar; shading interp; set(gca,'ydir','reverse'); caxis([-1 +1]);  title('T(z,lat) : Q(1.00)-Q(0.95)')
figure(4); pcolor(woo11-woo08); colorbar; shading interp; set(gca,'ydir','reverse'); caxis([-1 +1]);  title('T(z,lat) : Q(1.00)-Q(0.90)')
figure(5); pcolor(woo11-woo07); colorbar; shading interp; set(gca,'ydir','reverse'); caxis([-1 +1]);  title('T(z,lat) : Q(1.00)-Q(0.80)')
disp('ret to continue'); pause;

woo11 = reshape(allprof11.gas_1,101,72,64); woo11(woo11 < 0) = NaN; woo11 = squeeze(nanmean(woo11,2));
woo10 = reshape(allprof10.gas_1,101,72,64); woo10(woo10 < 0) = NaN; woo10 = squeeze(nanmean(woo10,2));
woo09 = reshape(allprof09.gas_1,101,72,64); woo09(woo09 < 0) = NaN; woo09 = squeeze(nanmean(woo09,2));
woo08 = reshape(allprof08.gas_1,101,72,64); woo08(woo08 < 0) = NaN; woo08 = squeeze(nanmean(woo08,2));
woo07 = reshape(allprof07.gas_1,101,72,64); woo07(woo07 < 0) = NaN; woo07 = squeeze(nanmean(woo07,2));
figure(1); pcolor(log10(woo11)); colorbar; shading interp; set(gca,'ydir','reverse');
figure(2); pcolor(woo11./woo10-1); colorbar; shading interp; set(gca,'ydir','reverse'); caxis([-1 +1]*0.05);  title('WV(z,lat) : Q(1.00)./Q(0.97)-1')
figure(3); pcolor(woo11./woo09-1); colorbar; shading interp; set(gca,'ydir','reverse'); caxis([-1 +1]*0.05);  title('WV(z,lat) : Q(1.00)./Q(0.95)-1')
figure(4); pcolor(woo11./woo08-1); colorbar; shading interp; set(gca,'ydir','reverse'); caxis([-1 +1]*0.05);  title('WV(z,lat) : Q(1.00)./Q(0.90)-1')
figure(5); pcolor(woo11./woo07-1); colorbar; shading interp; set(gca,'ydir','reverse'); caxis([-1 +1]*0.05);  title('WV(z,lat) : Q(1.00)./Q(0.80)-1')
disp('ret to continue'); pause;

woo11 = reshape(RH11,100,72,64); woo11(woo11 < 0) = NaN; woo11 = squeeze(nanmean(woo11,2));
woo10 = reshape(RH10,100,72,64); woo10(woo10 < 0) = NaN; woo10 = squeeze(nanmean(woo10,2));
woo09 = reshape(RH09,100,72,64); woo09(woo09 < 0) = NaN; woo09 = squeeze(nanmean(woo09,2));
woo08 = reshape(RH08,100,72,64); woo08(woo08 < 0) = NaN; woo08 = squeeze(nanmean(woo08,2));
woo07 = reshape(RH07,100,72,64); woo07(woo07 < 0) = NaN; woo07 = squeeze(nanmean(woo07,2));
figure(1); pcolor(woo11); colorbar; shading interp; set(gca,'ydir','reverse');
figure(2); pcolor(woo11-woo10); colorbar; shading interp; set(gca,'ydir','reverse'); caxis([-1 +1]);  title('RH(z,lat) : Q(1.00)-Q(0.97)')
figure(3); pcolor(woo11-woo09); colorbar; shading interp; set(gca,'ydir','reverse'); caxis([-1 +1]);  title('RH(z,lat) : Q(1.00)-Q(0.95)')
figure(4); pcolor(woo11-woo08); colorbar; shading interp; set(gca,'ydir','reverse'); caxis([-1 +1]);  title('RH(z,lat) : Q(1.00)-Q(0.90)')
figure(5); pcolor(woo11-woo07); colorbar; shading interp; set(gca,'ydir','reverse'); caxis([-1 +1]);  title('RH(z,lat) : Q(1.00)-Q(0.80)')
disp('ret to continue'); pause;

woo11 = reshape(co2ppmv11,98,72,64); woo11(woo11 < 0) = NaN; woo11 = squeeze(nanmean(woo11,2));
woo10 = reshape(co2ppmv10,98,72,64); woo10(woo10 < 0) = NaN; woo10 = squeeze(nanmean(woo10,2));
woo09 = reshape(co2ppmv09,98,72,64); woo09(woo09 < 0) = NaN; woo09 = squeeze(nanmean(woo09,2));
woo08 = reshape(co2ppmv08,98,72,64); woo08(woo08 < 0) = NaN; woo08 = squeeze(nanmean(woo08,2));
woo07 = reshape(co2ppmv07,98,72,64); woo07(woo07 < 0) = NaN; woo07 = squeeze(nanmean(woo07,2));
figure(1); pcolor(woo11); colorbar; shading interp; set(gca,'ydir','reverse');
figure(2); pcolor(woo11-woo10); colorbar; shading interp; set(gca,'ydir','reverse'); caxis([-1 +1]);  title('T(z,lat) : co2ppmv(1.00)-co2ppmv(0.97)')
figure(3); pcolor(woo11-woo09); colorbar; shading interp; set(gca,'ydir','reverse'); caxis([-1 +1]);  title('T(z,lat) : co2ppmv(1.00)-co2ppmv(0.95)')
figure(4); pcolor(woo11-woo08); colorbar; shading interp; set(gca,'ydir','reverse'); caxis([-1 +1]);  title('T(z,lat) : co2ppmv(1.00)-co2ppmv(0.90)')
figure(5); pcolor(woo11-woo07); colorbar; shading interp; set(gca,'ydir','reverse'); caxis([-1 +1]);  title('T(z,lat) : co2ppmv(1.00)-co2ppmv(0.80)')
disp('ret to continue'); pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
woo11 = reshape(allprof11.iceOD,72,64); woo11(woo11 < 0) = NaN; woo11 = squeeze(nanmean(woo11,1));
woo10 = reshape(allprof10.iceOD,72,64); woo10(woo10 < 0) = NaN; woo10 = squeeze(nanmean(woo10,1));
woo09 = reshape(allprof09.iceOD,72,64); woo09(woo09 < 0) = NaN; woo09 = squeeze(nanmean(woo09,1));
woo08 = reshape(allprof08.iceOD,72,64); woo08(woo08 < 0) = NaN; woo08 = squeeze(nanmean(woo08,1));
woo07 = reshape(allprof07.iceOD,72,64); woo07(woo07 < 0) = NaN; woo07 = squeeze(nanmean(woo07,1));
woo01 = reshape(allprof01.iceOD,72,64); woo01(woo01 < 0) = NaN; woo01 = squeeze(nanmean(woo01,1));
figure(1); plot(1:64,woo11,1:64,woo10,1:64,woo09,1:64,woo08,1:64,woo07,1:64,woo01); title('iceOD')

woo11 = reshape(allprof11.icefrac,72,64); woo11(woo11 < 0) = NaN; woo11 = squeeze(nanmean(woo11,1));
woo10 = reshape(allprof10.icefrac,72,64); woo10(woo10 < 0) = NaN; woo10 = squeeze(nanmean(woo10,1));
woo09 = reshape(allprof09.icefrac,72,64); woo09(woo09 < 0) = NaN; woo09 = squeeze(nanmean(woo09,1));
woo08 = reshape(allprof08.icefrac,72,64); woo08(woo08 < 0) = NaN; woo08 = squeeze(nanmean(woo08,1));
woo07 = reshape(allprof07.icefrac,72,64); woo07(woo07 < 0) = NaN; woo07 = squeeze(nanmean(woo07,1));
woo01 = reshape(allprof01.icefrac,72,64); woo01(woo01 < 0) = NaN; woo01 = squeeze(nanmean(woo01,1));
figure(2); plot(1:64,woo11,1:64,woo10,1:64,woo09,1:64,woo08,1:64,woo07,1:64,woo01); title('ice frac')

woo11 = reshape(allprof11.icetop,72,64); woo11(woo11 < 0) = NaN; woo11 = squeeze(nanmean(woo11,1));
woo10 = reshape(allprof10.icetop,72,64); woo10(woo10 < 0) = NaN; woo10 = squeeze(nanmean(woo10,1));
woo09 = reshape(allprof09.icetop,72,64); woo09(woo09 < 0) = NaN; woo09 = squeeze(nanmean(woo09,1));
woo08 = reshape(allprof08.icetop,72,64); woo08(woo08 < 0) = NaN; woo08 = squeeze(nanmean(woo08,1));
woo07 = reshape(allprof07.icetop,72,64); woo07(woo07 < 0) = NaN; woo07 = squeeze(nanmean(woo07,1));
woo01 = reshape(allprof01.icetop,72,64); woo01(woo01 < 0) = NaN; woo01 = squeeze(nanmean(woo01,1));
figure(3); plot(1:64,woo11,1:64,woo10,1:64,woo09,1:64,woo08,1:64,woo07,1:64,woo01); title('ice top (mb)'); set(gca,'ydir','reverse');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
woo11 = reshape(allprof11.waterOD,72,64); woo11(woo11 < 0) = NaN; woo11 = squeeze(nanmean(woo11,1));
woo10 = reshape(allprof10.waterOD,72,64); woo10(woo10 < 0) = NaN; woo10 = squeeze(nanmean(woo10,1));
woo09 = reshape(allprof09.waterOD,72,64); woo09(woo09 < 0) = NaN; woo09 = squeeze(nanmean(woo09,1));
woo08 = reshape(allprof08.waterOD,72,64); woo08(woo08 < 0) = NaN; woo08 = squeeze(nanmean(woo08,1));
woo07 = reshape(allprof07.waterOD,72,64); woo07(woo07 < 0) = NaN; woo07 = squeeze(nanmean(woo07,1));
woo01 = reshape(allprof01.waterOD,72,64); woo01(woo01 < 0) = NaN; woo01 = squeeze(nanmean(woo01,1));
figure(4); plot(1:64,woo11,1:64,woo10,1:64,woo09,1:64,woo08,1:64,woo07,1:64,woo01); title('waterOD')

woo11 = reshape(allprof11.waterfrac,72,64); woo11(woo11 < 0) = NaN; woo11 = squeeze(nanmean(woo11,1));
woo10 = reshape(allprof10.waterfrac,72,64); woo10(woo10 < 0) = NaN; woo10 = squeeze(nanmean(woo10,1));
woo09 = reshape(allprof09.waterfrac,72,64); woo09(woo09 < 0) = NaN; woo09 = squeeze(nanmean(woo09,1));
woo08 = reshape(allprof08.waterfrac,72,64); woo08(woo08 < 0) = NaN; woo08 = squeeze(nanmean(woo08,1));
woo07 = reshape(allprof07.waterfrac,72,64); woo07(woo07 < 0) = NaN; woo07 = squeeze(nanmean(woo07,1));
woo01 = reshape(allprof01.waterfrac,72,64); woo01(woo01 < 0) = NaN; woo01 = squeeze(nanmean(woo01,1));
figure(5); plot(1:64,woo11,1:64,woo10,1:64,woo09,1:64,woo08,1:64,woo07,1:64,woo01); title('water frac')

woo11 = reshape(allprof11.watertop,72,64); woo11(woo11 < 0) = NaN; woo11 = squeeze(nanmean(woo11,1));
woo10 = reshape(allprof10.watertop,72,64); woo10(woo10 < 0) = NaN; woo10 = squeeze(nanmean(woo10,1));
woo09 = reshape(allprof09.watertop,72,64); woo09(woo09 < 0) = NaN; woo09 = squeeze(nanmean(woo09,1));
woo08 = reshape(allprof08.watertop,72,64); woo08(woo08 < 0) = NaN; woo08 = squeeze(nanmean(woo08,1));
woo07 = reshape(allprof07.watertop,72,64); woo07(woo07 < 0) = NaN; woo07 = squeeze(nanmean(woo07,1));
woo01 = reshape(allprof01.watertop,72,64); woo01(woo01 < 0) = NaN; woo01 = squeeze(nanmean(woo01,1));
figure(6); plot(1:64,woo11,1:64,woo10,1:64,woo09,1:64,woo08,1:64,woo07,1:64,woo01); title('water top (mb)'); set(gca,'ydir','reverse');
disp('ret to continue'); pause;
