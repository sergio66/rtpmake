%% one per month, 19 years of AIRS data so 19x12 = 228 sets of data
%% one per month, 19 years of AIRS data so 20x12 = 240 sets of data

starttime = tic;
JOB = str2num(getenv('SLURM_ARRAY_TASK_ID'));  %% JOB = 1 : 4608
%JOB = 11 + 72*0;
%JOB = 11 + 72*3;
%JOB = 11 + 72*32
%JOB = 2710
%JOB = 3620

if JOB > 4608
  error('JOB = 1 .. 4608')
end

addpath /home/sergio/MATLABCODE/TIME
addpath /asl/matlib/rtptools/
addpath /asl/matlib/aslutil
addpath /asl/matlib/h4tools
addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/TIME
addpath /home/sergio/MATLABCODE/PLOTTER
addpath ../GRIB
addpath /home/sergio/MATLABCODE/matlib/clouds/sarta
addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS
addpath /home/sergio/MATLABCODE/CRODGERS_FAST_CLOUD

system_slurm_stats

addpath /asl/matlib/aslutil
addpath /asl/packages/time
addpath /asl/matlib/science

simulateYear = 2012;

load latB64.mat
rlat65 = latB2; rlon73 = -180 : 5 : +180;
rlon = -180 : 5 : +180;  rlat = latB2;
rlon = 0.5*(rlon(1:end-1)+rlon(2:end));
rlat = 0.5*(rlat(1:end-1)+rlat(2:end));

[Y,X] = meshgrid(rlat,rlon);
Y = Y(:); X = X(:);

% for JOB = 1 : 4608
%   tileXcenter = X(JOB);
%   tileYcenter = Y(JOB);
%   xindX = find(rlon73 < tileXcenter); xindX = xindX(end);
%   xindY = find(rlat65 < tileYcenter); yindY = xindY(end);
%   iaX(JOB) = xindX;
%   iaY(JOB) = yindY;
% end

tileXcenter = X(JOB);
tileYcenter = Y(JOB);
indX = find(rlon73 < tileXcenter); indX = indX(end);
indY = find(rlat65 < tileYcenter); indY = indY(end);

junk = read_netcdf_lls('/asl/models/era5_avg//2022/2022-01_sfc.nc');
[era5_Y,era5_X] = meshgrid(junk.latitude,junk.longitude);
era5_X = wrapTo180(era5_X);
[salti, landfrac] = usgs_deg10_dem(era5_Y,era5_X);

usethese = find(era5_X >= rlon73(indX) & era5_X < rlon73(indX+1) & era5_Y >= rlat65(indY) & era5_Y < rlat65(indY+1));
whos usethese

iDorA = -1;  %% asc
iDorA = +1;  %% desc

N = 29; N1 = 1;
firstORend = 0;
iPrint = -1;

yy0 = simulateYear; mm0 = 01; dd0 = 15;
thedateS(1,:) = [yy0 mm0 dd0];

[yy1,mm1,dd1] = addNdays(yy0,mm0,dd0,N,firstORend,iPrint);
dd1 = 15;
thedateE(1,:) = [yy1 mm1 dd1];

for ii = 2 : 12
  yy0 = yy1; mm0 = mm1; dd0 = dd1;
  dd1 = 15;
  thedateS(ii,:) = [yy0 mm0 dd0];  

  [yy1,mm1,dd1] = addNdays(yy0,mm0,dd0,N,firstORend,iPrint);
  dd1 = 15;
  thedateE(ii,:) = [yy1 mm1 dd1];

  fprintf(1,'ii=%4i   Start %4i/%2i/%2i  End %4i/%2i/%2i \n',ii,thedateS(ii,:),thedateE(ii,:))
end  

[thedateS thedateE]
[yyM,mmM,ddM] = addNdays(thedateS(1,1),thedateS(1,2),thedateS(1,3),8,firstORend,iPrint);
rtimex = utc2taiSergio(yyM,mmM,15,0.00);

fprintf(1,'JOB = %4i spans %4i/%2i/%2i to %4i/%2i/%2i both ends inclusive \n',JOB,[thedateS(1,:) thedateE(12,:)]);

[h,ha,p,pa] = rtpread('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/summary_17years_all_lat_all_lon_2002_2019_palts_startSept2002_CLEAR.rtp');

%% now we need to get overpass times and solzen angles, just set scanang to 22 deg
load('/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_For_HowardObs_TimeSeries/asc_desc_solzen_time_412_64x72.mat');
monitor_memory_whos;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% doing simulateYear/01 to simulateYear/12

JOBX = 5 + 6; %% remember we assumed AIRS went up  in September, look for eg June

%% from comment, see /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_For_HowardObs_TimeSeries/driver_loop_get_asc_desc_solzen_time.m
if iDorA > 0
  rlon   = thedata.rlon_desc(:,:,JOBX);     rlon = rlon(:)';
  rlat   = thedata.rlat_desc(:,:,JOBX);     rlat = rlat(:)';
  solzen = thedata.solzen_desc(:,:,JOBX);   solzen = solzen(:)';
  satzen = thedata.satzen_desc(:,:,JOBX);   satzen = satzen(:)';
  hour   = thedata.hour_desc(:,:,JOBX);     hour = hour(:)';
  
  rtime = utc2taiSergio(yyM*ones(size(hour)),mmM*ones(size(hour)),15*ones(size(hour)),hour);
    [xmooY,xmooM,xmooD,xmooH] = tai2utcSergio(rtimex);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%

elseif iDorA < 0
  rlon   = thedata.rlon_asc(:,:,JOBX);     rlon = rlon(:)';
  rlat   = thedata.rlat_asc(:,:,JOBX);     rlat = rlat(:)';
  solzen = thedata.solzen_asc(:,:,JOBX);   solzen = solzen(:)';
  satzen = thedata.satzen_asc(:,:,JOBX);   satzen = satzen(:)';
  hour   = thedata.hour_asc(:,:,JOBX);     hour = hour(:)';
  
  rtime = utc2taiSergio(yyM*ones(size(hour)),mmM*ones(size(hour)),15*ones(size(hour)),hour);
    [xmooY,xmooM,xmooD,xmooH] = tai2utcSergio(rtimex);
end

rlon = wrapTo180(rlon);

figure(1); scatter_coast(rlon,rlat,50,solzen); title('solzen'); colormap jet; 
figure(2); scatter_coast(rlon,rlat,50,satzen); title('satzen'); colormap jet
figure(3); scatter_coast(rlon,rlat,50,hour);   title('hour'); colormap jet
%% figure(4); scatter_coast(rlon,rlat,50,bt1231); title('BT1231'); colormap jet
figure(5); scatter_coast(rlon,rlat,50,rlon-p.rlon); title('rlon-p.rlon'); colormap jet
figure(6); scatter_coast(rlon,rlat,50,rlat-p.rlat); title('rlat-p.rlat'); colormap jet
figure(4); simplemap(era5_Y(usethese),era5_X(usethese),landfrac(usethese)); title('Landfrac')

whos usethese

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% dirs used by /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_For_HowardObs_TimeSeries
if iDorA > 0
  fout = ['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Center12months/DESC/' num2str(simulateYear,'%04d') '/era5_full12months_tile_center_' num2str(JOB,'%03d') '.mat']; %%% NOTE THIS IS DESC
  fout = ['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Center12months/DESC/' num2str(simulateYear,'%04d') '/era5_full12months_tile_center_' num2str(JOB,'%04d') '.mat']; %%% NOTE THIS IS DESC
elseif iDorA < 0
  fout = ['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Center12months/ASC/' num2str(simulateYear,'%04d') '/era5_full12months_tile_center_' num2str(JOB,'%03d') '.mat']; %%% NOTE THIS IS DESC
  fout = ['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Center12months/ASC/' num2str(simulateYear,'%04d') '/era5_full12months_tile_center_' num2str(JOB,'%04d') '.mat']; %%% NOTE THIS IS DESC
else
  iDorA
  error('need iDorA = +/- 1')
end

iDo = +1;
if exist(fout)
  fprintf(1,'JOB %3i : avg era timeseries file %s already exists \n',JOB,fout);
  iDo = -1;
end

%iDo = +1;  %% forgot to add wet bulb temp and RH
if iDorA > 0
  avgtime = squeeze(thedata.rtime_desc(36,32,:));
else
  avgtime = squeeze(thedata.rtime_asc(36,32,:));
end
[yyy,mmm,ddd,hhh] = tai2utcSergio(avgtime);

if iDo > 0

  for mmx = 1 : 12
    mootime = utc2taiSergio(simulateYear,mmx,15,0);
    JOBX = find(avgtime >= mootime,1);

    if iDorA > 0
      rlon   = thedata.rlon_desc(:,:,JOBX);     rlon = rlon(:)';
      rlat   = thedata.rlat_desc(:,:,JOBX);     rlat = rlat(:)';
      solzen = thedata.solzen_desc(:,:,JOBX);   solzen = solzen(:)';
      satzen = thedata.satzen_desc(:,:,JOBX);   satzen = satzen(:)';
      hour   = thedata.hour_desc(:,:,JOBX);     hour = hour(:)';
    else
      rlon   = thedata.rlon_asc(:,:,JOBX);      rlon = rlon(:)';
      rlat   = thedata.rlat_asc(:,:,JOBX);      rlat = rlat(:)';
      solzen = thedata.solzen_asc(:,:,JOBX);    solzen = solzen(:)';
      satzen = thedata.satzen_asc(:,:,JOBX);    satzen = satzen(:)';
      hour   = thedata.hour_asc(:,:,JOBX);      hour = hour(:)';
    end

    xpnew_ip.rlon = era5_X(usethese)';
    xpnew_ip.rlat = era5_Y(usethese)';
    xpnew_ip.rlon(xpnew_ip.rlon > +177.5) = +177.5;
    xpnew_ip.rlon(xpnew_ip.rlon < -177.5) = -177.5;
  
    %xpnew_ip.rtime = interpn(rlon,rlat,rtime,xpnew_ip.rlon,xpnew_ip.rlat);
    xpnew_ip.rtime  = griddata(rlon,rlat,hour,wrapTo180(double(xpnew_ip.rlon)),double(xpnew_ip.rlat))*3600 + mootime;
    xpnew_ip.solzen = griddata(rlon,rlat,solzen,wrapTo180(double(xpnew_ip.rlon)),double(xpnew_ip.rlat));
    xpnew_ip.satzen = griddata(rlon,rlat,satzen,wrapTo180(double(xpnew_ip.rlon)),double(xpnew_ip.rlat));
  
    bad = find(isnan(xpnew_ip.rtime) | isnan(xpnew_ip.solzen) | isnan(xpnew_ip.satzen)); 
    good = find(isfinite(xpnew_ip.rtime) & isfinite(xpnew_ip.solzen) & isfinite(xpnew_ip.satzen)); 
    if length(bad) > 0
      xpnew_ip.rtime(bad) = nanmean(xpnew_ip.rtime(good));
      xpnew_ip.solzen(bad) = nanmean(xpnew_ip.solzen(good));
      xpnew_ip.satzen(bad) = nanmean(xpnew_ip.satzen(good));
    end

    if mmx == 1
      pnew_ip = xpnew_ip;
    else
      pnew_ip.rlon = [pnew_ip.rlon xpnew_ip.rlon];
      pnew_ip.rlat = [pnew_ip.rlat xpnew_ip.rlat];
      pnew_ip.rtime = [pnew_ip.rtime xpnew_ip.rtime];
      pnew_ip.solzen = [pnew_ip.solzen xpnew_ip.solzen];
      pnew_ip.satzen = [pnew_ip.satzen xpnew_ip.satzen];
    end

  end

  [junkyy,junkmm,junkdd,junkhh] = tai2utcSergio(pnew_ip.rtime);
  localtime = junkhh + pnew_ip.rlon/15;
  localtime = utc2localtime(junkhh,pnew_ip.rlon);  
  figure(5); clf; scatter_coast(pnew_ip.rlon,pnew_ip.rlat,100,localtime); xlabel('lon'); ylabel('lat'); title('local solar time'); pause(0.1);

  clrfields = {'SP','SKT','10U','10V','TCC','CI','T','Q','O3'};
  cldfields = {'SP','SKT','10U','10V','TCC','CI','T','Q','O3',...
                   'CC','CIWC','CLWC'};

  pnew_ip0 = pnew_ip;
  hnew_ip = struct;
    hnew_ip.ptype = 0;
  [pnew_ip,hnew_ip] = fill_era5_monthly(pnew_ip,hnew_ip);

  [xyy,xmm,xdd,xhh] = tai2utcSergio(pnew_ip.rtime);        % <<<<<<<<<<<<<<<<<<<<<<<< for SdSM old time
  time_so_far = (xyy-2000) + ((xmm-1)+1)/12;
  co2ppm = 368 + 2.077*time_so_far;  %% 395.6933
  pnew_ip.co2ppm = co2ppm;
  fprintf(1,'CLIMATOLOGY co2ppm for FIRST %4i/%2i/%2i = %8.6f ppmv\n',xyy(1),xmm(1),xdd(1),pnew_ip.co2ppm(1));
  fprintf(1,'CLIMATOLOGY co2ppm for LAST  %4i/%2i/%2i = %8.6f ppmv\n',xyy(end),xmm(end),xdd(end),pnew_ip.co2ppm(end));
  figure(1); scatter_coast(pnew_ip.rlon,pnew_ip.rlat,50,pnew_ip.co2ppm); title('co2ppm')
  figure(2); scatter_coast(pnew_ip.rlon,pnew_ip.rlat,50,pnew_ip.stemp);  title('stemp')
  pause(0.1)

  fip = mktemp('fx.ip.rtp');
  fop = mktemp('fx.op.rtp');
  frp = mktemp('fx.rp.rtp');

  run_sarta.klayers_code = '/asl/packages/klayersV205/BinV201/klayers_airs';
  run_sarta.clear = +1;
  run_sarta.cloud = +1;
  run_sarta.cumsum = 9999;  %% larrabee likes this, puts clouds high so does well for DCC
  run_sarta.cumsum = -1;    %% this is "closer" to MRO but since cliuds are at centroid, does not do too well with DCC
  code0 = '/asl/packages/sartaV108/BinV201/sarta_apr08_m140_iceaggr_waterdrop_desertdust_slabcloud_hg3_wcon_nte';
  code1 = '/home/sergio/SARTA_CLOUDY/BinV201/sarta_apr08_m140x_iceGHMbaum_waterdrop_desertdust_slabcloud_hg3';
  code1 = '/home/sergio/SARTA_CLOUDY/BinV201/xsarta_apr08_m140_iceGHMbaum_waterdrop_desertdust_slabcloud_hg3';
  code1 = '/home/chepplew/gitLib/sarta/bin/airs_l1c_2834_cloudy_may19_prod_v3';
  run_sarta.sartaclear_code = code1;
  run_sarta.sartacloud_code = code1;
  run_sarta.co2ppm = co2ppm;

  pa = {{'profiles','rtime','seconds since 1993'}};
  ha = {{'header','hdf file','era_month'}};

  [pnew_ip.salti, pnew_ip.landfrac] = usgs_deg10_dem(pnew_ip.rlat,pnew_ip.rlon);

  addpath /home/sergio/MATLABCODE/matlib/rtp_prod2/emis/
  addpath /home/sergio/MATLABCODE/matlib/rtp_prod2/util/

  pnew_ip.rlon = wrapTo180(pnew_ip.rlon);
  [pnew_ip,pa] = rtp_add_emis(pnew_ip,pa);

  haha = load('/home/sergio/MATLABCODE/CRODGERS_FAST_CLOUD/h2645structure.mat');
  hnew_ip.nchan = haha.h.nchan;
  hnew_ip.ichan = haha.h.ichan;
  hnew_ip.vchan = haha.h.vchan;

  [p2] = driver_sarta_cloud_rtp(hnew_ip,ha,pnew_ip,pa,run_sarta);
  %pnew_ip.rcalc = p2.rcalc;
  %pnew_ip.sarta_rclearcalc = p2.sarta_rclearcalc;

  rtpwrite(fip,hnew_ip,ha,p2,pa)
  klayerser = ['!' run_sarta.klayers_code ' fin=' fip ' fout=' fop ' >& ugh']; eval(klayerser);
  [hnew_op,ha2,pnew_op,pa2] = rtpread(fop);
  pnew_op.rcalc            = p2.rcalc;
  pnew_op.sarta_rclearcalc = p2.sarta_rclearcalc;
  pnew_op.mmw = mmwater_rtp(hnew_op,pnew_op);

  figure(1); dbt = 200:1:320; plot(dbt,hist(rad2bt(1231,p2.rcalc(1520,:)),dbt)); title('hist(BT1231 cld)')
  figure(2); scatter_coast(p2.rlon,p2.rlat,50,rad2bt(1231,p2.rcalc(1520,:))); title('BT1231 cld')
  figure(3); scatter_coast(p2.rlon,p2.rlat,50,rad2bt(1231,p2.sarta_rclearcalc(1520,:))); title('BT1231 clr')

  %%%%%%%%%%%%%%%%%%%%%%%%%
  quants = [0 0.03 0.05 0.1 0.2 0.5 0.8 0.9 0.95 0.97 1.0];
  booQ = quantile(rad2bt(1231,p2.rcalc(1520,:)),quants)
  find_avg_Qprofiles

  globalavg  = convert_rtp_to_cloudOD(hnew_op,globalavg);
  globalQavg = convert_rtp_to_cloudOD(hnew_op,globalQavg);
  %%%%%%%%%%%%%%%%%%%%%%%%%
  find_avg_Qplots
  %%%%%%%%%%%%%%%%%%%%%%%%%
  quick_jacs
  %%%%%%%%%%%%%%%%%%%%%%%%%

  rmer = ['!/bin/rm '  fip ' ' fop ' ' frp]; eval(rmer);

  [hnew_ip,ha,pnew_ip,pa] = rtptrim_sartacloud(hnew_ip,ha,p2,pa);

  [Tw,Tw1km,Tdew,WBGT,RH,RH1km,colwater,TwSurf,RHSurf,TdewSurf] = layeramt2RH_wet_bulb_dew_point_temperature(hnew_op,pnew_op);
  pnew_op.Tw   = Tw;
  pnew_op.Tdew = Tdew;
  pnew_op.RH   = RH;
  pnew_op.mmw  = colwater;
  pnew_op.TwSurf   = TwSurf;
  pnew_op.TdewSurf = TdewSurf;
  pnew_op.RHSurf   = RHSurf;

  thedateS = thedateS(1,:);
  comment = 'see /home/sergio/MATLABCODE/RTPMAKE/CLUST_RTPMAKE/CLUSTMAKE_ERA5ERA5/driver_whole_monthly_ERA5_test.m';
  saver = ['save ' fout ' comment        hnew_ip ha pnew_ip pa     hnew_op ha2 pnew_op pa2 thedateS global* quants '];
  eval(saver)

  %figure(3); scatter_coast(pnew_ip.rlon,pnew_ip.rlat,50,rad2bt(1231,pnew_ip.robs1(1520,:)))
  figure(14); scatter_coast(pnew_ip.rlon,pnew_ip.rlat,50,rad2bt(1231,pnew_ip.rcalc(1520,:)));            title('allsky calc');
  figure(15); scatter_coast(pnew_ip.rlon,pnew_ip.rlat,50,rad2bt(1231,pnew_ip.sarta_rclearcalc(1520,:))); title('clrsky calc');
end

disp('now do missinglist_whole_12month_ERA5.m  and then cluster_driver_put_together_globalavg_profiles.m and then master_driver_put_together_globalavg_profiles')
stoptime = toc
fprintf(1,'took %8.6f minutes to run \n',(stoptime-starttime)/60);
