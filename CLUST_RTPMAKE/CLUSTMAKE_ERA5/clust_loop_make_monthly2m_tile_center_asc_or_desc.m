addpath /asl/matlib/rtptools/
addpath /asl/matlib/aslutil
addpath /asl/matlib/h4tools
addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/TIME
addpath /home/sergio/MATLABCODE/PLOTTER
addpath ../GRIB
addpath /home/sergio/MATLABCODE/matlib/clouds/sarta
addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS

%% one per month, 19 years of AIRS data so 19x12 = 228 sets of data
%% one per month, 20 years of AIRS data so 19x12 = 240 sets of data
%% one per month, 21 years of AIRS data so 19x12 = 252 sets of data
%JOB = 06  %% 2003/02
%JOB = 14  %% 2003/08

%% 23 per year so about 460 in 20 years .....
JOB = str2num(getenv('SLURM_ARRAY_TASK_ID'));
if length(JOB) == 0
  JOB = 240;
  JOB = 120;
end

system_slurm_stats

iDorA = +1;  %% desc
iDorA = -1;  %% asc

[h,ha,p,pa] = rtpread('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/summary_17years_all_lat_all_lon_2002_2019_palts_startSept2002_CLEAR.rtp');

get_dates_loop_make_monthly2m_tile_center_asc_or_desc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% dirs used by /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_For_HowardObs_TimeSeries
if iDorA > 0
  fout = ['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Center/DESC/era5_tile_center_monthly2m_' num2str(JOB,'%03d') '.mat']; %%% NOTE THIS IS DESC
elseif iDorA < 0
  fout = ['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Center/ASC/era5_tile_center_monthly2m_' num2str(JOB,'%03d') '.mat']; %%% NOTE THIS IS ASC
else
  iDorA
  error('need iDorA = +/- 1')
end

iDo = +1;
if exist(fout)
  fprintf(1,'JOB %3i : avg era5 timeseries file %s already exists \n',JOB,fout);
  iDo = -1;
end

%iDo = +1;  %% forgot to add wet bulb temp and RH
if iDo > 0
  % now fill in the fields see cloud_set_defaults_run_maker.m
  pnew_ip = p;
  hnew_ip = h;
  hnew_ip = rmfield(hnew_ip,'ngas');
  hnew_ip = rmfield(hnew_ip,'glist');
  hnew_ip = rmfield(hnew_ip,'gunit');

  pnew_ip.solzen = solzen;
  pnew_ip.satzen = satzen;
  pnew_ip.scanang = saconv(p.satzen,p.zobs);  
  pnew_ip.rtime   = rtime;

  pnew_ip = rmfield(pnew_ip,'stemp');
  pnew_ip = rmfield(pnew_ip,'ptemp');
  pnew_ip = rmfield(pnew_ip,'plevs');
  pnew_ip = rmfield(pnew_ip,'palts');
  pnew_ip = rmfield(pnew_ip,'mmw');
  pnew_ip = rmfield(pnew_ip,'gas_1');
  pnew_ip = rmfield(pnew_ip,'gas_2');
  pnew_ip = rmfield(pnew_ip,'gas_3');
  pnew_ip = rmfield(pnew_ip,'gas_4');
  pnew_ip = rmfield(pnew_ip,'gas_5');
  pnew_ip = rmfield(pnew_ip,'gas_6');
  pnew_ip = rmfield(pnew_ip,'gas_9');
  pnew_ip = rmfield(pnew_ip,'gas_12');

  pnew_ip = rmfield(pnew_ip,'cprtop');
  pnew_ip = rmfield(pnew_ip,'cprbot');
  pnew_ip = rmfield(pnew_ip,'cfrac');
  pnew_ip = rmfield(pnew_ip,'cngwat');
  pnew_ip = rmfield(pnew_ip,'cpsize');
  pnew_ip = rmfield(pnew_ip,'ctype');
  pnew_ip = rmfield(pnew_ip,'cprtop2');
  pnew_ip = rmfield(pnew_ip,'cprbot2');
  pnew_ip = rmfield(pnew_ip,'cfrac2');
  pnew_ip = rmfield(pnew_ip,'cngwat2');
  pnew_ip = rmfield(pnew_ip,'cpsize2');
  pnew_ip = rmfield(pnew_ip,'ctype2');
  pnew_ip = rmfield(pnew_ip,'cfrac12');

  pnew_ip = rmfield(pnew_ip,'rcalc');
  pnew_ip = rmfield(pnew_ip,'sarta_rclearcalc');

  clrfields = {'SP','SKT','10U','10V','TCC','CI','T','Q','O3'};
  cldfields = {'SP','SKT','10U','10V','TCC','CI','T','Q','O3',...
                   'CC','CIWC','CLWC'};

  pnew_ip0 = pnew_ip;
  [pnew_ip,hnew_ip] = fill_era5_monthly(pnew_ip,hnew_ip);
  if ~isfield(pnew_ip,'rh_2m') | ~isfield(pnew_ip,'ta_2m') | ~isfield(pnew_ip,'td_2m')
    [pnew_ip,hnew_ip] = fill_era5_monthly2m_mat(pnew_ip,hnew_ip,yyM,mmM);
    addpath /home/sergio/MATLABCODE/PLOTTER/
    addpath /home/sergio/MATLABCODE/COLORMAP
    addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS
    pnew_ip.rh2m =  airtemp_dewpointtemp_2_RH(pnew_ip.t2m,pnew_ip.d2m,2);
    scatter_coast(pnew_ip.rlon,pnew_ip.rlat,50,pnew_ip.stemp-pnew_ip.skt_2m); colormap(usa2); caxis([-15 +15]); title('Stemp - MonthlyAvg SKT')
    scatter_coast(pnew_ip.rlon,pnew_ip.rlat,50,pnew_ip.stemp-pnew_ip.t2m);    colormap(usa2); caxis([-15 +15]); title('Stemp - 2m AirTemp')
    scatter_coast(pnew_ip.rlon,pnew_ip.rlat,50,pnew_ip.stemp-pnew_ip.d2m);    colormap(usa2); caxis([-15 +15]); title('Stemp - 2m DewPointTemp')
    scatter_coast(pnew_ip.rlon,pnew_ip.rlat,50,pnew_ip.rh2m);                 colormap(jet); title('2m RH')
  end

  [xyy,xmm,xdd,xhh] = tai2utcSergio(pnew_ip.rtime);        % <<<<<<<<<<<<<<<<<<<<<<<< for SdSM old time
  time_so_far = (xyy-2000) + ((xmm-1)+1)/12;
  co2ppm = 368 + 2.077*time_so_far;  %% 395.6933
  pnew_ip.co2ppm = co2ppm;
  fprintf(1,'CLIMATOLOGY co2ppm for FIRST %4i/%2i/%2i = %8.6f ppmv\n',xyy(1),xmm(1),xdd(1),pnew_ip.co2ppm(1));
  fprintf(1,'CLIMATOLOGY co2ppm for LAST  %4i/%2i/%2i = %8.6f ppmv\n',xyy(end),xmm(end),xdd(end),pnew_ip.co2ppm(end));
  figure(1); scatter_coast(pnew_ip.rlon,pnew_ip.rlat,50,pnew_ip.co2ppm);
  figure(2); scatter_coast(pnew_ip.rlon,pnew_ip.rlat,50,pnew_ip.stemp); title('stemp')
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

  [p2] = driver_sarta_cloud_rtp(hnew_ip,ha,pnew_ip,pa,run_sarta);
  %pnew_ip.rcalc = p2.rcalc;
  %pnew_ip.sarta_rclearcalc = p2.sarta_rclearcalc;

  rtpwrite(fip,hnew_ip,ha,p2,pa)
  klayerser = ['!' run_sarta.klayers_code ' fin=' fip ' fout=' fop ' >& ugh']; eval(klayerser);
  [hnew_op,ha2,pnew_op,pa2] = rtpread(fop);
  pnew_op.sarta_rclearcalc = p2.sarta_rclearcalc;
  pnew_op.mmw = mmwater_rtp(hnew_op,pnew_op);

  [hnew_ip,ha,pnew_ip,pa] = rtptrim_sartacloud(hnew_ip,ha,p2,pa);

  rmer = ['!/bin/rm '  fip ' ' fop ' ' frp]; eval(rmer);

  [Tw,Tw1km,Tdew,WBGT,RH,RH1km,colwater,TwSurf,RHSurf,TdewSurf] = layeramt2RH_wet_bulb_dew_point_temperature(hnew_op,pnew_op);
  pnew_op.Tw   = Tw;
  pnew_op.Tdew = Tdew;
  pnew_op.RH   = RH;
  pnew_op.mmw  = colwater;
  pnew_op.TwSurf   = TwSurf;
  pnew_op.TdewSurf = TdewSurf;
  pnew_op.RHSurf   = RHSurf;
  
  pnew_op.d2m  = pnew_ip.d2m;
  pnew_op.t2m  = pnew_ip.t2m;
  pnew_op.rh2m = pnew_ip.rh2m;
  if isfield(pnew_ip,'skt_2m')
    pnew_op.skt_2m = pnew_ip.skt_2m; %% BIZARRE but remember that this comes from montly averages over all days in a month, see clust_make_monthlyavg_2m.m
  end

  thedateS = thedateS(JOB,:);
  comment = 'see /home/sergio/MATLABCODE/RTPMAKE/CLUST_RTPMAKE/CLUSTMAKE_ERA5ERA5/clust_loop_make_monthly_tile_center.m';
  saver = ['save ' fout ' comment        hnew_ip ha pnew_ip pa     hnew_op ha2 pnew_op pa2 thedateS'];
  saver
  eval(saver)

  %figure(3); scatter_coast(pnew_ip.rlon,pnew_ip.rlat,50,rad2bt(1231,pnew_ip.robs1(1520,:)))
  figure(4); scatter_coast(pnew_ip.rlon,pnew_ip.rlat,50,rad2bt(1231,pnew_ip.rcalc(1520,:)));            title('allsky calc');
  figure(5); scatter_coast(pnew_ip.rlon,pnew_ip.rlat,50,rad2bt(1231,pnew_ip.sarta_rclearcalc(1520,:))); title('clrsky calc');

  figure(6); scatter_coast(pnew_ip.rlon,pnew_ip.rlat,50,pnew_op.RHSurf); caxis([0 120])
  figure(6); scatter_coast(pnew_ip.rlon,pnew_ip.rlat,50,pnew_op.RHSurf-pnew_op.rh2m); caxis([-20 +20]); colormap(usa2); title('RH : mycalc-ERA5sfc = RHsurf-rh2m')

  if isfield(pnew_ip,'skt_2m')
    figure(7); scatter_coast(pnew_ip.rlon,pnew_ip.rlat,50,pnew_op.stemp-pnew_op.skt_2m); caxis([-0.1 +0.1]); colormap(usa2); title('Seeing effects of daily avg on stemp')
  end

end
