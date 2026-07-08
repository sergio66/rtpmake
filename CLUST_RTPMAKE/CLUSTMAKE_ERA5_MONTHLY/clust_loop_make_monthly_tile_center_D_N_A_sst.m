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

%% 12 per year so about 240 in 20 years .....
JOB = str2num(getenv('SLURM_ARRAY_TASK_ID'));
if length(JOB) == 0
  JOB = 240;
  JOB = 120;
  JOB = 001;
  JOB = 169;
end


system_slurm_stats  %% does rng('shuffle')

iDorA = +1;   %% desc, use 1.30 am for all
iDorA = +10;  %% desc, use 1.30 pm for all, random pt about tile center

iDorA = -1;   %% asc,  use 1.30 pm for all
iDorA = -10;  %% asc,  use 1.30 pm for all, random pt about tile center

iDorA =  100; %% asc+desc, random pt about tile center
iDorA =  0;   %% asc+desc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% dirs used by /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_For_HowardObs_TimeSeries
if iDorA == +1
  fout = ['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Center/DAY_NIGHT_BOTH/NIGHT/era5_tile_center_DmonthlySST_' num2str(JOB,'%03d') '.mat']; %%% NOTE THIS IS DESC
elseif iDorA == +10
  fout = ['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Center/DAY_NIGHT_BOTH/NIGHT/randompt_era5_tile_center_DmonthlySST_' num2str(JOB,'%03d') '.mat']; %%% NOTE THIS IS DESC

elseif iDorA == -1
  fout = ['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Center/DAY_NIGHT_BOTH/DAY/era5_tile_center_NmonthlySST_' num2str(JOB,'%03d') '.mat']; %%% NOTE THIS IS ASC
elseif iDorA == -10
  fout = ['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Center/DAY_NIGHT_BOTH/DAY/randomptera5_tile_center_NmonthlySST_' num2str(JOB,'%03d') '.mat']; %%% NOTE THIS IS ASC

elseif iDorA == 0
  fout = ['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Center/DAY_NIGHT_BOTH/BOTH/era5_tile_center_BmonthlySST_' num2str(JOB,'%03d') '.mat']; %%% NOTE THIS IS ASC and DESC
elseif iDorA == 100
  fout = ['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Center/DAY_NIGHT_BOTH/BOTH/randomptera5_tile_center_BmonthlySST_' num2str(JOB,'%03d') '.mat']; %%% NOTE THIS IS ASC and DESC
else
  iDorA
  error('need iDorA = +/- 1,0')
end

%% /home/sergio/MATLABCODE/check_all_jobs_done('/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Center/DAY_NIGHT_BOTH/BOTH/randomptera5_tile_center_BmonthlySST_',240,'.mat')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[h,ha,p,pa] = rtpread('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/summary_17years_all_lat_all_lon_2002_2019_palts_startSept2002_CLEAR.rtp');

get_dates_loop_make_monthly2m_tile_center_asc_or_desc

iDo2m = -1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iDo = +1;
if exist(fout)
  fprintf(1,'JOB %3i : avg era5 timeseries file %s already exists \n',JOB,fout);
  iDo = -1;
  error('output file exists')
end

if iDorA == 1 | iDorA == +10
  iLoop = 4; iaOff = [-4.5 -1.5 +1.5 +4.5];  %% 1.30 am +/- these hours
elseif iDorA == -1 | iDorA == -10
  iLoop = 4; iaOff = [-4.5 -1.5 +1.5 +4.5];  %% 1.30 pm +/- these hours
elseif iDorA == 0 | iDorA == 100
  iLoop = 8; iaOff = [-10.5 -7.5 -4.5 -1.5 +1.5 +4.5 +7.5 +10.5];  %% 1.30 pm +/- these hourss
else
  error('iDorA')
end

%iDo = +1;  %% forgot to add wet bulb temp and RH
if iDo > 0
  % now fill in the fields see cloud_set_defaults_run_maker.m
  for iL = 1 : iLoop

    xpnew_ip = p;
    hnew_ip = h;

    hnew_ip = rmfield(hnew_ip,'ngas');
    hnew_ip = rmfield(hnew_ip,'glist');
    hnew_ip = rmfield(hnew_ip,'gunit');
  
    xpnew_ip.solzen = solzen;
    xpnew_ip.satzen = satzen;
    xpnew_ip.scanang = saconv(p.satzen,p.zobs);  

    xpnew_ip.rtime0  = rtime;
    xpnew_ip.rtime   = rtime + iaOff(iL)*3600;
    fprintf(1,'iL = %2i of %2i   avg rtime0 = %8.6e  avg rtime = %8.6e  avg diffrtime = rtime0-rtime = %8.6e \n',iL,iLoop,mean(xpnew_ip.rtime0),mean(xpnew_ip.rtime),mean(xpnew_ip.rtime0)-mean(xpnew_ip.rtime))

    if abs(iDorA) == 10 | iDorA == 100
%{
%%% https://en.wikipedia.org/wiki/Irwin%E2%80%93Hall_distribution : sum of N uniform distr over [a,b] is an Irwin-Hall distribution

rng('shuffle'); %% done by system_slurm_stats
%rng('default'); %% keeps things reproducible

q1 = (rand(1,4608)-0.5)*2.5*2;
q2 = (rand(1,4608)-0.5)*2.5*2;
q3 = (rand(1,4608)-0.5)*2.5*2;
q4 = (rand(1,4608)-0.5)*2.5*2;
q5 = (rand(1,4608)-0.5)*2.5*2;
q6 = (rand(1,4608)-0.5)*2.5*2;
q7 = (rand(1,4608)-0.5)*2.5*2;
q8 = (rand(1,4608)-0.5)*2.5*2;
hist(q1)
hist(q2)
hist(q3)
hist(q4)
hist((q1+q2+q3+q4+q5+q6+q7+q8)/8)  %% IS A NORMAL DISTR if shuffle
%}

      xpnew_ip.rlon0 = rlon;
      xpnew_ip.rlat0 = rlat;
      %%% https://en.wikipedia.org/wiki/Irwin%E2%80%93Hall_distribution : sum of N uniform distr over [a,b] is an Irwin-Hall distribution
      xpnew_ip.rlon  = rlon + (rand(size(rlon))-0.5)*2.50*2;  %% since all  tiles are 5 deg wide, should get nice uniform distribution                  hist(xpnew_ip.rlon0 - pnew_ip.rlon)
      xpnew_ip.rlat  = rlat + (rand(size(rlat))-0.5)*1.50*2;  %% since most tiles are 3 deg wide, some 5 deg wide, will get mostly uniform distribution hist(xpnew_ip.rlat0 - pnew_ip.rlat) but with noticeable tails
    end

    xpnew_ip = rmfield(xpnew_ip,'stemp');
    xpnew_ip = rmfield(xpnew_ip,'ptemp');
    xpnew_ip = rmfield(xpnew_ip,'plevs');
    xpnew_ip = rmfield(xpnew_ip,'palts');
    xpnew_ip = rmfield(xpnew_ip,'mmw');
    xpnew_ip = rmfield(xpnew_ip,'gas_1');
    xpnew_ip = rmfield(xpnew_ip,'gas_2');
    xpnew_ip = rmfield(xpnew_ip,'gas_3');
    xpnew_ip = rmfield(xpnew_ip,'gas_4');
    xpnew_ip = rmfield(xpnew_ip,'gas_5');
    xpnew_ip = rmfield(xpnew_ip,'gas_6');
    xpnew_ip = rmfield(xpnew_ip,'gas_9');
    xpnew_ip = rmfield(xpnew_ip,'gas_12');
  
    xpnew_ip = rmfield(xpnew_ip,'cprtop');
    xpnew_ip = rmfield(xpnew_ip,'cprbot');
    xpnew_ip = rmfield(xpnew_ip,'cfrac');
    xpnew_ip = rmfield(xpnew_ip,'cngwat');
    xpnew_ip = rmfield(xpnew_ip,'cpsize');
    xpnew_ip = rmfield(xpnew_ip,'ctype');
    xpnew_ip = rmfield(xpnew_ip,'cprtop2');
    xpnew_ip = rmfield(xpnew_ip,'cprbot2');
    xpnew_ip = rmfield(xpnew_ip,'cfrac2');
    xpnew_ip = rmfield(xpnew_ip,'cngwat2');
    xpnew_ip = rmfield(xpnew_ip,'cpsize2');
    xpnew_ip = rmfield(xpnew_ip,'ctype2');
    xpnew_ip = rmfield(xpnew_ip,'cfrac12');
  
    xpnew_ip = rmfield(xpnew_ip,'rcalc');
    xpnew_ip = rmfield(xpnew_ip,'sarta_rclearcalc');
  
    clrfields = {'SP','SKT','10U','10V','TCC','CI','T','Q','O3'};
    cldfields = {'SP','SKT','10U','10V','TCC','CI','T','Q','O3',...
                     'CC','CIWC','CLWC'};
  
    xpnew_ip0 = xpnew_ip;
    [xpnew_ip,hnew_ip] = fill_era5_monthly(xpnew_ip,hnew_ip);
    
    if (~isfield(xpnew_ip,'rh_2m') | ~isfield(xpnew_ip,'ta_2m') | ~isfield(xpnew_ip,'td_2m')) & iDo2m > 0
      [xpnew_ip,hnew_ip] = fill_era5_monthly2m_mat(xpnew_ip,hnew_ip,yyM,mmM);
      addpath /home/sergio/MATLABCODE/PLOTTER/
      addpath /home/sergio/MATLABCODE/COLORMAP
      addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS
      xpnew_ip.rh2m =  airtemp_dewpointtemp_2_RH(xpnew_ip.t2m,xpnew_ip.d2m,2);
      scatter_coast(xpnew_ip.rlon,xpnew_ip.rlat,50,xpnew_ip.stemp-xpnew_ip.skt_2m); colormap(usa2); caxis([-15 +15]); title('Stemp - MonthlyAvg SKT')
      scatter_coast(xpnew_ip.rlon,xpnew_ip.rlat,50,xpnew_ip.stemp-xpnew_ip.t2m);    colormap(usa2); caxis([-15 +15]); title('Stemp - 2m AirTemp')
      scatter_coast(xpnew_ip.rlon,xpnew_ip.rlat,50,xpnew_ip.stemp-xpnew_ip.d2m);    colormap(usa2); caxis([-15 +15]); title('Stemp - 2m DewPointTemp')
      scatter_coast(xpnew_ip.rlon,xpnew_ip.rlat,50,xpnew_ip.rh2m);                 colormap(jet); title('2m RH')
    end
  
    [xyy,xmm,xdd,xhh] = tai2utcSergio(xpnew_ip.rtime);        % <<<<<<<<<<<<<<<<<<<<<<<< for SdSM old time
    time_so_far = (xyy-2000) + ((xmm-1)+1)/12;
    co2ppm = 368 + 2.077*time_so_far;  %% 395.6933
    xpnew_ip.co2ppm = co2ppm;
    fprintf(1,'CLIMATOLOGY co2ppm for FIRST %4i/%2i/%2i = %8.6f ppmv\n',xyy(1),xmm(1),xdd(1),xpnew_ip.co2ppm(1));
    fprintf(1,'CLIMATOLOGY co2ppm for LAST  %4i/%2i/%2i = %8.6f ppmv\n',xyy(end),xmm(end),xdd(end),xpnew_ip.co2ppm(end));
    figure(1); scatter_coast(xpnew_ip.rlon,xpnew_ip.rlat,50,xpnew_ip.co2ppm);
    figure(2); scatter_coast(xpnew_ip.rlon,xpnew_ip.rlat,50,xpnew_ip.stemp); title('stemp')
    pause(0.1)

   junker = ['xp' num2str(iL) ' = xpnew_ip;']; eval(junker);
  end

  whos xp*
  if iDorA == 0 | iDorA == 100
    pnew_ip = average_Nprofiles(iLoop,xp1,xp2,xp3,xp4,xp5,xp6,xp7,xp8);
  else
    pnew_ip = average_Nprofiles(iLoop,xp1,xp2,xp3,xp4);
  end

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
 
  if iDo2m > 0
    pnew_op.d2m  = pnew_ip.d2m;
    pnew_op.t2m  = pnew_ip.t2m;
    pnew_op.rh2m = pnew_ip.rh2m;
    if isfield(pnew_ip,'skt_2m')
      pnew_op.skt_2m = pnew_ip.skt_2m; %% BIZARRE but remember that this comes from montly averages over all days in a month, see clust_make_monthlyavg_2m.m
    end
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
  if iDo2m > 0
    figure(6); scatter_coast(pnew_ip.rlon,pnew_ip.rlat,50,pnew_op.RHSurf-pnew_op.rh2m); caxis([-20 +20]); colormap(usa2); title('RH : mycalc-ERA5sfc = RHsurf-rh2m')
    if isfield(pnew_ip,'skt_2m')
      figure(7); scatter_coast(pnew_ip.rlon,pnew_ip.rlat,50,pnew_op.stemp-pnew_op.skt_2m); caxis([-0.1 +0.1]); colormap(usa2); title('Seeing effects of daily avg on stemp')
    end
  end

end
