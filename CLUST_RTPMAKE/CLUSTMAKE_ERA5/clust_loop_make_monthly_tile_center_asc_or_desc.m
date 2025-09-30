addpath /asl/matlib/rtptools/
addpath /asl/matlib/aslutil
addpath /asl/matlib/h4tools
addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/TIME
addpath /home/sergio/MATLABCODE/PLOTTER
addpath ../GRIB
addpath /home/sergio/MATLABCODE/matlib/clouds/sarta
addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS
addpath /home/sergio/MATLABCODE/COLORMAP

%% run with                  sbatch -p cpu2021 --exclude=  --array=220-240 sergio_matlab_jobB.sbatch 4
%% run with                  sbatch -p cpu2021 --exclude=  --array=241-264 sergio_matlab_jobB.sbatch 4
%% check success using eg    iaFound = check_all_jobs_done('/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Center/DESC_WithOLR/randomptera5_tile_center_monthly_',240,'.mat');  %% 20 years 2002/09-2022/08
%% check success using eg    iaFound = check_all_jobs_done('/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Center/DESC_WithOLR/randomptera5_tile_center_monthly_',264,'.mat');  %% 22 years 2002/09-2024/08
%% check success using eg    iaFound = check_all_jobs_done('/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Center/ASC_WithOLR/era5_tile_center_monthly_',264,'.mat');           %% 22 years 2002/09-2024/08
%% check success using eg    iaFound = check_all_jobs_done('/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Center/DESC_WithOLR/era5_tile_center_monthly_',264,'.mat');          %% 22 years 2002/09-2024/08

%% one per month, 19 years of AIRS data so 19x12 = 228 sets of data
%% one per month, 20 years of AIRS data so 20x12 = 240 sets of data
%% one per month, 21 years of AIRS data so 20x12 = 252 sets of data
%% one per month, 23 years of AIRS data so 20x12 = 276 sets of data
JOB = str2num(getenv('SLURM_ARRAY_TASK_ID'));
if length(JOB) == 0
  disp('no input JOB; being set to nYear * 12')
  JOB = 276;
end

%{
liststr = '/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Center/DESC_WithOLR/era5_tile_center_monthly_';
iaFound = check_all_jobs_done(liststr,252);
liststr = '/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Center/ASC_WithOLR/era5_tile_center_monthly_';
iaFound = check_all_jobs_done(liststr,252);
%}

%JOB = 235
%JOB = 009

system_slurm_stats

%%% iDorA = -90;  %% asc,  use 1.30 pm for all, use hottest 10% CANNOT DO, see clust_loop_make_monthly_tile_273points.m instead
%%% iDorA = +90;  %% desc, use 1.30 am for all, use hottest 10% CANNOT DO, see clust_loop_make_monthly_tile_273points.m instead

iDorA = -10;  %% asc,  use 1.30 pm for all, random pt about tile center
iDorA = +10;  %% desc, use 1.30 am for all, random pt about tile center  DONE THIS
iDorA = +1;   %% desc, use 1.30 am for all, tile center
iDorA = -1;   %% asc,  use 1.30 pm for all, tile center

iDo2m = -1;
iDo2m = +1;

[h,ha,p,pa] = rtpread('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/summary_17years_all_lat_all_lon_2002_2019_palts_startSept2002_CLEAR.rtp');

iNewOrOld = +1;
if iNewOrOld < 0
  get_dates_loop_make_monthly2m_tile_center_asc_or_desc_16days   %% this goes in 16 day intervals, kinda silly
else  
  get_dates_loop_make_monthly2m_tile_center_asc_or_desc_yy_mm    %% this goes in yy/mm intervals, more sensible
end

iOLR = +1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% dirs used by /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_For_HowardObs_TimeSeries
%% taki /asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/
%% chip /asl/s1/sergio/alldata/MakeAvgObsStats2002_2020_startSept2002_v3/

fdir0 = '/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/';
fdir0 = '/asl/s1/sergio/alldata/MakeAvgObsStats2002_2020_startSept2002_v3/';

if iDorA == +1
  fout = [fdir0 '/TimeSeries/ERA5/Tile_Center/DESC/era5_tile_center_monthly_' num2str(JOB,'%03d') '.mat']; %%% NOTE THIS IS DESC
  if iOLR > 0
    fout = [fdir0 '/TimeSeries/ERA5/Tile_Center/DESC_WithOLR/era5_tile_center_monthly_' num2str(JOB,'%03d') '.mat']; %%% NOTE THIS IS DESC
  end
elseif iDorA == -1
  fout = [fdir0 '/TimeSeries/ERA5/Tile_Center/ASC/era5_tile_center_monthly_' num2str(JOB,'%03d') '.mat']; %%% NOTE THIS IS DESC
  if iOLR > 0
    fout = [fdir0 '/TimeSeries/ERA5/Tile_Center/ASC_WithOLR/era5_tile_center_monthly_' num2str(JOB,'%03d') '.mat']; %%% NOTE THIS IS DESC
  end
elseif iDorA == +10
  fout = [fdir0 '/TimeSeries/ERA5/Tile_Center/DESC/randomptera5_tile_center_monthly_' num2str(JOB,'%03d') '.mat']; %%% NOTE THIS IS DESC
  if iOLR > 0
    fout = [fdir0 '/TimeSeries/ERA5/Tile_Center/DESC_WithOLR/randomptera5_tile_center_monthly_' num2str(JOB,'%03d') '.mat']; %%% NOTE THIS IS DESC      DONE THIS DONE THIS
  end
elseif iDorA == -10
  fout = [fdir0 '/TimeSeries/ERA5/Tile_Center/ASC/randomptera5_tile_center_monthly_' num2str(JOB,'%03d') '.mat']; %%% NOTE THIS IS DESC
  if iOLR > 0
    fout = [fdir0 '/TimeSeries/ERA5/Tile_Center/ASC_WithOLR/randomptera5_tile_center_monthly_' num2str(JOB,'%03d') '.mat']; %%% NOTE THIS IS DESC
  end
% elseif iDorA == +90
%   fout = [fdir0 '/TimeSeries/ERA5/Tile_Center/DESC/hottest_10percent_era5_tile_center_monthly_' num2str(JOB,'%03d') '.mat']; %%% NOTE THIS IS DESC
%   if iOLR > 0
%     fout = [fdir0 '/TimeSeries/ERA5/Tile_Center/DESC_WithOLR/hottest_10percent_era5_tile_center_monthly_' num2str(JOB,'%03d') '.mat']; %%% NOTE THIS IS DESC
%   end
% elseif iDorA == -90
%   fout = [fdir0 '/TimeSeries/ERA5/Tile_Center/ASC/hottest_10percent_era5_tile_center_monthly_' num2str(JOB,'%03d') '.mat']; %%% NOTE THIS IS DESC
%   if iOLR > 0
%     fout = [fdir0 '/TimeSeries/ERA5/Tile_Center/ASC_WithOLR/hottest_10percent_era5_tile_center_monthly_' num2str(JOB,'%03d') '.mat']; %%% NOTE THIS IS DESC
%   end
else
  iDorA
  error('need iDorA = +/- 1, +/- 10')
end

iDo = +1;
if exist(fout)
  fprintf(1,'JOB %3i : avg era timeseries file %s already exists \n',JOB,fout);
  iDo = -1;
  error('output file exists')
end

iDo = +1;  %% forgot to add RH and wet bulb temp, temp at 2m
if iDo > 0
  % now fill in the fields see cloud_set_defaults_run_maker.m
  pnew_ip = p;
  hnew_ip = h;
  hnew_ip = rmfield(hnew_ip,'ngas');
  hnew_ip = rmfield(hnew_ip,'glist');
  hnew_ip = rmfield(hnew_ip,'gunit');

  if abs(iDorA) == 10 | iDorA == 100
    pnew_ip.rlon0 = rlon;
    pnew_ip.rlat0 = rlat;
    pnew_ip.rlon  = rlon + (rand(size(rlon))-0.5)*2.50*2;  %% since all  tiles are 5 deg wide, should get nice uniform distribution                  hist(pnew_ip.rlon0 - pnew_ip.rlon)
    pnew_ip.rlat  = rlat + (rand(size(rlat))-0.5)*1.50*2;  %% since most tiles are 3 deg wide, some 5 deg wide, will get mostly uniform distribution hist(pnew_ip.rlat0 - pnew_ip.rlat) but with noticeable tails
  end

  pnew_ip.solzen = solzen;
  pnew_ip.satzen = satzen;
  pnew_ip.scanang = saconv(p.satzen,p.zobs);  
  pnew_ip.rtime   = rtime;
  [yywah,mmwah,ddwah,hhwah] = tai2utcSergio(rtime);
  [yywah,mmwah,ddwah,hhwah] = tai2utcSergio(mean(rtime));
  fprintf(1,' JOB = %3i   rtime = %10.6f  yywah/mmwah/ddwah = %4i/%2i/%2i \n',JOB,mean(rtime),yywah,mmwah,ddwah)
  
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
  cldfields = {'SP','SKT','10U','10V','TCC','CI','T','Q','O3','CC','CIWC','CLWC'};

  if iOLR > 0
    clrfields = {'SP','SKT','10U','10V','TCC','CI','T','Q','O3','OLR','OLRCS'};
    cldfields = {'SP','SKT','10U','10V','TCC','CI','T','Q','O3','CC','CIWC','CLWC','OLR','OLRCS'};
  end

  pnew_ip0 = pnew_ip;
  [pnew_ip,hnew_ip,~,iOLR] = fill_era5_monthly(pnew_ip,hnew_ip,[],iOLR);

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
  %code0 = '/asl/packages/sartaV108/BinV201/sarta_apr08_m140_iceaggr_waterdrop_desertdust_slabcloud_hg3_wcon_nte';
  %code1 = '/home/sergio/SARTA_CLOUDY/BinV201/sarta_apr08_m140x_iceGHMbaum_waterdrop_desertdust_slabcloud_hg3';
  %code1 = '/home/sergio/SARTA_CLOUDY/BinV201/xsarta_apr08_m140_iceGHMbaum_waterdrop_desertdust_slabcloud_hg3';
  %code1 = '/home/chepplew/gitLib/sarta/bin/airs_l1c_2834_cloudy_may19_prod_v3';
  code1 = '/home/sergio/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/JACvers/bin/jac_airs_l1c_2834_cloudy_jan25_H2020';
  run_sarta.sartaclear_code = code1;
  run_sarta.sartacloud_code = code1;
  run_sarta.co2ppm = co2ppm;

  if isfield(pnew_ip,'olr_clr')
    pnew_ip.olr_clr = abs(pnew_ip.olr_clr);
    pnew_ip.olr     = abs(pnew_ip.olr);
    pnew_ip.ilr_clr = abs(pnew_ip.ilr_clr);
    pnew_ip.ilr     = abs(pnew_ip.ilr);
  end

  [p2] = driver_sarta_cloud_rtp(hnew_ip,ha,pnew_ip,pa,run_sarta);
  %pnew_ip.rcalc = p2.rcalc;
  %pnew_ip.sarta_rclearcalc = p2.sarta_rclearcalc;

  rtpwrite(fip,hnew_ip,ha,p2,pa)
  klayerser = ['!' run_sarta.klayers_code ' fin=' fip ' fout=' fop ' >& ugh']; eval(klayerser);
  [hnew_op,ha2,pnew_op,pa2] = rtpread(fop);
  pnew_op.rcalc            = p2.rcalc;
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

  if isfield(pnew_ip,'olr_clr')
    pnew_op.olr_clr = pnew_ip.olr_clr;
    pnew_op.olr     = pnew_ip.olr;
    pnew_op.ilr_clr = pnew_ip.ilr_clr;
    pnew_op.ilr     = pnew_ip.ilr;
    figure(6); scatter_coast(pnew_ip.rlon,pnew_ip.rlat,50,pnew_ip.olr);     title('allsky OLR W/m2'); colormap jet
    figure(7); scatter_coast(pnew_ip.rlon,pnew_ip.rlat,50,pnew_ip.olr_clr); title('clrsky OLR W/m2'); colormap jet

    figure(8); scatter_coast(pnew_ip.rlon,pnew_ip.rlat,50,pnew_ip.ilr);     title('allsky ILR W/m2'); colormap jet
    figure(9); scatter_coast(pnew_ip.rlon,pnew_ip.rlat,50,pnew_ip.ilr_clr); title('clrsky ILR W/m2'); colormap jet

    figure(8); scatter_coast(pnew_ip.rlon,pnew_ip.rlat,50,pnew_ip.ilr_clr); title('clrsky ILR W/m2'); colormap jet
    figure(9); scatter_coast(pnew_ip.rlon,pnew_ip.rlat,50,-pnew_ip.ilr_clr + 5.67e-8 * pnew_op.stemp.^(4)); title('clrsky ILR W/m2 with STEMP^4'); colormap jet
      caxis([0 450]-25); %% see fig 24 in Wang, K., and R. E. Dickinson (2013), Rev. Geophys., 51, 150–185, doi:10.1002/rog.20009.
                         %% Global atmospheric downward longwave radiation at the surface from ground-based observations, satellite retrievals, and reanalyses, 

    bonk = 800:1:1200; 
    clear kabonk
    for iii = 1 : length(pnew_ip.stemp)
      kabonk(:,iii) = ttorad(bonk,pnew_ip.stemp(iii)); 
    end
    kabonk = sum(kabonk,1) * mean(diff(bonk))/1000;
    kabonk = kabonk * 3.0 .* (1-pnew_op.mmw/max(pnew_op.mmw)); %% make the remoeved ILR at window a function of column water; 
                                                               %% less if there is a lot of water (opque), more if there is little water (transperent)
    figure(17); scatter_coast(pnew_ip.rlon,pnew_ip.rlat,50,pnew_ip.stemp.^4 * 5.67e-8 - kabonk);     title('expected ILR W/m2'); colormap jet
      caxis([50 450])
      caxis([0 450]-25); %% see fig 24 in Wang, K., and R. E. Dickinson (2013), Rev. Geophys., 51, 150–185, doi:10.1002/rog.20009.
                         %% Global atmospheric downward longwave radiation at the surface from ground-based observations, satellite retrievals, and reanalyses, 
  end

  if isfield(pnew_ip,'d2m') & iDo2m > 0
    pnew_op.d2m = pnew_ip.d2m;
    pnew_op.t2m = pnew_ip.t2m;
    figure(10); scatter_coast(pnew_ip.rlon,pnew_ip.rlat,50,pnew_op.d2m);      title('ERA5   dew point 2m (K)'); colormap jet; caxis([200 320])
    figure(11); scatter_coast(pnew_ip.rlon,pnew_ip.rlat,50,pnew_op.t2m);      title('ERA5   air temp  2m (K)'); colormap jet; caxis([200 320])
    figure(12); scatter_coast(pnew_ip.rlon,pnew_ip.rlat,50,pnew_op.TdewSurf); title('Sergio dew point 2m (K)'); colormap jet; caxis([200 320])
    figure(13); scatter_coast(pnew_ip.rlon,pnew_ip.rlat,50,pnew_op.d2m-pnew_op.TdewSurf); title('(ERA5-Sergio) dew point 2m (K)'); colormap(usa2); caxis([-1 +1]*5)

    e2a = 6.1079*exp(17.269 * (pnew_op.d2m-273.13)./(237.3 + (pnew_op.d2m-273.13)));            %% from "Understanding variations in downwelling longwave radiation using Brutsaert’s equation"
                                                                              %% Earth Syst. Dynam., 14, 1363–1374, 2023 https://doi.org/10.5194/esd-14-1363-2023, eqn 6
    figure(14); scatter_coast(pnew_ip.rlon,pnew_ip.rlat,50,e2a);              title('ERA5   vapor pressure 2m (mb)'); colormap jet; 

    ecs = 1.24*(e2a./pnew_op.t2m).^(1/7);                                     %% from "Understanding variations in downwelling longwave radiation using Brutsaert’s equation"
    figure(15); scatter_coast(pnew_ip.rlon,pnew_ip.rlat,50,ecs);              title('Brutsaert Emissivity')

    Rld = 5.67e-8 * ecs .* pnew_op.t2m.^(4);                                  %% Earth Syst. Dynam., 14, 1363–1374, 2023 https://doi.org/10.5194/esd-14-1363-2023, eqn 1,2
    figure(16); scatter_coast(pnew_ip.rlon,pnew_ip.rlat,50,Rld);              title('Brutsaert ILR W/m2')
      caxis([50 450])
      caxis([0 450]-25); %% see fig 24 in Wang, K., and R. E. Dickinson (2013), Rev. Geophys., 51, 150–185, doi:10.1002/rog.20009.
                         %% Global atmospheric downward longwave radiation at the surface from ground-based observations, satellite retrievals, and reanalyses, 
  end

  thedateS = thedateS(JOB,:);
  comment = 'see /home/sergio/MATLABCODE/RTPMAKE/CLUST_RTPMAKE/CLUSTMAKE_ERA5ERA5/clust_loop_make_monthly_tile_center.m';
  saver = ['save ' fout ' comment        hnew_ip ha pnew_ip pa     hnew_op ha2 pnew_op pa2 thedateS'];

  eval(saver)

  %figure(3); scatter_coast(pnew_ip.rlon,pnew_ip.rlat,50,rad2bt(1231,pnew_ip.robs1(1520,:)))
  figure(4); scatter_coast(pnew_ip.rlon,pnew_ip.rlat,50,rad2bt(1231,pnew_ip.rcalc(1520,:)));            title('allsky calc');
  figure(5); scatter_coast(pnew_ip.rlon,pnew_ip.rlat,50,rad2bt(1231,pnew_ip.sarta_rclearcalc(1520,:))); title('clrsky calc');

  fprintf(1,'%5i profiles saved to %s \n',length(pnew_op.stemp),fout)
end
