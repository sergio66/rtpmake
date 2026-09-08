%{
MLS data in /asl/xfs3/mls
MLS data in /umbc/xfs3/strow/asl/mls
%}

mlsdata = '/asl/xfs3/mls/';                  vers = '004'; xyvers = 'v04';
mlsdata = '/umbc/xfs3/strow/asl/mls/V004/';  vers = '004'; xyvers = 'v04';

mlsdata = '/umbc/xfs3/strow/asl/mls/V005/';  vers = '005'; xyvers = 'v05';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% run this with sbatch -p cpu2024  --array=1-216 sergio_matlab_chip.sbatch 1          for mls 2004-2022
% run this with sbatch -p cpu2024  --array=1-252 sergio_matlab_chip.sbatch 1          for mls 2004-2025
% run this with sbatch -p high_mem --array=1-276 sergio_matlab_jobB.sbatch 1
%   brought data down 2004 to 2025
%   monthly steps since 2004 to 2025 = 21 x 12 = 252 so anything past 253->blah is not gonna be done

JOB = str2num(getenv('SLURM_ARRAY_TASK_ID'));
if length(JOB) == 0
  JOB = 172;
  JOB = 1;  
end  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dirout = '/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/';
dirout = '/home/sergio/git/oem_climate_jacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/DATAObsStats_StartSept2002_CORRECT_LatLon/';
%% data saved to eg /home/sergio/git/oem_climate_jacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/DATAObsStats_StartSept2002_CORRECT_LatLon//TimeSeries/MLS/Tile_Center/*.mat | wc -l

%{
liststr = '/home/sergio/git/oem_climate_jacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/DATAObsStats_StartSept2002_CORRECT_LatLon//TimeSeries/MLS/Tile_Center/mls_tile_center_monthly_timestep_';
check_all_jobs_done(liststr,252,'.mat');

jobs_not_done.sc : sbatch   --array=172,237,244-245,247,250 sergio_matlab_chip.sbatch
>> [yyuse([172 237 244 245 247 250]); mmuse([172 237 244 245 247 250])]'
        2018          12
        2024           5
        2024          12
        2025           1
        2025           3
        2025           6

slurm-199004_237.out:139:bad gas_1
slurm-199004_244.out:139:bad gas_1
slurm-199004_247.out:139:bad gas_1
slurm-199004_250.out:139:bad gas_1

see do_gridded_interpolant.m
slurm-199004_172.out:120:  junkT      1x0                 0  double
slurm-199004_172.out:130:  junkO      1x0                 0  double
slurm-199004_172.out:125:  junkW      1x0                 0  double

slurm-199004_245.out:120:  junkT      1x43              344  double
slurm-199004_245.out:130:  junkO      1x42              336  double
slurm-199004_245.out:125:  junkW      1x0                 0  double
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath0
addpath /home/sergio/git/matlabcode/matlibSergio/matlib/rtp_prod2_Aug11_2020/emis/

addpath ../COMMON_SETTINGS
set_path_to_execs

fip = mktempS('fx.ip.rtp');
fop = mktempS('fx.op.rtp');
frp = mktempS('fx.rp.rtp');

iNumYears = 20;
iNumYears = 23;

% see /umbc/rs/pi_sergio/WorkDirDec2025/oem_climate_code/AIRS_gridded_STM_May2021_trendsonlyCLR/read_fileMean17years.m
if iNumYears == 20
  fileMean17years = '/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/summary_20years_all_lat_all_lon_2002_2022_monthlyERA5.rp.rtp';
  fileMean17years = '/home/sergio/git/kcarta_gen/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/summary_23years_all_lat_all_lon_2002_2025_monthlyERA5.op.rtp';   %% this is YY = 23
  yySE = [2004 2021];  %% should do this
  yySE = [2004 2022];  %% but only have this data  
elseif iNumYears == 23
  fileMean17years = '/home/sergio/git/kcarta_gen/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/summary_23years_all_lat_all_lon_2002_2025_monthlyERA5.op.rtp';
  fileMean17years = '/home/sergio/git/kcarta_gen/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/summary_23years_all_lat_all_lon_2002_2025_monthlyERA5.op.rtp';   %% this is YY = 23
  yySE = [2004 2021];  %% should do this
  yySE = [2004 2025];  %% but only have this data    
end

[h,ha,p,pa] = rtpread(fileMean17years);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yyuse = [2004 2004 2004 2004];
yyuse = ones(1,4)*yySE(1); 
mmuse = [09   10   11   12  ];

%% now go from 2005 to 20XY
for yy = yySE(1)+1 : yySE(2)-1
  yjunk = yy*ones(1,12);
  mjunk = [1 2 3 4 5 6 7 8 9 10 11 12];
  yyuse = [yyuse yjunk];
  mmuse = [mmuse mjunk];
end
yy = 2020;    yjunk = yy*ones(1,8);
yy = yySE(2); yjunk = yy*ones(1,8);
           mjunk = [1 2 3 4 5 6 7 8];
  yyuse = [yyuse yjunk];
  mmuse = [mmuse mjunk];
[1:length(yyuse); yyuse; mmuse]';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
whos yyuse mmuse

%% just pull out the one for YY(ii) MM(ii) that you want
for ii = JOB
  [ii yyuse(ii) mmuse(ii)]  
  p.rtime = ones(size(p.rtime))*utc2taiSergio(yyuse(ii),mmuse(ii),15,12.0);

  fout = [dirout '/TimeSeries/MLS/Tile_Center/mls_tile_center_monthly_timestep_' num2str(ii,'%03d') '.mat']; 
  if exist(fout)
    fprintf(1,'%s already exists ... not saving \n',fout) 
    error('oioi')
  end

  [hL3,haL3,pL3,paL3] = rtp_airsL3_climatologyFAST_rtime(p.rtime,p.rlon,p.rlat,p.landfrac,-1);
  hL3.ngas = 2;
  hL3.glist = [1 3]';
  hL3.gunit = [20 12]';  %% this is native
  hL3.gunit = [10 10]';  %% change to ppmv so we can tack on MLS
  hL3.nchan = 2645;
  h2645 = load('/home/sergio/git/umbc_singlefootprint/h2645structure.mat');
  hL3.ichan = h2645.h.ichan;
  hL3.vchan = h2645.h.vchan;
  pL3 = rmfield(pL3,'gas_5');
  pL3 = rmfield(pL3,'gas_6');
  pL3.gas_1 = toppmv(pL3.plevs,pL3.ptemp,pL3.gas_1,18,20);
  pL3.gas_3 = toppmv(pL3.plevs,pL3.ptemp,pL3.gas_3,48,12);
  pL3.zobs  = p.zobs;
  pL3.satzen  = p.satzen;
  pL3.solzen  = p.solzen;
  pL3.scanang = p.scanang;
  pL3.salti  = p.salti;

  mls_T  = [mlsdata 'ML3MBT_' vers '/MLS-Aura_L3MB-Temperature_' xyvers '-23-c02_' num2str(yyuse(ii)) '.nc'];
  mls_W  = [mlsdata 'ML3MBH2O_' vers '/MLS-Aura_L3MB-H2O_' xyvers '-23-c02_' num2str(yyuse(ii)) '.nc'];
  mls_O3 = [mlsdata 'ML3MBO3_' vers '/MLS-Aura_L3MB-O3_' xyvers '-23-c02_' num2str(yyuse(ii)) '.nc'];

  mls_T  = [mlsdata 'ML3MBT_' vers '/MLS-Aura_L3MB-Temperature_' xyvers '*' num2str(yyuse(ii)) '.nc'];
    bonk = dir(mls_T);
    mls_T = [mlsdata 'ML3MBT_' vers '/' bonk.name];
  mls_W  = [mlsdata 'ML3MBH2O_' vers '/MLS-Aura_L3MB-H2O_' xyvers '*' num2str(yyuse(ii)) '.nc'];
    bonk = dir(mls_W);
    mls_W = [mlsdata 'ML3MBH2O_' vers '/' bonk.name];
  mls_O3 = [mlsdata 'ML3MBO3_' vers '/MLS-Aura_L3MB-O3_' xyvers '*' num2str(yyuse(ii)) '.nc'];
    bonk = dir(mls_O3);
    mls_O3 = [mlsdata 'ML3MBO3_' vers '/' bonk.name];
    
  fprintf(1,' T  = %s \n',mls_T)
  fprintf(1,' W  = %s \n',mls_W)
  fprintf(1,' O3 = %s \n',mls_O3)

  [xaT,xaW,xaO3] = mls_reader_L3(mls_T,mls_W,mls_O3,mmuse(ii));
  %% T,O3  works from level 08 to 55 (261 to 0,001 mb)
  %% WV works from level 07 to 55 (316 to 0,001 mb)

  mls_plev = xaT.Temperature_PressureGrid.lev;
  mls_lon  = xaT.Temperature_PressureGrid.lon;
  mls_lat  = xaT.Temperature_PressureGrid.lat;
  [mls_LAT,mls_LON] = ndgrid(mls_lat,mls_lon);

  %% see  ~/MATLABCODE/CONVERT_GAS_UNITS/toppmv.m
  %% elseif iGasUnitIN == 12    %%%in vmr (volume mixing ratio)
  %%  PPMV = VMR*1E+6
  %%  y = q*1e6;

  xaO3
  xaT
  xaW

%{
size(xaT.Temperature_PressureGrid.value)
    45    72    55    12
%}

  aT = squeeze(xaT.Temperature_PressureGrid.value(:,:,:,mmuse(ii)));
  aW = squeeze(xaW.H2O_PressureGrid.value(:,:,:,mmuse(ii))) *1e6;
  aO3 = squeeze(xaO3.O3_PressureGrid.value(:,:,:,mmuse(ii))) *1e6;;

  do_gridded_interpolant

  %% do_scattered_interpolant %%% VERY SLOW  XXXXXXXXXXXXX

  %% 250 mb
  figure(1); scatter_coast(pL3.rlon,pL3.rlat,50,pL3.ptemp(16,:)); title('AIRS L3 ptemp 250 mb'); cx = caxis;
  figure(2); scatter_coast(pL3.rlon,pL3.rlat,50,aTnew(08,:));     title('MLS  L3 ptemp 250 mb'); caxis(cx);

  figure(3); scatter_coast(pL3.rlon,pL3.rlat,50,pL3.gas_1(16,:)); title('AIRS L3 WV ppmv 250 mb'); cx = caxis;
  figure(4); scatter_coast(pL3.rlon,pL3.rlat,50,aWnew(08,:));     title('MLS  L3 WV ppmv 250 mb'); caxis(cx);

  figure(5); scatter_coast(pL3.rlon,pL3.rlat,50,pL3.gas_3(16,:)); title('AIRS L3 O3 ppmv 250 mb'); cx = caxis;
  figure(6); scatter_coast(pL3.rlon,pL3.rlat,50,aO3new(08,:));    title('MLS  L3 O3 ppmv 250 mb'); caxis(cx);

  %% 100 mb
  figure(1); scatter_coast(pL3.rlon,pL3.rlat,50,pL3.ptemp(13,:)); title('AIRS L3 ptemp 100 mb'); cx = caxis;
  figure(2); scatter_coast(pL3.rlon,pL3.rlat,50,aTnew(13,:));     title('MLS  L3 ptemp 100 mb'); caxis(cx);

  figure(3); scatter_coast(pL3.rlon,pL3.rlat,50,pL3.gas_1(13,:)); title('AIRS L3 WV ppmv 100 mb'); cx = caxis;
  figure(4); scatter_coast(pL3.rlon,pL3.rlat,50,aWnew(13,:));     title('MLS  L3 WV ppmv 100 mb'); caxis(cx);

  figure(5); scatter_coast(pL3.rlon,pL3.rlat,50,pL3.gas_3(13,:)); title('AIRS L3 O3 ppmv 100 mb'); cx = caxis;
  figure(6); scatter_coast(pL3.rlon,pL3.rlat,50,aO3new(13,:));    title('MLS  L3 O3 ppmv 100 mb'); caxis(cx);

  %% 050 mb
  figure(1); scatter_coast(pL3.rlon,pL3.rlat,50,pL3.ptemp(11,:)); title('AIRS L3 ptemp 050 mb'); cx = caxis;
  figure(2); scatter_coast(pL3.rlon,pL3.rlat,50,aTnew(17,:));     title('MLS  L3 ptemp 050 mb'); caxis(cx);

  figure(3); scatter_coast(pL3.rlon,pL3.rlat,50,pL3.gas_1(11,:)); title('AIRS L3 WV ppmv 050 mb'); cx = caxis;
  figure(4); scatter_coast(pL3.rlon,pL3.rlat,50,aWnew(17,:));     title('MLS  L3 WV ppmv 050 mb'); caxis(cx);

  figure(5); scatter_coast(pL3.rlon,pL3.rlat,50,pL3.gas_3(11,:)); title('AIRS L3 O3 ppmv 050 mb'); cx = caxis;
  figure(6); scatter_coast(pL3.rlon,pL3.rlat,50,aO3new(17,:));    title('MLS  L3 O3 ppmv 050 mb'); caxis(cx);

  %% replace levels 16 to 24 (250 to 1 mb) from AIRS L3 with levels 08 to 55 (261 to 0.0001 mb) from MLS
  booMLS  = find(mls_plev >= 0.005 & mls_plev <= 270); mls_plev([booMLS(1) booMLS(end)])
  booMLS  = find(mls_plev >= 0.05 & mls_plev <= 270); mls_plev([booMLS(1) booMLS(end)])
  airs_plev = pL3.plevs(:,1);
  booAIRS = find(airs_plev >= 0.005 & airs_plev <= 260); airs_plev([booAIRS(1) booAIRS(end)])

  new_airs_plev = [flipud(mls_plev(booMLS))' airs_plev(booAIRS(end)+1:end)']; 
  pL3new = pL3;
  pL3new.ptime = pL3new.rtime;
  junk1_to_24  = 24-pL3.nlevs;
  pL3new.nlevs = length(booAIRS(end)+1:24) + length(booMLS) - junk1_to_24;
  pL3new.plevs = new_airs_plev' * ones(1,length(pL3.stemp));
  pL3new.ptemp = [flipud(aTnew(booMLS,:)); pL3.ptemp(booAIRS(end)+1:24,:)];
  pL3new.gas_1 = [flipud(aWnew(booMLS,:)); pL3.gas_1(booAIRS(end)+1:24,:)];
  pL3new.gas_3 = [flipud(aO3new(booMLS,:)); pL3.gas_3(booAIRS(end)+1:24,:)];

  %%%%%%%%%%%%%%%%%%%%%%%%%
  if length(find(isnan(pL3new.ptemp))) > 0
    error('bad ptemp');
  elseif length(find(isnan(pL3new.gas_1))) > 0
    error('bad gas_1');
  elseif length(find(isnan(pL3new.gas_3))) > 0
    error('bad gas_3');
  end
  
  [zzz,rrr] = find(pL3new.gas_1 < 0 | pL3new.gas_3 < 0 | pL3new.ptemp < 150);
  rrr = unique(rrr);
  for jjj = 1 : length(rrr)
    plevs = pL3new.plevs(:,rrr(jjj));
    ptemp = pL3new.ptemp(:,rrr(jjj));
    gas_1 = pL3new.gas_1(:,rrr(jjj));
    gas_3 = pL3new.gas_3(:,rrr(jjj));
    zgood = find(ptemp > 150 & gas_1 >= 0 & gas_3 >= 0);
    zbad = find(ptemp <= 150 | gas_1 < 0 | gas_3 < 0);
    ptemp(zbad) = interp1(log(plevs(zgood)),ptemp(zgood),log(plevs(zbad)),[],'extrap');  
    gas_1(zbad) = exp(interp1(log(plevs(zgood)),log(gas_1(zgood)),log(plevs(zbad)),[],'extrap'));  
    gas_3(zbad) = exp(interp1(log(plevs(zgood)),log(gas_3(zgood)),log(plevs(zbad)),[],'extrap'));  
    pL3new.ptemp(:,rrr(jjj)) = ptemp;
    pL3new.gas_1(:,rrr(jjj)) = gas_1;
    pL3new.gas_3(:,rrr(jjj)) = gas_3;
  end

  rrr = find(isnan(pL3new.stemp));
  for jjj = 1 : length(rrr)
    plevs = pL3new.plevs(:,rrr(jjj));
    ptemp = pL3new.ptemp(:,rrr(jjj));
    zgood = find(ptemp > 150);
    stemp = interp1(log(plevs(zgood)),ptemp(zgood),log(pL3new.spres(rrr(jjj))),[],'extrap');  
    pL3new.stemp(:,rrr(jjj)) = stemp;
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%

  figure(4);  
  wonk = pL3new.ptemp; wonk(wonk < 0) = NaN; semilogy(nanmean(wonk,2),mean(pL3new.plevs,2)); set(gca,'ydir','reverse'); ylim([min(new_airs_plev) max(new_airs_plev)])
  figure(5)
  wonk = pL3new.gas_1; wonk(wonk < 0) = NaN; loglog(nanmean(wonk,2),mean(pL3new.plevs,2)); set(gca,'ydir','reverse'); ylim([min(new_airs_plev) max(new_airs_plev)])
  figure(6)
  wonk = pL3new.gas_3; wonk(wonk < 0) = NaN; loglog(nanmean(wonk,2),mean(pL3new.plevs,2)); set(gca,'ydir','reverse'); ylim([min(new_airs_plev) max(new_airs_plev)])

  psavejunk = p;
  pasavejunk = pa; 
  
  p = pL3new;
  pa = paL3;
  h = hL3;
  ha = haL3;

  p.rlon = wrapTo180(p.rlon);
  [p,pa] = rtp_add_emis(p,pa);
  [p,bady0] = fix_nan_emis(p);
  
  pnew_ip = p;
  hnew_ip = hL3;
  hnew_ip.pfields = 1;

  if isfield(pnew_ip,'palts')
    pnew_ip = rmfield(pnew_ip,'palts');
  end
  if isfield(pnew_ip,'mmw')
    pnew_ip = rmfield(pnew_ip,'mmw');
  end
  
  [xyy,xmm,xdd,xhh] = tai2utcSergio(pnew_ip.rtime);        % <<<<<<<<<<<<<<<<<<<<<<<< for SdSM old time
  time_so_far = (xyy-2000) + ((xmm-1)+1)/12;
  co2ppm = 368 + 2.077*time_so_far;  %% 395.6933
  pnew_ip.co2ppm = co2ppm;
  fprintf(1,'CLIMATOLOGY co2ppm for FIRST %4i/%2i/%2i = %8.6f ppmv\n',xyy(1),xmm(1),xdd(1),pnew_ip.co2ppm(1));
  fprintf(1,'CLIMATOLOGY co2ppm for LAST  %4i/%2i/%2i = %8.6f ppmv\n',xyy(end),xmm(end),xdd(end),pnew_ip.co2ppm(end));
  figure(7); scatter_coast(pnew_ip.rlon,pnew_ip.rlat,50,pnew_ip.co2ppm);
  figure(8); scatter_coast(pnew_ip.rlon,pnew_ip.rlat,50,pnew_ip.stemp); title('stemp')
  pause(0.1)

  run_sarta.klayers_code = '/asl/packages/klayersV205/BinV201/klayers_airs';
  %run_sarta.klayers_code = '/home/sergio/KLAYERS/BinV201/klayers_airs_wetwater_140levs';
  code0 = '/asl/packages/sartaV108/BinV201/sarta_apr08_m140_iceaggr_waterdrop_desertdust_slabcloud_hg3_wcon_nte';
  code1 = '/home/sergio/SARTA_CLOUDY/BinV201/sarta_apr08_m140x_iceGHMbaum_waterdrop_desertdust_slabcloud_hg3';
  code1 = '/home/sergio/SARTA_CLOUDY/BinV201/xsarta_apr08_m140_iceGHMbaum_waterdrop_desertdust_slabcloud_hg3';
  code1 = '/home/chepplew/gitLib/sarta/bin/airs_l1c_2834_cloudy_may19_prod_v3';

  run_sarta.klayers_code = klayers;
  code0 = sartaCld;
  code1 = sartaCld;
  
  run_sarta.clear = +1;
  run_sarta.cloud = +1;
  run_sarta.cumsum = 9999;  %% larrabee likes this, puts clouds high so does well for DCC
  run_sarta.cumsum = -1;    %% this is "closer" to MRO but since cliuds are at centroid, does not do too well with DCC
  run_sarta.sartaclear_code = code1;
  run_sarta.sartacloud_code = code1;
  run_sarta.co2ppm = co2ppm;

  rtpwrite(fip,hnew_ip,ha,pnew_ip,pa);  
  klayerser = ['!' run_sarta.klayers_code '      fin=' fip ' fout=' fop ' >& ugh']; eval(klayerser);
  [hnew_op,ha2,pnew_op,pa2] = rtpread(fop);
  pnew_op = make_rtp_plays(pnew_op);
  figure(9); clf; 
  jjj = 1141;  
    semilogy(pnew_ip.ptemp(:,jjj),pnew_ip.plevs(:,jjj),pnew_op.ptemp(:,jjj),pnew_op.plays(:,jjj)); set(gca,'ydir','reverse'); ylim([1e-2 pnew_op.spres(jjj)])
    loglog(pnew_ip.gas_1(:,jjj),pnew_ip.plevs(:,jjj),pnew_op.gas_1(:,jjj)/1e20,pnew_op.plays(:,jjj)); set(gca,'ydir','reverse'); ylim([1e-2 pnew_op.spres(jjj)])

  sartaer   = ['!' run_sarta.sartacloud_code '   fin=' fop ' fout=' frp ' >& ugh']; eval(sartaer);
  [hnew_op,ha2,pnew_op,pa2] = rtpread(frp);
  pnew_op = make_rtp_plays(pnew_op);

  %pnew_op.rcalc            = p2.rcalc;
  %pnew_op.sarta_rclearcalc = p2.sarta_rclearcalc;
  pnew_op.mmw = mmwater_rtp(hnew_op,pnew_op);

  %[hnew_ip,ha,pnew_ip,pa] = rtptrim_sartacloud(hnew_ip,ha,p2,pa);
  rmer = ['!/bin/rm '  fip ' ' fop ' ' frp]; eval(rmer);

  [Tw,Tw1km,Tdew,WBGT,RH,RH1km,colwater,TwSurf,RHSurf,TdewSurf] = layeramt2RH_wet_bulb_dew_point_temperature(hnew_op,pnew_op);
  %if isfield(pnew_ip,'d2m') & isfield(pnew_ip,'t2m')
  %  pnew_op.d2m = pnew_ip.d2m;
  %  pnew_op.t2m = pnew_ip.t2m;
  %  pnew_op.rh2m = airtemp_dewpointtemp_2_RH(pnew_op.t2m,pnew_op.d2m);
  %end
  pnew_op.Tw   = Tw;
  pnew_op.Tdew = Tdew;
  pnew_op.RH   = RH;
  pnew_op.mmw  = colwater;
  pnew_op.TwSurf   = TwSurf;
  pnew_op.TdewSurf = TdewSurf;
  pnew_op.RHSurf   = RHSurf;

  pnew_ip.rh = convert_humidity (pnew_ip.plevs*100,pnew_ip.ptemp,pnew_ip.gas_1,'specific humidity','relative humidity');

  yyuseII = yyuse(ii);
  mmuseII = mmuse(ii);
  comment = 'see /home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/driver_compute_mls_profile_rtpfiles.m';
  comment = 'see /home/sergio/MATLABCODE/RTPMAKE/CLUST_RTPMAKE/CLUSTMAKE_MLS/clust_compute_mls_profile_rtpfiles.m';
  saver = ['save ' fout ' comment        hnew_ip ha pnew_ip pa     hnew_op ha2 pnew_op pa2    yyuseII mmuseII'];
  if ~exist(fout)
    [yyuseII mmuseII]
    fprintf(1,'saving %s \n',fout)     
    eval(saver)
  else
    fprintf(1,'%s already exists ... not saving \n',fout) 
  end

  figure(1); scatter_coast(pnew_op.rlon,pnew_op.rlat,50,pnew_op.stemp); colormap jet; title('stemp (K)');
  figure(2); scatter_coast(pnew_op.rlon,pnew_op.rlat,50,pnew_op.mmw); colormap jet; title('col water (mm)');
  figure(3); scatter_coast(pnew_op.rlon,pnew_op.rlat,50,pnew_op.spres); colormap jet; title('spres (mb)');
  figure(4); scatter_coast(pnew_op.rlon,pnew_op.rlat,50,rad2bt(1231,pnew_op.rcalc(1520,:))); colormap jet; title('BT1231 (K)');

  i1419 = find(hnew_ip.vchan >= 1419,1); 
  figure(4); scatter_coast(pnew_op.rlon,pnew_op.rlat,50,rad2bt(1419,pnew_op.rcalc(i1419,:))); colormap jet; title('BT1419 (K)');

  [zzz,rrr] = find(isnan(pnew_op.rcalc));
  figure(3); plot(hnew_ip.vchan,std(rad2bt(hnew_ip.vchan,pnew_op.rcalc)'))
  figure(3); plot(hnew_ip.vchan,mean(rad2bt(hnew_ip.vchan,pnew_op.rcalc)'))
  pause(0.1)

end

rmer = ['!/bin/rm ' fip ' ' fop ' ' frp];
eval(rmer);
disp('now go do trends in /home/sergio/git/oem_climate_code/FIND_NWP_MODEL_TRENDS/driver_computeMLS_monthly_trends.m')

