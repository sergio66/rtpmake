JOB = str2num(getenv('SLURM_ARRAY_TASK_ID'));

%% 388 timesteps from 2002/09 to 2019/08 which can use ERA
%% after that need ECMWF

%JOB = 387  %% already done  7/30/2019 till 8/14/2019
%JOB = 388  %% alreadt done  8/15/2019 till 8/31/2019
%JOB = 38  %% 
%JOB = 261

addpath /asl/matlib/rtptools/
addpath /asl/matlib/aslutil
addpath /asl/matlib/h4tools
addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/TIME
addpath /home/sergio/MATLABCODE/PLOTTER
addpath ../GRIB
addpath /home/sergio/MATLABCODE/matlib/clouds/sarta
addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS

system_slurm_stats

switchERAtoECM = utc2taiSergio(2019,09,01,12);  %% this is when ERA ends

iDorA = -1;  %% asc
iDorA = +1;  %% desc

iFIll_ERA_Interp = -1; %% just use nearest in time
iFIll_ERA_Interp = +1; %% default

N = 16; N = 15; N1 = 1;
  firstORend = 0;
  iPrint = -1;

yy0 = 2002; mm0 = 09; dd0 = 01;
thedateS(1,:) = [yy0 mm0 dd0];
rtimeS(1) = utc2taiSergio(yy0,mm0,dd0,12);

[yy1,mm1,dd1] = addNdays(yy0,mm0,dd0,N,firstORend,iPrint);
thedateE(1,:) = [yy1 mm1 dd1];
rtimeE(1) = utc2taiSergio(yy1,mm1,dd1,12);

for ii = 2 : JOB
  yy0 = yy1; mm0 = mm1; dd0 = dd1;

  %% go forward by 1
  [yy1,mm1,dd1] = addNdays(yy0,mm0,dd0,N1,firstORend,iPrint);
  yy0 = yy1; mm0 = mm1; dd0 = dd1;
  thedateS(ii,:) = [yy0 mm0 dd0];  
  rtimeS(ii) = utc2taiSergio(yy0,mm0,dd0,12);

  [yy1,mm1,dd1] = addNdays(yy0,mm0,dd0,N,firstORend,iPrint);
  thedateE(ii,:) = [yy1 mm1 dd1];
  rtimeE(ii) = utc2taiSergio(yy1,mm1,dd1,12);

  if switchERAtoECM > rtimeE(ii) 
    fprintf(1,' ERA ii=%4i   Start %4i/%2i/%2i  End %4i/%2i/%2i rtimeS/E = %12.10e %12.10e\n',ii,thedateS(ii,:),thedateE(ii,:),rtimeS(ii),rtimeE(ii))
  elseif switchERAtoECM <= rtimeE(ii) 
    fprintf(1,' ECM ii=%4i   Start %4i/%2i/%2i  End %4i/%2i/%2i rtimeS/E = %12.10e %12.10e\n',ii,thedateS(ii,:),thedateE(ii,:),rtimeS(ii),rtimeE(ii))
  end
end  

[thedateS thedateE]
[yyM,mmM,ddM] = addNdays(thedateS(JOB,1),thedateS(JOB,2),thedateS(JOB,3),8,firstORend,iPrint);

fprintf(1,'JOB = %3i spans %4i/%2i/%2i to %4i/%2i/%2i both ends inclusive \n',JOB,[thedateS(JOB,:) thedateE(JOB,:)]);
fprintf(1,'          midpoint %4i/%2i/%2i \n',yyM,mmM,ddM);

[h,ha,p,pa] = rtpread('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/summary_17years_all_lat_all_lon_2002_2019_palts_startSept2002_CLEAR.rtp');

%% now we need to get overpass times and solzen angles, just set scanang to 22 deg
%% VERY VERY IMPORTANT : see driver_fix_thedata_asc_desc_solzen_time_412_64x72.m
oldwrong = load('/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_For_HowardObs_TimeSeries/asc_desc_solzen_time_412_64x72.mat');  %% this had THREE skips of data so bad offsets
if iDorA == 1
  load('/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_For_HowardObs_TimeSeries/asc_desc_solzen_time_412_64x72_fix_desc.mat');
elseif iDorA == -1
  load('/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_For_HowardObs_TimeSeries/asc_desc_solzen_time_412_64x72_fix_asc.mat');
end

monitor_memory_whos;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% from comment, see /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_For_HowardObs_TimeSeries/driver_loop_get_asc_desc_solzen_time.m
if iDorA > 0
  meanrtime = nanmean(squeeze(nanmean(thedata.rtime_desc,1)),1);
  minrtime  = nanmin(squeeze(nanmin(thedata.rtime_desc,[],1)),[],1);
  maxrtime  = nanmax(squeeze(nanmax(thedata.rtime_desc,[],1)),[],1);
  rlon   = thedata.rlon_desc(:,:,JOB);     rlon = rlon(:)';
  rlat   = thedata.rlat_desc(:,:,JOB);     rlat = rlat(:)';
  solzen = thedata.solzen_desc(:,:,JOB);   solzen = solzen(:)';
  satzen = thedata.satzen_desc(:,:,JOB);   satzen = satzen(:)';
  hour   = thedata.hour_desc(:,:,JOB);     hour = hour(:)';
  bt1231 = thedata.bt1231_desc(:,:,JOB,2); bt1231 = bt1231(:)';
  rtime  = thedata.rtime_desc(:,:,JOB);    rtime = rtime(:)';
elseif iDorA < 0
  meanrtime = nanmean(squeeze(nanmean(thedata.rtime_asc,1)),1);
  minrtime  = nanmin(squeeze(nanmin(thedata.rtime_asc,[],1)),[],1);
  maxrtime  = nanmax(squeeze(nanmax(thedata.rtime_asc,[],1)),[],1);
  rlon   = thedata.rlon_asc(:,:,JOB);     rlon = rlon(:)';
  rlat   = thedata.rlat_asc(:,:,JOB);     rlat = rlat(:)';
  solzen = thedata.solzen_asc(:,:,JOB);   solzen = solzen(:)';
  satzen = thedata.satzen_asc(:,:,JOB);   satzen = satzen(:)';
  hour   = thedata.hour_asc(:,:,JOB);     hour = hour(:)';
  bt1231 = thedata.bt1231_asc(:,:,JOB,2); bt1231 = bt1231(:)';
  rtime  = thedata.rtime_asc(:,:,JOB);    rtime = rtime(:)';
end

%% now check rtimeS < meanrtime < rtimeE
if iDorA > 0
  for ii = 1 : JOB
    tt=ii; junk=thedata.rtime_desc(:,:,tt); junk=min(junk(:)); [junkYY,junkMM,junkDD,junkHH] = tai2utcSergio(junk); 
    fprintf(1,' %3i  tS tData tE  %4i/%2i/%2i     %4i/%2i/%2i     %4i/%2i/%2i     deltaRtime(Max-Junk) = %8.6f\n',ii,thedateS(tt,:),junkYY,junkMM,junkDD,thedateE(tt,:),(rtimeE(ii)-junk)/86400)
  end
else
  for ii = 1 : JOB
    tt=ii; junk=thedata.rtime_asc(:,:,tt); junk=min(junk(:)); [junkYY,junkMM,junkDD,junkHH] = tai2utcSergio(junk); 
    fprintf(1,' %3i  tS tData tE  %4i/%2i/%2i     %4i/%2i/%2i     %4i/%2i/%2i     deltaRtime(Max-Junk) = %8.6f\n',ii,thedateS(tt,:),junkYY,junkMM,junkDD,thedateE(tt,:),(rtimeE(ii)-junk)/86400)
  end
end

boo = rtimeE-rtimeS; boo = find(boo < 0);
if length(boo) > 0
  boo
  error('oops found some rtimeE-rtimeS < 0');
end
boo = meanrtime(1:JOB)-rtimeS; boo = find(boo < 0);
if length(boo) > 0
  boo
  error('oops found some meanrtime(1:JOB)-rtimeS < 0');
end
boo = rtimeE-meanrtime(1:JOB); boo = find(boo < 0);
if length(boo) > 0
  boo
  error('oops found some rtimeE-meanrtime(1:JOB) < 0');
end

[yyMean,mmMean,ddMean,hhMean] = tai2utcSergio(meanrtime);
[yyMin,mmMin,ddMin,hhMin] = tai2utcSergio(minrtime);
[yyMax,mmMax,ddMax,hhMax] = tai2utcSergio(maxrtime);
figure(1); plot(yyMean);
figure(2); plot(mmMean);
figure(3); plot(ddMean);
figure(4); plot(hhMean);
fprintf(1,'yy/mm/dd hh for JOB %3i min  = %4i/%2i/%2i %8.6f \n',JOB,yyMin(JOB),mmMin(JOB),ddMin(JOB),hhMin(JOB))
fprintf(1,'                        mean = %4i/%2i/%2i %8.6f \n',yyMean(JOB),mmMean(JOB),ddMean(JOB),hhMean(JOB))
fprintf(1,'                        max  = %4i/%2i/%2i %8.6f \n',yyMax(JOB),mmMax(JOB),ddMax(JOB),hhMax(JOB))

figure(1); scatter_coast(rlon,rlat,50,solzen); colormap jet
figure(2); scatter_coast(rlon,rlat,50,satzen); colormap jet
figure(3); scatter_coast(rlon,rlat,50,hour); colormap jet
figure(4); scatter_coast(rlon,rlat,50,bt1231); colormap jet
figure(5); scatter_coast(rlon,rlat,50,rlon-p.rlon); colormap jet
figure(6); scatter_coast(rlon,rlat,50,rlat-p.rlat); colormap jet

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% dirs used by /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_For_HowardObs_TimeSeries
if abs(iFIll_ERA_Interp) ~= 1
  iFIll_ERA_Interp
  error('need iFIll_ERA_Interp == +/- 1')
end

if iFIll_ERA_Interp == +1
  if iDorA > 0
    fout = ['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA/Tile_Center/DESC/era_tile_center_timestep_' num2str(JOB,'%03d') '.mat'];
  elseif iDorA < 0
    fout = ['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA/Tile_Center/ASC/era_tile_center_timestep_' num2str(JOB,'%03d') '.mat'];
  end
else
  if iDorA > 0
    fout = ['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA/Tile_Center/DESC/era_tile_center_timestep_' num2str(JOB,'%03d') '_closestINtime.mat'];
  elseif iDorA < 0
    fout = ['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA/Tile_Center/ASC/era_tile_center_timestep_' num2str(JOB,'%03d') '_closestINtime.mat'];
  end
end

iDo = +1;
if exist(fout)
  fprintf(1,'JOB %3i : avg era timeseries file %s already exists \n',JOB,fout);
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
  hnew_ip = rmfield(hnew_ip,'ptype');
  hnew_ip = rmfield(hnew_ip,'pfields');

  pnew_ip.solzen = solzen;
  pnew_ip.satzen = satzen;
  pnew_ip.scanang = saconv(p.satzen,p.zobs);  
  pnew_ip.rtime   = rtime;

  pnew_ip = rmfield(pnew_ip,'stemp');
  pnew_ip = rmfield(pnew_ip,'ptemp');
  pnew_ip = rmfield(pnew_ip,'plevs');
  pnew_ip = rmfield(pnew_ip,'palts'); %% VERY IMPORTANT, new Aug 2021
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

  %[pnew_ip,hnew_ip] = fill_era(pnew_ip,hnew_ip);
  %[pnew_ip,hnew_ip] = fill_era_interp(pnew_ip,hnew_ip);
  %[pnew_ip,hnew_ip] = fill_era(pnew_ip,hnew_ip);
  hnew_ip.pfields = 1;   %% so say there is no profile and no rcalc
  if rtimeE(JOB) < switchERAtoECM
    %% 2002/09 to 2019/08
    if iFIll_ERA_Interp == +1
      [pnew_ip,hnew_ip] = fill_era_interp(pnew_ip,hnew_ip);
    elseif iFIll_ERA_Interp == -1
      [pnew_ip,hnew_ip] = fill_era(pnew_ip,hnew_ip);   
    end
  else
    %% 2019/09 onwards
    %[pnew_ip,hnew_ip] = fill_ecmwf_interp(pnew_ip,hnew_ip);
    [pnew_ip,hnew_ip] = fill_ecmwf(pnew_ip,hnew_ip);
  end

  [xyy,xmm,xdd,xhh] = tai2utcSergio(pnew_ip.rtime);        % <<<<<<<<<<<<<<<<<<<<<<<< for SdSM old time
  time_so_far = (xyy-2000) + ((xmm-1)+1)/12;
  co2ppm = 368 + 2.077*time_so_far;  %% 395.6933
  pnew_ip.co2ppm = co2ppm;
  fprintf(1,'CLIMATOLOGY co2ppm for FIRST %4i/%2i/%2i = %8.6f ppmv\n',xyy(1),xmm(1),xdd(1),pnew_ip.co2ppm(1));
  fprintf(1,'CLIMATOLOGY co2ppm for LAST  %4i/%2i/%2i = %8.6f ppmv\n',xyy(end),xmm(end),xdd(end),pnew_ip.co2ppm(end));
  scatter_coast(pnew_ip.rlon,pnew_ip.rlat,50,pnew_ip.co2ppm);

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
  pnew_op.rcalc = p2.rcalc;
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

  comment = 'see /home/sergio/MATLABCODE/RTPMAKE/CLUST_RTPMAKE/CLUSTMAKE_ERA/clust_loop_make16day_tile_center.m';
  saver = ['save ' fout ' comment        hnew_ip ha pnew_ip pa     hnew_op ha2 pnew_op pa2 yyM mmM ddM '];
  eval(saver)
  fprintf(1,'saving to %s \n',fout)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
%% for files made before Aug 5, 2021 have to fix filenames since there was a dumb 16day timestep offset (missing data) for steps 169, 394, 409
dOrig = '/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA/Tile_Center/DESC/'
dTemp = '/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA/Tile_Center/DESC/TEMP/'
fname = 'era_tile_center_timestep_';
if ~exist(dTemp)
  mker = ['!mkdir ' dTemp];
  fprintf(1,'making sTemp %s \n',mker)
  eval(mker);
end

%% move 169-393, upping index by one, so will have 170-394   NO 169, and will put 394 into the above segment -- SO NO 169 394 .. have to run those
for ii = 169:393
  mver = ['!/bin/mv ' dOrig '/' fname num2str(ii,'%03d') '.mat ' dTemp '/' fname num2str(ii+1,'%03d') '.mat'];
  fprintf(1,'%s \n',mver)
  eval(mver)
end

%move all renamed files back from temp dir to current dir 
mver = ['!/bin/mv ' dTemp '/' fname '*.mat  ' dOrig '/.']
eval(mver)

rmer = ['!/bin/rmdir ' dTemp];
eval(rmer)
%}
