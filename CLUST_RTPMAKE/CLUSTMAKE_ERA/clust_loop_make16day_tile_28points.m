%% takes 3 minutes per rtp set of 4608 profiles, with 28 grid points per day, that is 90 minutes ... with 16 days that is 24 hours on strowinteract

JOB = str2num(getenv('SLURM_ARRAY_TASK_ID'));

%% 388 timesteps from 2002/09 to 2019/08 which can use ERA
%% after that need ECMWF

%JOB = 387  %% already done  7/30/2019 till 8/14/2019
%JOB = 388  %% alreadt done  8/15/2019 till 8/31/2019
%JOB = 389  %% 

%% 23/year so 365/16 * 17 ~ 387 since we go 2002/09 to 2019/08  use ERA from 2002/09-2019/08
%% 23/year so 365/16 * 18 ~ 410 since we go 2002/09 to 2020/08  use ECM from 2019/09-2020/08
%% 23/year so 365/16 * 19 ~ 433 since we go 2002/09 to 2021/08  use ECM from 2019/09-2021/08

addpath /asl/matlib/rtptools/
addpath /asl/matlib/aslutil
addpath /asl/matlib/h4tools
addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/TIME
addpath /home/sergio/MATLABCODE/PLOTTER
addpath ../GRIB
addpath /home/sergio/MATLABCODE/matlib/clouds/sarta

system_slurm_stats

switchERAtoECM = utc2taiSergio(2019,09,01,12);  %% this is when ERA ends

iDorA = +1;  %% desc

iFIll_ERA_Interp = +1; %% default
iFIll_ERA_Interp = -1; %% just use nearest in time

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

%junk = [(1:JOB)' thedateS thedateE rtimeS' rtimeE']; fprintf(1,'%3i %4i/%2i/%2i %4i/%2i/%2i %8.6e %8.6e \n',junk')

fprintf(1,'JOB = %3i spans %4i/%2i/%2i to %4i/%2i/%2i both ends inclusive \n',JOB,[thedateS(JOB,:) thedateE(JOB,:)]);

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

%% this gives the 7x4 points in every tile
load('/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_For_HowardObs_TimeSeries/theERAgridptslist_ERA_tile_points.mat')
for ii = 1 : 72
  for jj = 1 : 64
    eraX(ii,jj) = theERAgridptslist.adjusted_xvalues{ii,jj}(1);
    eraY(ii,jj) = theERAgridptslist.adjusted_yvalues{ii,jj}(1);

    eraX(ii,jj) = theERAgridptslist.adjusted_xvalues{ii,jj}(4); %% 1 2 3 <4> 5 6 7 == center
    eraY(ii,jj) = theERAgridptslist.adjusted_yvalues{ii,jj}(2); %% 1 2 3 4 ~ center
  end
end

if iDorA < 0
  tilecenterX = squeeze(thedata.rlon_asc(:,:,1));
  tilecenterY = squeeze(thedata.rlat_asc(:,:,1));
else
  tilecenterX = squeeze(thedata.rlon_desc(:,:,1));
  tilecenterY = squeeze(thedata.rlat_desc(:,:,1));
end
plot(tilecenterX,tilecenterY,'b.',eraX,eraY,'ro')
plot(tilecenterX-eraX,tilecenterY-eraY,'b.')
plot(eraX(:,21)- tilecenterX(:,21))
plot(eraY(21,:)- tilecenterY(21,:))

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
else
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
for ii = 1 : JOB
  tt=ii; junk=thedata.rtime_desc(:,:,tt); junk=min(junk(:)); [junkYY,junkMM,junkDD,junkHH] = tai2utcSergio(junk); 
  fprintf(1,' %3i  tS tData tE  %4i/%2i/%2i     %4i/%2i/%2i     %4i/%2i/%2i     deltaRtime(Max-Junk) = %8.6f\n',ii,thedateS(tt,:),junkYY,junkMM,junkDD,thedateE(tt,:),(rtimeE(ii)-junk)/86400)
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

for ddd = 16 : -1 : 1      %% 16 day averages
  eeeXY_cnt = 0;

  for  eeeX = 1 : 7   %% 7 X ERA grid points
    for eeeY = 1 : 4  %% 4Y ERA grid points

      eeeXY_cnt =   eeeXY_cnt + 1;
      %% dirs used by /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_For_HowardObs_TimeSeries
      if iDorA > 0
        dout = ['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA/Tile_16day/DESC/Day' num2str(ddd,'%02d') '/ERAindex' num2str(eeeXY_cnt,'%02d')  '/'];
      else
        dout = ['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA/Tile_16day/ASC/Day' num2str(ddd,'%02d') '/ERAindex' num2str(eeeXY_cnt,'%02d')  '/'];
      end
      if ~exist(dout)
        mker = ['!mkdir -p ' dout];
        eval(mker)
      end

      for ii = 1 : 72
        for jj = 1 : 64
          eraX(ii,jj) = theERAgridptslist.adjusted_xvalues{ii,jj}(1);
          eraY(ii,jj) = theERAgridptslist.adjusted_yvalues{ii,jj}(1);
  
          eraX(ii,jj) = theERAgridptslist.adjusted_xvalues{ii,jj}(4); %% 1 2 3 <4> 5 6 7 == center
          eraY(ii,jj) = theERAgridptslist.adjusted_yvalues{ii,jj}(2); %% 1 2 3 4 ~ center

          eraX(ii,jj) = theERAgridptslist.adjusted_xvalues{ii,jj}(eeeX); 
          eraY(ii,jj) = theERAgridptslist.adjusted_yvalues{ii,jj}(eeeY); 
        end
      end
  
      if iFIll_ERA_Interp == +1
        fout = [dout '/era_tile_X_' num2str(eeeX) '_Y_' num2str(eeeY)  '_day_' num2str(ddd,'%02d') '_individual_timestep_' num2str(JOB,'%03d') '.mat'];
      elseif iFIll_ERA_Interp == -1
        fout = [dout '/era_tile_X_' num2str(eeeX) '_Y_' num2str(eeeY)  '_day_' num2str(ddd,'%02d') '_individual_timestep_' num2str(JOB,'%03d') '_closestINtime.mat'];
      end

      iDo = +1;    
      if exist(fout)
        fprintf(1,'JOB %3i : individual era timeseries file %s already exists \n',JOB,fout);
        iDo = -1;
      end

      if iDo > 0
        % now fill in the fields see cloud_set_defaults_run_maker.m
        pnew_ip = p;
        hnew_ip = h;
        hnew_ip = rmfield(hnew_ip,'ngas');
        hnew_ip = rmfield(hnew_ip,'glist');
        hnew_ip = rmfield(hnew_ip,'gunit');
        hnew_ip = rmfield(hnew_ip,'ptype');
        hnew_ip = rmfield(hnew_ip,'pfields');
      
        pnew_ip.rlon   = eraX(:)';
        pnew_ip.rlat   = eraY(:)';

        pnew_ip.solzen = solzen;
        pnew_ip.satzen = satzen;
        pnew_ip.scanang = saconv(p.satzen,p.zobs);  
        pnew_ip.rtime   = rtime;

        figure(3); plot(p.rlon,p.rlat,'.')
        figure(4); plot(pnew_ip.rlon,pnew_ip.rlat,'.')
        figure(5); plot(p.rlon-pnew_ip.rlon,p.rlat-pnew_ip.rlat,'.')

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
        hnew_ip.pfields = 0;   %% so say there is no profile and no rcalc
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
      
        comment = 'see /home/sergio/MATLABCODE/RTPMAKE/CLUST_RTPMAKE/CLUSTMAKE_ERA/clust_loop_make16day_tile_28points.m';
        saver = ['save ' fout ' comment        hnew_ip ha pnew_ip pa     hnew_op ha2 pnew_op pa2'];
        saver = ['save ' fout ' comment                                  hnew_op ha2 pnew_op pa2'];
        eval(saver)
      end
    end     %% loop over Y ERA grid point index
  end       %% loop over X ERA grid point index
end         %% loop over 16 days

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
%% for files made before Aug 5, 2021 have to fix filenames since there was a dumb 16day timestep offset (missing data) for steps 169, 394, 409
[sergio@strow-interact CLUSTMAKE_ERA]$ ls /asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA/Tile_16day/DESC/Day01/ERAindex*/era_tile*_day_01_individual_timestep_291.mat
/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA/Tile_16day/DESC/Day01/ERAindex01/era_tile_X_1_Y_1_day_01_individual_timestep_291.mat
/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA/Tile_16day/DESC/Day01/ERAindex02/era_tile_X_1_Y_2_day_01_individual_timestep_291.mat
/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA/Tile_16day/DESC/Day01/ERAindex03/era_tile_X_1_Y_3_day_01_individual_timestep_291.mat
/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA/Tile_16day/DESC/Day01/ERAindex04/era_tile_X_1_Y_4_day_01_individual_timestep_291.mat
/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA/Tile_16day/DESC/Day01/ERAindex05/era_tile_X_2_Y_1_day_01_individual_timestep_291.mat
/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA/Tile_16day/DESC/Day01/ERAindex06/era_tile_X_2_Y_2_day_01_individual_timestep_291.mat
/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA/Tile_16day/DESC/Day01/ERAindex07/era_tile_X_2_Y_3_day_01_individual_timestep_291.mat
/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA/Tile_16day/DESC/Day01/ERAindex08/era_tile_X_2_Y_4_day_01_individual_timestep_291.mat
/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA/Tile_16day/DESC/Day01/ERAindex09/era_tile_X_3_Y_1_day_01_individual_timestep_291.mat
/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA/Tile_16day/DESC/Day01/ERAindex10/era_tile_X_3_Y_2_day_01_individual_timestep_291.mat
/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA/Tile_16day/DESC/Day01/ERAindex11/era_tile_X_3_Y_3_day_01_individual_timestep_291.mat
/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA/Tile_16day/DESC/Day01/ERAindex12/era_tile_X_3_Y_4_day_01_individual_timestep_291.mat
/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA/Tile_16day/DESC/Day01/ERAindex13/era_tile_X_4_Y_1_day_01_individual_timestep_291.mat
/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA/Tile_16day/DESC/Day01/ERAindex14/era_tile_X_4_Y_2_day_01_individual_timestep_291.mat
/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA/Tile_16day/DESC/Day01/ERAindex15/era_tile_X_4_Y_3_day_01_individual_timestep_291.mat
/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA/Tile_16day/DESC/Day01/ERAindex16/era_tile_X_4_Y_4_day_01_individual_timestep_291.mat
/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA/Tile_16day/DESC/Day01/ERAindex17/era_tile_X_5_Y_1_day_01_individual_timestep_291.mat
/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA/Tile_16day/DESC/Day01/ERAindex18/era_tile_X_5_Y_2_day_01_individual_timestep_291.mat
/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA/Tile_16day/DESC/Day01/ERAindex19/era_tile_X_5_Y_3_day_01_individual_timestep_291.mat
/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA/Tile_16day/DESC/Day01/ERAindex20/era_tile_X_5_Y_4_day_01_individual_timestep_291.mat
/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA/Tile_16day/DESC/Day01/ERAindex21/era_tile_X_6_Y_1_day_01_individual_timestep_291.mat
/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA/Tile_16day/DESC/Day01/ERAindex22/era_tile_X_6_Y_2_day_01_individual_timestep_291.mat
/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA/Tile_16day/DESC/Day01/ERAindex23/era_tile_X_6_Y_3_day_01_individual_timestep_291.mat
/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA/Tile_16day/DESC/Day01/ERAindex24/era_tile_X_6_Y_4_day_01_individual_timestep_291.mat
/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA/Tile_16day/DESC/Day01/ERAindex25/era_tile_X_7_Y_1_day_01_individual_timestep_291.mat
/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA/Tile_16day/DESC/Day01/ERAindex26/era_tile_X_7_Y_2_day_01_individual_timestep_291.mat
/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA/Tile_16day/DESC/Day01/ERAindex27/era_tile_X_7_Y_3_day_01_individual_timestep_291.mat
/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA/Tile_16day/DESC/Day01/ERAindex28/era_tile_X_7_Y_4_day_01_individual_timestep_291.mat

for dd = 1 : 16
  for ee = 1 : 28
    dOrig = ['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA/Tile_16day/DESC/Day' num2str(dd,'%02d') '/ERAindex' num2str(ee,'%02d') '/'];
    dTemp = '/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA/Tile_16day/DESC/TEMP/'
    thedir = dir([dOrig '/era_tile_*_day_' num2str(dd,'%02d') '_individual_timestep_*.mat']);
    junkname = thedir(1).name;
    junkjunk = strfind(junkname,'_');    
    fname = junkname(1:junkjunk(end));
    fprintf(1,'dd = %2i ee = %2i fname = %s \n',dd,ee,fname);

    %% recall 4 lat NWP points, 7 lon NWP points   ee=1 X=1,Y=1     ee=2 X=1,Y=2    ee=3 X=1,Y=3    ee=4 X=1,Y=4        ee=5 X=2,Y=1
    %% so eg if ee=23,dd=1   fname=era_tile_X_6_Y_3_day_01_individual_timestep_     ee=23 means X_6_Y_3_  so ee = (X-1)*4+Y = (6-1)*4+3 = 20+3 = 23
    
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
  end
end
%}
