%% takes 3 minutes per rtp set of 4608 profiles, with 13*21 grid points per day, that is 900 minutes ... with 16 days that is 13 hours on strowinteract

%% if run.sarta == -1      then do    sbatch -p high_mem --array=1-240  sergio_matlab_jobB.sbatch 5    so can do complete time series for trending SKT
%% if run.sarta == 9999    then do    sbatch -p high_mem --array=120    sergio_matlab_jobB.sbatch 5    so can compare to histogram in Foig 8 of trends.tex

JOB = str2num(getenv('SLURM_ARRAY_TASK_ID'));  %% about 20years x 12 months = 240 months
if length(JOB) == 0
  JOB = 240;
  JOB = 120;
  JOB = 001;
end

addpath /asl/matlib/rtptools/
addpath /asl/matlib/aslutil
addpath /asl/matlib/h4tools
addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/TIME
addpath /home/sergio/MATLABCODE/PLOTTER
addpath ../GRIB
addpath /home/sergio/MATLABCODE/matlib/clouds/sarta

system_slurm_stats

iDorA = +1;  %% desc

iDo2m = -1;

iAllChan_or_1231 = -1;  %% only do one channel since this is my "clear BT1231 Q0.90 "

[h,ha,p,pa] = rtpread('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/summary_17years_all_lat_all_lon_2002_2019_palts_startSept2002_CLEAR.rtp');
if iAllChan_or_1231 < 0
  [h,p] = subset_rtp_allcloudfields(h,p,[],1291,[]);
end

get_dates_loop_make_monthly2m_tile_center_asc_or_desc

%% this gives the 13x21 points in every tile
%% /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_For_HowardObs_TimeSeries/test_read_eraI.m
%%   save /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_For_HowardObs_TimeSeries/theECM5gridptslist_ERA5_tile_points.mat theECM5gridptslist comment
load('/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_For_HowardObs_TimeSeries/theECMgridptslist_ECM_tile_points.mat');
for ii = 1 : 72
  for jj = 1 : 64
    eraX(ii,jj) = theECMgridptslist.adjusted_xvalues{ii,jj}(1);
    eraY(ii,jj) = theECMgridptslist.adjusted_yvalues{ii,jj}(1);

    eraX(ii,jj) = theECMgridptslist.adjusted_xvalues{ii,jj}(10); %% 1 2 3 .. 10 .... 21 == center
    eraY(ii,jj) = theECMgridptslist.adjusted_yvalues{ii,jj}(7); %% 1 2 3 .. 7  .....13 ~ center
  end
end
tilecenterX = squeeze(thedata.rlon_asc(:,:,1));
tilecenterY = squeeze(thedata.rlat_asc(:,:,1));
plot(tilecenterX,tilecenterY,'b.',eraX,eraY,'ro')
plot(tilecenterX-eraX,tilecenterY-eraY,'b.')
plot(eraX(:,21)- tilecenterX(:,21))
plot(eraY(21,:)- tilecenterY(21,:))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ddd = 1 : 1      %% monthly day averages, as we only get one data point per month
  eeeXY_cnt = 0;

  for  eeeX = 1 : 21   %% 7X ERA grid points
    for eeeY = 1 : 12  %% 4Y ERA grid points

      eeeXY_cnt =   eeeXY_cnt + 1;
      %% dirs used by /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_For_HowardObs_TimeSeries
      %% dout = ['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_16day/Day' num2str(ddd,'%02d') '/ERAindex' num2str(eeeXY_cnt,'%02d')  '/'];
      dout = ['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Monthly/ERAindex' num2str(eeeXY_cnt,'%02d')  '/'];
      if ~exist(dout)
        mker = ['!mkdir -p ' dout];
        eval(mker)
      end

      for ii = 1 : 72
        for jj = 1 : 64
          eraX(ii,jj) = theECMgridptslist.adjusted_xvalues{ii,jj}(1);
          eraY(ii,jj) = theECMgridptslist.adjusted_yvalues{ii,jj}(1);
  
          eraX(ii,jj) = theECMgridptslist.adjusted_xvalues{ii,jj}(4); %% 1 2 3 <4> 5 6 7 == center
          eraY(ii,jj) = theECMgridptslist.adjusted_yvalues{ii,jj}(2); %% 1 2 3 4 ~ center

          eraX(ii,jj) = theECMgridptslist.adjusted_xvalues{ii,jj}(eeeX); 
          eraY(ii,jj) = theECMgridptslist.adjusted_yvalues{ii,jj}(eeeY); 
        end
      end
  
      run_sarta.klayers_code = '/asl/packages/klayersV205/BinV201/klayers_airs';
      run_sarta.clear = +1;
      run_sarta.cloud = +1;
      run_sarta.cumsum = -1;    %% this is "closer" to MRO but since cliuds are at centroid, does not do too well with DCC  ... did this June 2024
      run_sarta.cumsum = 9999;  %% larrabee likes this, puts clouds high so does well for DCC  try this for timestep 12-0 = Fig 8 in trends.tex

      code0 = '/asl/packages/sartaV108/BinV201/sarta_apr08_m140_iceaggr_waterdrop_desertdust_slabcloud_hg3_wcon_nte';
      code1 = '/home/sergio/SARTA_CLOUDY/BinV201/sarta_apr08_m140x_iceGHMbaum_waterdrop_desertdust_slabcloud_hg3';
      code1 = '/home/sergio/SARTA_CLOUDY/BinV201/xsarta_apr08_m140_iceGHMbaum_waterdrop_desertdust_slabcloud_hg3';
      code1 = '/home/chepplew/gitLib/sarta/bin/airs_l1c_2834_cloudy_may19_prod_v3';
      run_sarta.sartaclear_code = code1;
      run_sarta.sartacloud_code = code1;

      %% fout = [dout '/era_tile_X_' num2str(eeeX) '_Y_' num2str(eeeY)  '_day_' num2str(ddd,'%02d') '_individual_timestep_' num2str(JOB,'%03d') '.mat'];
      if run_sarta.cumsum == -1
        fout = [dout '/era_tile_X_' num2str(eeeX) '_Y_' num2str(eeeY)  '_individual_timestep_' num2str(JOB,'%03d') '.mat'];
      elseif run_sarta.cumsum == 9999
        fout = [dout '/era_tile_X_' num2str(eeeX) '_Y_' num2str(eeeY)  '_individual_timestep_' num2str(JOB,'%03d') '_9999.mat'];
      end

      if exist(fout)
        fprintf(1,'JOB %3i : individual era timeseries file %s already exists \n',JOB,fout);
      else
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

        pnew_ip.rlon   = eraX(:)';
        pnew_ip.rlat   = eraY(:)';
        figure(3); plot(p.rlon,p.rlat,'.')
        figure(4); plot(pnew_ip.rlon,pnew_ip.rlat,'.')
        figure(5); plot(p.rlon-pnew_ip.rlon,p.rlat-pnew_ip.rlat,'.')

        pnew_ip = rmfield(pnew_ip,'stemp');
        pnew_ip = rmfield(pnew_ip,'ptemp');
        pnew_ip = rmfield(pnew_ip,'plevs');
        pnew_ip = rmfield(pnew_ip,'palts');
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
        %[pnew_ip,hnew_ip] = fill_era5_daily(pnew_ip,hnew_ip);      
        [pnew_ip,hnew_ip] = fill_era5_monthly(pnew_ip,hnew_ip);

        [xyy,xmm,xdd,xhh] = tai2utcSergio(pnew_ip.rtime);        % <<<<<<<<<<<<<<<<<<<<<<<< for SdSM old time
        time_so_far = (xyy-2000) + ((xmm-1)+1)/12;
        co2ppm = 368 + 2.077*time_so_far;  %% 395.6933
        pnew_ip.co2ppm = co2ppm;
        run_sarta.co2ppm = co2ppm;
        fprintf(1,'CLIMATOLOGY co2ppm for FIRST %4i/%2i/%2i = %8.6f ppmv\n',xyy(1),xmm(1),xdd(1),pnew_ip.co2ppm(1));
        fprintf(1,'CLIMATOLOGY co2ppm for LAST  %4i/%2i/%2i = %8.6f ppmv\n',xyy(end),xmm(end),xdd(end),pnew_ip.co2ppm(end));
        scatter_coast(pnew_ip.rlon,pnew_ip.rlat,50,pnew_ip.co2ppm);
      
printarray([min(pnew_ip.rlon) max(pnew_ip.rlon) min(pnew_ip.rlat) max(pnew_ip.rlat)],'in clust_loop_make_monthly_tile_273points.m  : min/max rlon  min.max rlat')
printarray([min(pnew_ip.plon) max(pnew_ip.plon) min(pnew_ip.plat) max(pnew_ip.plat)],'in clust_loop_make_monthly_tile_273points.m  : min/max plon  min.max plat')

        fip = mktemp('fx.ip.rtp');
        fop = mktemp('fx.op.rtp');
        frp = mktemp('fx.rp.rtp');
      
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
      
        comment = 'see /home/sergio/MATLABCODE/RTPMAKE/CLUST_RTPMAKE/CLUSTMAKE_ERA/clust_loop_make16day_tile_28points.m';
        saver = ['save ' fout ' comment                                  hnew_op ha2 pnew_op pa2'];

        comment = 'see /home/sergio/MATLABCODE/RTPMAKE/CLUST_RTPMAKE/CLUSTMAKE_ERA5/clust_loop_make_monthly_tile_273points.m';
        saver = ['save -v7.3 ' fout ' comment        hnew_ip ha pnew_ip pa     hnew_op ha2 pnew_op pa2'];

        eval(saver)
      end
    end     %% loop over Y ERA grid point index
  end       %% loop over X ERA grid point index
end         %% loop over the one day/month

disp('now call << driver_find_hottest_10percent_from_ERA5clearcalcs >> then << clust_find_hottest_10percent_from_ERA5clearcalcs >> ')
disp('now call << driver_find_hottest_10percent_from_ERA5clearcalcs >> then << clust_find_hottest_10percent_from_ERA5clearcalcs >> ')
disp('now call << driver_find_hottest_10percent_from_ERA5clearcalcs >> then << clust_find_hottest_10percent_from_ERA5clearcalcs >> ')
