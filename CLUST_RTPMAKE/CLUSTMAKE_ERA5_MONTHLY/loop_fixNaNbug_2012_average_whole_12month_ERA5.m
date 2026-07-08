%% one per month, 19 years of AIRS data so 19x12 = 228 sets of data
%% one per month, 19 years of AIRS data so 20x12 = 240 sets of data

%starttime = tic;
ticcStart0 = clock;

simulateYear = 2012;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% dirs used by /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_For_HowardObs_TimeSeries

iDorA = +1;

if iDorA > 0
  fin = ['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Center12months/DESC/' num2str(simulateYear,'%04d') '/era5_full12months_tile_center_' num2str(JOB,'%04d') '.mat']; %%% NOTE THIS IS DESC
  fout = ['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Center12months/DESC/' num2str(simulateYear,'%04d') '/FixedNAN/fixedNANavgbug_era5_full12months_tile_center_' num2str(JOB,'%04d') '.mat']; %%% NOTE THIS IS DESC
elseif iDorA < 0
  fin = ['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Center12months/ASC/' num2str(simulateYear,'%04d') '/era5_full12months_tile_center_' num2str(JOB,'%04d') '.mat']; %%% NOTE THIS IS DESC
  fout = ['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Center12months/ASC/' num2str(simulateYear,'%04d') '/FixedNAN/fixedNANavgbug_era5_full12months_tile_center_' num2str(JOB,'%04d') '.mat']; %%% NOTE THIS IS ASC
else
  iDorA
  error('need iDorA = +/- 1')
end

iDoX = +1;
if ~exist(fin)
  fprintf(1,'JOB %3i : avg era timeseries file %s DNE \n',JOB,fin);
  iDoX = -1;
end

iDo = +1;
if exist(fout)
  fprintf(1,'JOB %3i : avg era timeseries file %s already exists \n',JOB,fout);
  iDo = -1;
end

if iDo > 0 & iDoX > 0

  loader = ['load ' fin];
  eval(loader)
  p2 = pnew_op;
  %%%%%%%%%%%%%%%%%%%%%%%%%
  quants = [0 0.03 0.05 0.1 0.2 0.5 0.8 0.9 0.95 0.97 1.0];
  booQ = quantile(rad2bt(1231,p2.rcalc(1520,:)),quants)
  find_avg_Qprofiles

  globalavg  = convert_rtp_to_cloudOD(hnew_op,globalavg);
  globalQavg = convert_rtp_to_cloudOD(hnew_op,globalQavg);
  %%%%%%%%%%%%%%%%%%%%%%%%%
  find_avg_Qplots
  %%%%%%%%%%%%%%%%%%%%%%%%%
  %quick_jacs
  %%%%%%%%%%%%%%%%%%%%%%%%%

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
  if ~exist(fout)
    eval(saver)
  end

  %figure(3); scatter_coast(pnew_ip.rlon,pnew_ip.rlat,50,rad2bt(1231,pnew_ip.robs1(1520,:)))
  figure(14); scatter_coast(pnew_ip.rlon,pnew_ip.rlat,50,rad2bt(1231,pnew_ip.rcalc(1520,:)));            title('allsky calc');
  figure(15); scatter_coast(pnew_ip.rlon,pnew_ip.rlat,50,rad2bt(1231,pnew_ip.sarta_rclearcalc(1520,:))); title('clrsky calc');

  fprintf(1,'JOB = %4i fout = %s \n',JOB,fout)
end

disp('now do missinglist_whole_12month_ERA5.m  and then cluster_driver_put_together_globalavg_profiles.m and then master_driver_put_together_globalavg_profiles')
%stoptime = toc;
ticcEndF = clock;
%fprintf(1,'took %8.6f minutes to run \n',(stoptime-starttime)/60);
fprintf(1,'elapsed time : took %8.6f minutes to run \n',etime(ticcEndF,ticcStart0)/60);
monitor_memory_whos;
