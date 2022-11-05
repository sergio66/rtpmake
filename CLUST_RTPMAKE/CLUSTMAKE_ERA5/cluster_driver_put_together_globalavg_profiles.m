addpath /home/sergio/MATLABCODE/CRODGERS_FAST_CLOUD/
addpath /asl/matlab2012/rtptoolsV201/
addpath /home/sergio/MATLABCODE/CRODGERS_FAST_CLOUD/

simulateYear = 2012;

JOB = str2num(getenv('SLURM_ARRAY_TASK_ID'));  %% JOB = 1 : 64
% JOB = 32

iOffSet = (JOB-1)*72;

for ii = 1 : 72
  fprintf(1,'.')
  iaProfile(ii) = ii + iOffSet;
  fname = ['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Center12months/DESC/' num2str(simulateYear,'%04d') '/era5_full12months_tile_center_' num2str(ii + iOffSet,'%04d') '.mat'];
  loader = ['load ' fname];
  eval(loader)
  [~,junk01] = subset_rtp_allcloudfields(hnew_op,globalavg,[],[],01);
  [~,junk02] = subset_rtp_allcloudfields(hnew_op,globalavg,[],[],02);
  [~,junk03] = subset_rtp_allcloudfields(hnew_op,globalavg,[],[],03);
  [~,junk04] = subset_rtp_allcloudfields(hnew_op,globalavg,[],[],04);
  [~,junk05] = subset_rtp_allcloudfields(hnew_op,globalavg,[],[],05);
  [~,junk06] = subset_rtp_allcloudfields(hnew_op,globalavg,[],[],06);
  [~,junk07] = subset_rtp_allcloudfields(hnew_op,globalavg,[],[],07);
  [~,junk08] = subset_rtp_allcloudfields(hnew_op,globalavg,[],[],08);
  [~,junk09] = subset_rtp_allcloudfields(hnew_op,globalavg,[],[],09);
  [~,junk10] = subset_rtp_allcloudfields(hnew_op,globalavg,[],[],10);
  [~,junk11] = subset_rtp_allcloudfields(hnew_op,globalavg,[],[],11);
  if ii == 1
    hout = hnew_op;
    prof01 = junk01;
    prof02 = junk02;
    prof03 = junk03;
    prof04 = junk04;
    prof05 = junk05;
    prof06 = junk06;
    prof07 = junk07;
    prof08 = junk08;
    prof09 = junk09;
    prof10 = junk10;
    prof11 = junk11;
  else
    [~,prof01] = cat_rtp(hnew_op,prof01,hnew_op,junk01);
    [~,prof02] = cat_rtp(hnew_op,prof02,hnew_op,junk02);
    [~,prof03] = cat_rtp(hnew_op,prof03,hnew_op,junk03);
    [~,prof04] = cat_rtp(hnew_op,prof04,hnew_op,junk04);
    [~,prof05] = cat_rtp(hnew_op,prof05,hnew_op,junk05);
    [~,prof06] = cat_rtp(hnew_op,prof06,hnew_op,junk06);
    [~,prof07] = cat_rtp(hnew_op,prof07,hnew_op,junk07);
    [~,prof08] = cat_rtp(hnew_op,prof08,hnew_op,junk08);
    [~,prof09] = cat_rtp(hnew_op,prof09,hnew_op,junk09);
    [~,prof10] = cat_rtp(hnew_op,prof10,hnew_op,junk10);
    [~,prof11] = cat_rtp(hnew_op,prof11,hnew_op,junk11);
  end
end

comment = 'see ~/MATLABCODE/RTPMAKE/CLUST_RTPMAKE/CLUSTMAKE_ERA5/cluster_driver_put_together_globalavg_profiles.m';
saver = ['save /asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Center12months/DESC/' num2str(simulateYear,'%04d') '/era5_full12months_latbin_' num2str(JOB,'%02d') '_tile_center_profilesQcumulative_1_11.mat hnew_op prof* comment'];
eval(saver)

