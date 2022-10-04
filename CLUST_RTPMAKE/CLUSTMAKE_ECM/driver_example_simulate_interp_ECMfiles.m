%% see eg ~/MATLABCODE/WetBulbTemperatures/make_usa_august_timeseries.m

iMake = +1;
iMake = -1;
if iMake > 0
  addpath /home/sergio/MATLABCODE/RTPMAKE/CLUST_RTPMAKE/CLUSTMAKE_ECM/
  %% look at MATLABCODE/RTPMAKE/CLUST_RTPMAKE/CLUSTMAKE_ECM/clustbatch_make_ecmcloudrtp_sergio_sarta_filelist_interp.m
  %% look at MATLABCODE/RTPMAKE/CLUST_RTPMAKE/CLUSTMAKE_ECM/cloud_set_defaults_run_maker_interp_analysis.m

  yymmdd0 = [2020 08 23];
  iaGlist = 241;
  iPertTCC = +1;  %% use tcc model 1 (best so far)
  iPertTCC = -1;  %% use default tcc in ECMWF  <<<<<<<<<<<<<<<<< DEFAULT >>>>>>>>>>>>>>

  iSlabCld_CumSumStrowORGeorge = +1; %% strow,  cumsum 9999, cloud at PEAK of wgt fcn <<<< DEFAULT >>>>>>>
  iSlabCld_CumSumStrowORGeorge = -1; %% aumann, cumsum -1,   cloud at mean of cld profile

  iTimeOffset = 000;     %% time offset in minutes
  cloud_set_defaults_run_maker_interp_analysis %% this interps the 00,06,12,18 analysis

  iTimeOffset = 060;     %% time offset in minutes
  cloud_set_defaults_run_maker_interp_analysis %% this interps the 00,06,12,18 analysis

  iTimeOffset = 120;     %% time offset in minutes
  cloud_set_defaults_run_maker_interp_analysis %% this interps the 00,06,12,18 analysis

  iTimeOffset = -60;     %% time offset in minutes
  cloud_set_defaults_run_maker_interp_analysis %% this interps the 00,06,12,18 analysis

  iTimeOffset = -120;     %% time offset in minutes
  cloud_set_defaults_run_maker_interp_analysis %% this interps the 00,06,12,18 analysis
end
