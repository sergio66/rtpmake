%% have to hardcode  YYMMDD before starting, loops over granules

% local running to test
% clustcmd -L clust_make_ecmcloudrtp_sergio_sarta_YYMMDD_loopGG_2834.m 001:240
%
% otherwise when happy
% clustcmd -q medium -n 64 -p 4 clust_make_ecmcloudrtp_sergio_sarta_YYMMDD_loopGG_2834.m 001:240
%
% or
% clustcmd -q long_contrib -n 64 -p 4 clust_make_ecmcloudrtp_sergio_sarta_YYMMDD_loopGG_2834.m 001:240

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
JOB = str2num(getenv('SLURM_ARRAY_TASK_ID'));
%JOB = 215

%% March 11, 2011 is my good example
%yymmdd0  = input('Enter [YYYY MM DD] : ');
%iaGlist = input('Enter [GranStart GranList] : ');

yymmdd0  = [2011 03 11];  ddLoop = [];  %% for JPL set 1
yymmdd0  = [2011 07 11];  ddLoop = [];  %% for JPL set 2

yymmdd0  = [2018 06 29];  ddLoop = [];  %% for Ruben Delgado chesapeake data and testing IASI-NG
yymmdd0  = [2022 01 15];  ddLoop = [13 14 15 16 17 18 19 20 21 22];  %% <<<<<<<<<<<<<<<<<< CAN LOOP OVER DAYS >>>>>>>>>>>>>>>>>>>>>>    %% ECMWF says ATMS shows gravity waves from Tonga

yymmdd0  = [2018 10 31];  ddLoop = [];  %% for Alan Geer set 1
yymmdd0  = [2018 11 01];  ddLoop = [];  %% for Alan Geer set 2

if length(ddLoop) == 0
  ddLoop = yymmdd0(3);
end

for iiddloop = 1 : length(ddLoop)
  yymmdd0 = [yymmdd0(1:2) ddLoop(iiddloop)]

  iaGlist  = 001 : 240;
  iaGlist = iaGlist(JOB);
  
  iPertTCC = +1;  %% use tcc model 1 (best so far)
  iPertTCC = -1;  %% use default tcc in ECMWF  <<<<<<<<<<<<<<<<< DEFAULT >>>>>>>>>>>>>>
  
  iSlabCld_CumSumStrowORGeorge = +1; %% strow,  cumsum 9999, cloud at PEAK of wgt fcn <<<< DEFAULT >>>>>>>
  iSlabCld_CumSumStrowORGeorge = -1; %% aumann, cumsum -1,   cloud at mean of cld profile
  
  if iPertTCC <= 0
    cloud_set_defaults_run_maker
  else
    cloud_set_defaults_run_maker_pertTCC
  end
end
