%% have to hardcode  YYMMDD before starting, loops over granules

% local running to test
% clustcmd -L clust_make_eracloudrtp_sergio_sarta_YYMMDD_loopGG.m 001:240
%
% otherwise when happy
% clustcmd -q medium -n 64 clust_make_eracloudrtp_sergio_sarta_YYMMDD_loopGG.m 001:240
%
% or
% clustcmd -q long_contrib -n 64 clust_make_eracloudrtp_sergio_sarta_YYMMDD_loopGG.m 001:240
%
% sbatch --array=1-240 sergio_matlab_jobB.sbatch

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

JOB = str2num(getenv('SLURM_ARRAY_TASK_ID'));

%% March 11, 2011 is my good example
%yymmdd0  = input('Enter [YYYY MM DD] : ');
%iaGlist = input('Enter [GranStart GranList] : ');

yymmdd0  = [2011 03 11];
yymmdd0  = [2013 08 27];   %% for CrIS high res, do AIRS for fun

yymmdd0  = [2010 08 01];   %% lotsa AIRS CO from Russian fires
yymmdd0  = [2010 08 14];   %% lotsa AIRS CO from Russian fires
yymmdd0  = [2010 07 21];   %% lotsa AIRS CO from Russian fires

yymmdd0  = [2018 09 06];   %% steve is testing uniform clear for CRIS so need AIRS, 208 day 249

yymmdd0  = [2017 01 17];   %% NLTE from Zhenglong .. January!
yymmdd0  = [2017 07 17];   %% NLTE from Zhenglong
yymmdd0  = [2018 06 29];   %% NLTE from Zhenglong, now I am puzzled why ERA is so bad compared to ECM

iaGlist  = JOB;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iPertTCC = +1;  %% use tcc model 1 (best so far)
iPertTCC = 0;   %% do clear only!!!!
iPertTCC = -1;  %% use default tcc in ECMWF  <<<<<<<<<<<<<<<<< DEFAULT >>>>>>>>>>>>>>

iSlabCld_CumSumStrowORGeorge = +1; %% strow,  cumsum 9999, cloud at PEAK of wgt fcn <<<< DEFAULT >>>>>>>
iSlabCld_CumSumStrowORGeorge = -1; %% aumann, cumsum -1,   cloud at mean of cld profile

cloud_set_defaults_run_maker
