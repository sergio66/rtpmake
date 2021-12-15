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

yymmdd0  = [2018 10 31];   %% George/Alan
yymmdd0  = [2018 11 01];   %% George/Alan

iaGlist  = JOB;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iPertTCC = +1;  %% use tcc model 1 (best so far)
iPertTCC = 0;   %% do clear only!!!!
iPertTCC = -1;  %% use default tcc in ECMWF  <<<<<<<<<<<<<<<<< DEFAULT >>>>>>>>>>>>>>

iSlabCld_CumSumStrowORGeorge = +1; %% strow,  cumsum 9999, cloud at PEAK of wgt fcn <<<< DEFAULT >>>>>>>
iSlabCld_CumSumStrowORGeorge = -1; %% aumann, cumsum -1,   cloud at mean of cld profile

cloud_set_defaults_run_maker
