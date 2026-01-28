%% run with
%% sbatch --array=N1-N2 --output='testslurm' sergio_matlab_jobB.sbatch
%% sbatch --array=1-48 sergio_matlab_jobB.sbatch 
%% N1 = 1, N2 = number of files to be processed

JOB = str2num(getenv('SLURM_ARRAY_TASK_ID'));  %% JOB = 1--72 lonbins, then loop over the 10 files in set_filelist_chirptile.m
if length(JOB) == 0
  JOB = 1;
  JOB = 30;
  JOB = 24;
  JOB = 37;
end
warning('off', 'MATLAB:imagesci:hdfeos:removalWarningHDFSW');
fprintf(1,'JOB = %3i \n',JOB);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% specify text file which has YY MM DD GG lst that needs to be processed 
set_filelist_chirptile

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
len_dir0 = length(dir0);

%iaGlist = dir([dir0{JOB} '/*.nc']);   %% orig ... choose a dir0{JOB},  loop over 1:72 lonbins per JOB
iaGlist = 1 : 72;                     %% new ...  chose a lonbin_JOB), loop over thedir0 list

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iPertTCC = +1;  %% use tcc model 1 (best so far)
iPertTCC = -1;  %% use default tcc in ECMWF  <<<<<<<<<<<<<<<<< DEFAULT >>>>>>>>>>>>>>

iSlabCld_CumSumStrowORGeorge = +1; %% strow,  cumsum 9999, cloud at PEAK of wgt fcn <<<< DEFAULT >>>>>>>
iSlabCld_CumSumStrowORGeorge = -1; %% aumann, cumsum -1,   cloud at mean of cld profile

if iPertTCC <= 0
  iTimeOffset = 060;     %% time offset in minutes
  iTimeOffset = 120;     %% time offset in minutes
  iTimeOffset = -060;    %% time offset in minutes
  iTimeOffset = -120;    %% time offset in minutes

  iTimeOffset = 000;     %% time offset in minutes

  cloud_set_defaults_run_maker_chirptiles %% this interps the 00,06,12,18 analysis
else
  error('enh wazzup doc, not yet possible to do this')
end

