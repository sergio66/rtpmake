%% run with
%% sbatch --array=N1-N2 --output='testslurm' sergio_matlab_jobB.sbatch
%% sbatch --array=1-48 sergio_matlab_jobB.sbatch 
%% N1 = 1, N2 = number of files to be processed

%% note sergio_matlab_jobB.sbatch reminds you that
%% echo cmd line arg = 6, making one whole day interp ERA, SYMBOLIC LINK to clustbatch_make_eracloudrtp_sergio_sarta_filelist_interp_YYMMDD_loopGG.m
%% ie 
%% ls -lt clustbatch_eracloudrtp_sarta_filelist_interp_YYMMDD_loopGG.m
%%     clustbatch_eracloudrtp_sarta_filelist_interp_YYMMDD_loopGG.m -> clustbatch_make_eracloudrtp_sergio_sarta_filelist_interp_YYMMDD_loopGG.m

addpath /home/sergio/MATLABCODE
system_slurm_stats

JOB = str2num(getenv('SLURM_ARRAY_TASK_ID'));
%JOB = 9

warning('off', 'MATLAB:imagesci:hdfeos:removalWarningHDFSW');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% specify text file which has YY MM DD GG lst that needs to be processed
set_filelist

thefilelist = load(filelist);

yymmdd0  = [2011 03 11];
yymmdd0  = thefilelist(JOB,1:3);

iaGlist = 001 : 240;

%disp('problem with G151 on 2005/01/30')
%iaGlist = 152 : 240;
%iaGlist = 151 : 151


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iPertTCC = +1;  %% use tcc model 1 (best so far)
iPertTCC = -1;  %% use default tcc in ECMWF  <<<<<<<<<<<<<<<<< DEFAULT >>>>>>>>>>>>>>

iSlabCld_CumSumStrowORGeorge = +1; %% strow,  cumsum 9999, cloud at PEAK of wgt fcn <<<< DEFAULT >>>>>>>
iSlabCld_CumSumStrowORGeorge = -1; %% aumann, cumsum -1,   cloud at mean of cld profile

%% for ixx = 1 : length(iaGlist)

  if iPertTCC <= 0
    iTimeOffset = 060;     %% time offset in minutes
    iTimeOffset = 120;     %% time offset in minutes
    iTimeOffset = -060;    %% time offset in minutes
    iTimeOffset = -120;    %% time offset in minutes
  
    iTimeOffset = 000;     %% time offset in minutes
  
    cloud_set_defaults_run_maker_interp_analysis %% this interps the 00,06,12,18 analysis
  else
    error('enh wazzup doc, not yet possible to interp this : cloud_set_defaults_run_maker_pertTCC')
  end
  
%% end
