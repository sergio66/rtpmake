%% run with
%% sbatch --array=N1-N2 --output='testslurm' sergio_matlab_jobB.sbatch
%% sbatch --array=1-48 sergio_matlab_jobB.sbatch 
%% N1 = 1, N2 = number of files to be processed

JOB = str2num(getenv('SLURM_ARRAY_TASK_ID'));
%JOB = 1
warning('off', 'MATLAB:imagesci:hdfeos:removalWarningHDFSW');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% specify text file which has YY MM DD GG lst that needs to be processed 
set_filelist

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

thefilelist = load(filelist);
thefilelist = thefilelist(JOB,:)

yymmdd0  = thefilelist(1:3); %% YY MM DD
iaGlist  = thefilelist(4);   %% granule

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iPertTCC = +1;  %% use tcc model 1 (best so far)
iPertTCC = -1;  %% use default tcc in ECMWF  <<<<<<<<<<<<<<<<< DEFAULT >>>>>>>>>>>>>>

iSlabCld_CumSumStrowORGeorge = +1; %% strow,  cumsum 9999, cloud at PEAK of wgt fcn <<<< DEFAULT >>>>>>>
iSlabCld_CumSumStrowORGeorge = -1; %% aumann, cumsum -1,   cloud at mean of cld profile

if iPertTCC <= 0
  cloud_set_defaults_run_maker  %% this is just the nearest ECMWF file, whether 00,06,12,18 analysis or 03,09,15,21 forecast
  %% cloud_set_defaults_run_maker_interp_analysis %% this interps the 00,06,12,18 analysis
else
  cloud_set_defaults_run_maker_pertTCC
end

