%% run with
%% sbatch --array=N1-N2 --output='testslurm' sergio_matlab_jobB.sbatch
%% sbatch --array=1-48 sergio_matlab_jobB.sbatch 
%% N1 = 1, N2 = number of files to be processed

%% this is for TESTING ERIC MADY AI using a randomly made set of FOVS

addpath /home/sergio/MATLABCODE
system_slurm_stats

JOB = str2num(getenv('SLURM_ARRAY_TASK_ID'));    %% this is irrelevent since I assign MM DD GG randomly

warning('off', 'MATLAB:imagesci:hdfeos:removalWarningHDFSW');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% specify text file which has YY MM DD GG list that needs to be processed
%% set_filelist

%% thefilelist = load(filelist);

rng('shuffle')
pause(floor(JOB/10))
anana = clock; rng(ceil(anana(6)*10000000));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yy = 2004;  %%%%%%%%%%%%%%%%%%%%%%%%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
yy = 2005;  %%%%%%%%%%%%%%%%%%%%%%%%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mm = floor(rand(1,1)*12); 
if mm < 1
  mm = 1;
elseif mm > 12
  mm = 12;
end

daysINmonth = [31 28 31 30 31 30 31 31 30 31 30 31];
if mod(yy,4) == 0
  daysINmonth(2) = 29;
end
daysINmonth = daysINmonth(mm);
dd = floor(rand(1,1)*daysINmonth); 
if dd < 1
  dd = 1;
elseif dd > daysINmonth
  dd = daysINmonth;
end

gg = floor(rand(1,1)*240); 
if gg < 1
  gg = 1;
elseif gg > 240
  gg = 240;
end

yymmdd0  = [yy mm dd];
iaGlist  = gg;

fprintf(1,'will be making data for this random MM YY GG = %4i/%2i/%2i G%3i \n',[yymmdd0 iaGlist])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iPertTCC = +1;  %% use tcc model 1 (best so far)
iPertTCC = -1;  %% use default tcc in ECMWF  <<<<<<<<<<<<<<<<< DEFAULT >>>>>>>>>>>>>>

iSlabCld_CumSumStrowORGeorge = +1; %% strow,  cumsum 9999, cloud at PEAK of wgt fcn <<<< DEFAULT >>>>>>>
iSlabCld_CumSumStrowORGeorge = -1; %% aumann, cumsum -1,   cloud at mean of cld profile

iRandomMMDDGG = +1;

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

