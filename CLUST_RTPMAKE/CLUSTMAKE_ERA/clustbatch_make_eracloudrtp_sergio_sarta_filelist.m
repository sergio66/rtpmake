%% run with
%% sbatch --array=N1-N2 --output='testslurm' sergio_matlab_jobB.sbatch
%% N1 = 1, N2 = number of files to be processed

%% specify text file which has YY MM DD GG lst that needs to be processed 
set_filelist
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

JOB = str2num(getenv('SLURM_ARRAY_TASK_ID'));
%JOB = 37
%JOB = 3

thefilelist = load(filelist);
thefilelist = thefilelist(JOB,:);

yymmdd0  = thefilelist(1:3); %% YY MM DD
iaGlist  = thefilelist(4);   %% granule

%yymmdd0 = [2016 01 18]; iaGlist = 180;   %% steve has made this. 
                                          %% looks same as mine (colwater = 50 mm water in tropics) 
                                          %% ../CLUSTMAKE_ERA5/ERAorECM/allfov_era_airicrad_day_2016018_180.rtp
%yymmdd0 = [2014 01 29]; iaGlist = 16; iaGlist = 32;  iaGlist = 31; %% crazy col water amount, 120 mm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iPertTCC = +1;  %% use tcc model 1 (best so far)
iPertTCC = -1;  %% use default tcc in ECMWF  <<<<<<<<<<<<<<<<< DEFAULT >>>>>>>>>>>>>>

iSlabCld_CumSumStrowORGeorge = +1; %% strow,  cumsum 9999, cloud at PEAK of wgt fcn <<<< DEFAULT >>>>>>>
iSlabCld_CumSumStrowORGeorge = -1; %% aumann, cumsum -1,   cloud at mean of cld profile

cloud_set_defaults_run_maker
