%% see clust_loop_make_monthly_tile_273points.m
addpath /asl/matlib/aslutil
addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/PLOTTER

check_all_jobs_done('/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Monthly/summary',24,'.mat');

JOB = str2num(getenv('SLURM_ARRAY_TASK_ID'));  %% about 240 months, each JOB does 10 so 24 of these for 20 years
if length(JOB) == 0
  JOB = 24;
  JOB = 12;
  JOB = 001;
end

numJOBS = 10;

existence   = zeros(numJOBS,21,12);  %% Months x Xpts x Ypts;
if ~exist('loadedfiles')
  loadedfiles = zeros(numJOBS,21,12);     %% Months x Xpts x Ypts;
  bt1231clr   = nan(numJOBS,21,12,4608);  %% Months x Xpts x Ypts;
  bt1231cld   = nan(numJOBS,21,12,4608);  %% Months x Xpts x Ypts;
  stemp       = nan(numJOBS,21,12,4608);  %% Months x Xpts x Ypts;
  landfrac    = nan(numJOBS,21,12,4608);  %% Months x Xpts x Ypts;
end

%% remember a 3 x 5 CHIRP tile has about 240 (turns out to be 252) ERA5 grid points per tile
%% so there is /asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Monthly/ERAindex1 to /asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Monthly/ERAindex252
%% each of the 252 directories has 240 timesteps!!!!!
find_hottest_10percent_from_ERA5clearcalc_progress_JOB   %% first plot the progress ... there should be about (21 x 12 = 252) mat files made per month
pause(1)


%% load them in, process them (average and hottest 10 percent) 10 timesteps at a go, then save into /asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Monthly/summary_*
load_in_unreadfiles_10percent_from_ERA5clearcalc_progress_JOB

disp('when all 240/10 = 20 files done, run driver_load_in_summary_10percent_from_ERA5clearcalc.m')
disp('when all 240/10 = 20 files done, run driver_load_in_summary_10percent_from_ERA5clearcalc.m')
disp('when all 240/10 = 20 files done, run driver_load_in_summary_10percent_from_ERA5clearcalc.m')
