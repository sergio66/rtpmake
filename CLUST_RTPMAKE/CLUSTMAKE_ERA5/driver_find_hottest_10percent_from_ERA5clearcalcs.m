%% see clust_loop_make_monthly_tile_273points.m
addpath /asl/matlib/aslutil
addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/PLOTTER

run_sarta.cumsum = -1;    %% this is "closer" to MRO but since cliuds are at centroid, does not do too well with DCC  ... did this June 2024
run_sarta.cumsum = 9999;  %% larrabee likes this, puts clouds high so does well for DCC  try this for timestep 12-0 = Fig 8 in trends.tex

if run_sarta.cumsum == -1
  num_months = 240;  %% doing the complete timeseries, so we can do SKT trends
elseif run_sarta.cumsum == 9999
  num_months = 1;  %% just doing 2012/08 so we can compare to histogram Fig 8 in trends.tex
end

existence   = zeros(num_months,21,12);  %% Months x Xpts x Ypts;
if ~exist('loadedfiles')
  loadedfiles = zeros(num_months,21,12);     %% Months x Xpts x Ypts;
  bt1231clr   = nan(num_months,21,12,4608);  %% Months x Xpts x Ypts;
  bt1231cld   = nan(num_months,21,12,4608);  %% Months x Xpts x Ypts;
  stemp       = nan(num_months,21,12,4608);  %% Months x Xpts x Ypts;
  landfrac    = nan(num_months,21,12,4608);  %% Months x Xpts x Ypts;
end

find_hottest_10percent_from_ERA5clearcalc_progress   %% first plot the progress ... there should be about (21 x 12 = 252) mat files made per month
pause(1)

load_in_unreadfiles_10percent_from_ERA5clearcalc_progress
