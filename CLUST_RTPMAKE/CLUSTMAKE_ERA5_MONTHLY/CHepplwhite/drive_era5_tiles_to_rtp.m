% drive_era5_monthly_tiles_to_rtp
%
% test tile:cloudy [67,35] lonbin,latbin

addpath /home/chepplew/myLib/matlib

% Standard anomaly date range
sdate = '2002/09/01';
edate = '2019/08/31';

% Get the slurm array job and assign current job.
% Current RANGE: 72*64 = 4608

sindex = str2num(getenv('SLURM_ARRAY_TASK_ID'));
%sindex=3416;

% check valid range
if(sindex < 1 | sindex > 4608) 
  error('slurm task number out of range');
  return;
end

thisJob = sindex;
disp(['this job: ' num2str(thisJob)]);

% Lookup table for job configuration (one less than tile boundaries!)
lonbins = [1:72];
latbins = [1:64];
[X,Y] = meshgrid(lonbins,latbins);

% job=1 to start lonbin=1; latbin=1
% eg cloudy test tile: ix = find(X == 67 & Y == 35) -> 4529

lonbin  = X(thisJob);
latbin  = Y(thisJob);

disp(['Processing lonbin: ' num2str(lonbin) ' latbin: ' num2str(latbin)])

% call the load script
disp(' ***  calling load_era5_monthly *** ')
[era5] = load_era5_monthly(sdate, edate, lonbin, latbin);

% call the raw to tile value processor
disp(' ***  calling get_era5_tile_center_vals  ***')
[scen, smean, stdv] = get_era5_tile_center_vals(era5);

% call the rates calculator
disp(' ***  calling get_era5_tile_rates  *** ')
[xfits]   = get_era5_tile_rates(era5, scen,smean);


% call the rtp saver
disp(' ***  calling  write)era5_tile_vals_to_rtp  ***')
[iok] = write_era5_tile_vals_to_rtp(era5, scen);

iok
