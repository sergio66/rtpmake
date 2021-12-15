function drive_era_interim_rates()

%
% Driver file for mast script: load_era_interim.m
% Given date range and airs tile lon/lat pair.
% e.g. sbatch --array=0-72 run_era_interim_rates.sh
%

addpath /home/chepplew/myLib/matlib       % load_era_interim

% Define place to save results
dhome = '/home/chepplew/data/rates_anomalies/tiled/era_interim/';
if(~exist(dhome)) mkdir(dhome); end

% Set the start and end date (use defaults)
date1 = [];
date2 = [];

% Get the slurm array job and assign to lonBin for fixed latBin.
sindex  = str2num(getenv('SLURM_ARRAY_TASK_ID'));   % 0-19999
thisBin = sindex+1;
if(~ismember(thisBin,[1:72]));     % tile boundary + 1
  error('Bin number out of range')
  return;
end
%thisBin = 64;

% Set the AIRS.tile by {lon,lat}Bin
% latBin can run from 2:64
lonBin = 64;
lonBin = thisBin;
latBin = 34;

% file name to save:
savfn = [dhome 'lon' sprintf('%02d',lonBin) '_lat' ...
         sprintf('%02d',latBin) '.mat'];

% Call the master script
disp(['Starting master script ' char(datetime('now'))]);
disp(['lon,lat Bins = ' num2str(lonBin) ', ' num2str(latBin)])

  [era,fit] = load_era_interim([],[],lonBin,latBin);

disp(['saving data: ' savfn])
  save(savfn,'era','fit');
fprintf(1, '*** Task run end %s\n', char(datetime('now')));
