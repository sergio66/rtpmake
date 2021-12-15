function drive_tiles_to_daily()
%
% Driver file for master script: proc_tiles_to_daily.m
% Given the Zone, batch drives on year number
% e.g. sbatch --array=0-17 run_iles_to_daily.sh
%

% Set thisYear
all_years = [2003:2020];
thisYear = 2004;

% Set rzone (1:S.P. 64:N.P.)
thisZone = 1;

% Get the slurm array job and assign to year nummber [1:18].
sindex = str2num(getenv('SLURM_ARRAY_TASK_ID'));   % 0-19999
thisYear = all_years(sindex+1);
%thisYear = all_years(1);

% Call the master script
disp(['Starting master script ' char(datetime('now'))]);
disp(['thisYear = ' num2str(thisYear)])

  proc_tiles_to_daily(thisYear, thisZone);

fprintf(1, '*** Task run end %s\n', char(datetime('now')));
