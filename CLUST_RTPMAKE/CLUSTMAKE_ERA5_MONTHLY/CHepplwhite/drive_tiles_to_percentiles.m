function drive_tiles_to_percentiles()
%
% Driver file for master script: proc_tiles_to_percentiles.m
% Given Latitude band (rzone) batch controls year.
% e.g. sbatch --array=0-17 run_tiles_to_percentiles.sh
%

% Set thisYear
all_years = [2003:2020];
thisYear = 2003;

% Set rzone (1:S.P. 64:N.P.)
thisZone = 33;

% Get the slurm array job and assign to year number [1:18].
sindex = str2num(getenv('SLURM_ARRAY_TASK_ID'));   % 0-19999
thisYear = all_years(sindex+1);
%thisSet = 1;

% Call the master script
disp(['Starting master script ' char(datetime('now'))]);
disp(['thisYear = ' num2str(thisYear)])

  proc_tiles_to_percentiles(thisYear, thisZone);

fprintf(1, '*** Task run end %s\n', char(datetime('now')));
