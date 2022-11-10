JOB0 = str2num(getenv('SLURM_ARRAY_TASK_ID'));  %% JOB0 = 1 : 64
% JOB0 = 1

if JOB0 > 64
  error('JOB0 = 1 .. 64')
end

addpath /home/sergio/MATLABCODE/TIME
addpath /asl/matlib/rtptools/
addpath /asl/matlib/aslutil
addpath /asl/matlib/h4tools
addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/TIME
addpath /home/sergio/MATLABCODE/PLOTTER
addpath ../GRIB
addpath /home/sergio/MATLABCODE/matlib/clouds/sarta
addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS
addpath /home/sergio/MATLABCODE/CRODGERS_FAST_CLOUD
addpath /home/sergio/MATLABCODE/COLORMAP

system_slurm_stats

addpath /asl/matlib/aslutil
addpath /asl/packages/time
addpath /asl/matlib/science

for ii = 1 : 72
  JOB = ii + (JOB0-1)*72;
  loop_fixNaNbug_2012_average_whole_12month_ERA5
end
