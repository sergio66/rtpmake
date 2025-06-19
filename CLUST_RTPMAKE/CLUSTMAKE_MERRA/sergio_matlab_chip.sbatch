#!/bin/bash

# no run not done  this with sbatch -p 2021 --array=1-X --output='testslurm' sergio_matlab_chip.sbatch
# no run not done  this with sbatch -p 2021 --array=1-Y sergio_matlab_chip.sbatch 1 for filelist
# no run not done  this with sbatch -p 2021 --array=1-Z --output='/dev/null' sergio_matlab_chip.sbatch 
# no run not done  this with sbatch -p 2021 --array=1-240 sergio_matlab_chip.sbatch 2    for loop over a day
# yes          run this with sbatch -p 2021 --array=1-228 sergio_matlab_chip.sbatch 4    for NERRA5 tile center avg, desc, 19 yearsx12months/yr = 228
# no run not done  this with sbatch -p 2021 --array=1-228 sergio_matlab_chip.sbatch 5    for tile with 7x4=28 MERRA grid points/tile, 16 days .. set quue to medium, time to 24 hours
# no run not done  this with sbatch -p 2021 --array=1-240 sergio_matlab_chip.sbatch 6    for clust_make_monthlyavg_2m.m
#                            sbatch -p 2021 --array=1-240 sergio_matlab_chip.sbatch 7    for clustbatch_make_merra_surfaceT.m

#  Name of the job:
#SBATCH --job-name=CLUST_MAKE_MERRA_RTP

# requeue; if held type scontrol release JOBID
#SBATCH --requeue

#  N specifies that 1 job step is to be allocated per instance of matlab
#SBATCH -N1

#  This specifies the number of cores per matlab session will be
#available for parallel jobs
#SBATCH --cpus-per-task 1

#SBATCH --account=pi_strow

#  Specify the desired partition develop/batch/prod
##SBATCH --partition=batch
## can shorthand this as -p high_mem on the command line
##SBATCH --partition=2018

#SBATCH --partition=2021
#SBATCH --cluster=chip-cpu

########################################################################
#  Specify the qos and run time (format:  dd-hh:mm:ss)

##SBATCH --qos=medium
##SBATCH --time=3:59:00

##SBATCH --qos=medium
##SBATCH --time=1:59:00

#SBATCH --qos=short
#SBATCH --time=0:59:00

#  This is in MB, very aggressive but I have been running outta memory
#SBATCH --mem-per-cpu=12000

########################################################################
### check the number of input arguments
if [ $# -gt 0 ]; then
  echo "Your command line contains $# arguments"
elif [ $# -eq 0 ]; then
  echo "Your command line contains no arguments"
fi
########################################################################

if [[ "$1" -eq "" ]]; then
  echo "cmd line arg = DNE, filelist run"
  #srun /usr/ebuild/installs/software/MATLAB/2023b/bin/matlab -nodisplay -r "clust_make_merracloudrtp_sergio_DCC_filelist; exit"
  srun /usr/ebuild/installs/software/MATLAB/2023b/bin/matlab -nodisplay -r "clustbatch_make_merracloudrtp_sergio_sarta_filelist;exit"
elif [[ "$1" -eq "1" ]]; then
  echo "cmd line arg = 1, filelist run"
  #srun /usr/ebuild/installs/software/MATLAB/2023b/bin/matlab -nodisplay -r "clust_make_merracloudrtp_sergio_DCC_filelist; exit"
  srun /usr/ebuild/installs/software/MATLAB/2023b/bin/matlab -nodisplay -r "clustbatch_make_merracloudrtp_sergio_sarta_filelist;exit"
elif [[ "$1" -eq "2" ]]; then
  echo "cmd line arg = 2, choose day, loop over all granules in the day"
  srun /usr/ebuild/installs/software/MATLAB/2023b/bin/matlab -nodisplay -r "clust_make_merracloudrtp_sergio_sarta_YYMMDD_loopGG; exit"
elif [[ "$1" -eq "3" ]]; then
  echo "cmd line arg = 3, making uniform clear from allfov"
  srun /usr/ebuild/installs/software/MATLAB/2023b/bin/matlab -nodisplay -r "clustmake_uniformclear_from_clearcalc; exit"
elif [[ "$1" -eq "4" ]]; then
  echo "cmd line arg = 4, making tile MERRA5 monthly avg"
  srun /usr/ebuild/installs/software/MATLAB/2023b/bin/matlab -nodisplay -r "clust_loop_make_monthly_tile_center_asc_or_desc; exit"
  #srun /usr/ebuild/installs/software/MATLAB/2023b/bin/matlab -nodisplay -r "clust_loop_make_monthly2m_tile_center_asc_or_desc; exit"
elif [[ "$1" -eq "5" ]]; then
  echo "cmd line arg = 5, making tile MERRA5 monthly individiual days/MERRA grid points"
  srun /usr/ebuild/installs/software/MATLAB/2023b/bin/matlab -nodisplay -r "clust_loop_make_monthly_tile_28points; exit"
elif [[ "$1" -eq "6" ]]; then
  echo "cmd line arg = 6, making tile MERRA5 monthly averaged 2m data"
  srun /usr/ebuild/installs/software/MATLAB/2023b/bin/matlab -nodisplay -r "clust_make_monthlyavg_2m; exit"
elif [[ "$1" -eq "7" ]]; then
  echo "cmd line arg = 7, making tile clustbatch_make_merra_surfaceT.m for trends paper"
  srun /usr/ebuild/installs/software/MATLAB/2023b/bin/matlab -nodisplay -r "clustbatch_make_merra2_monthly_surfaceT; exit"
fi
