#!/bin/bash

# run this with sbatch -p high_mem --array=1-X --output='testslurm' sergio_matlab_jobB.sbatch
# run this with sbatch -p high_mem --array=1-Y sergio_matlab_jobB.sbatch 1 for filelist
# run this with sbatch -p high_mem --array=1-Z --output='/dev/null' sergio_matlab_jobB.sbatch 
# run this with sbatch -p high_mem --array=1-240  sergio_matlab_jobB.sbatch 2    for loop over a day
# run this with sbatch -p high_mem --array=1-228  sergio_matlab_jobB.sbatch 4    for tile center avg, desc, 19 yearsx12months/yr = 228
# run this with sbatch -p high_mem --array=1-228  sergio_matlab_jobB.sbatch 5    for tile with 7x4=28 ERA grid points/tile, 16 days .. set quue to medium, time to 24 hours
# run this with sbatch -p high_mem --array=1-240  sergio_matlab_jobB.sbatch 6    for clust_make_monthlyavg_2m.m
#
# run this with sbatch -p high_mem --array=1-4608 sergio_matlab_jobB.sbatch 7    for clustmake_whole_12month_ERA5.m
# run this with sbatch -p high_mem --array=1-64   sergio_matlab_jobB.sbatch 8    for cluster_driver_put_together_globalavg_profiles.m


#  Name of the job:
#SBATCH --job-name=CLUST_MAKE_ERA_RTP
#  N specifies that 1 job step is to be allocated per instance of
#matlab
#SBATCH -N1

#  This specifies the number of cores per matlab session will be
#available for parallel jobs
#SBATCH --cpus-per-task 1

#  Specify the desired partition develop/batch/prod
#SBATCH --partition=batch

#  Specify the qos and run time (format:  dd-hh:mm:ss)
# for all command line arguments 0,1,2,3,4 except selection 5
#SBATCH --qos=short+
#SBATCH --time=0:59:00 

# for command line argument selection 5
##SBATCH --qos=medium+
##SBATCH --time=23:59:00 

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
  #srun matlab -nodisplay -r "clust_make_eracloudrtp_sergio_DCC_filelist; exit"
  srun matlab -nodisplay -r "clustbatch_make_eracloudrtp_sergio_sarta_filelist;exit"
elif [[ "$1" -eq "1" ]]; then
  echo "cmd line arg = 1, filelist run"
  #srun matlab -nodisplay -r "clust_make_eracloudrtp_sergio_DCC_filelist; exit"
  srun matlab -nodisplay -r "clustbatch_make_eracloudrtp_sergio_sarta_filelist;exit"
elif [[ "$1" -eq "2" ]]; then
  echo "cmd line arg = 2, choose day, loop over all granules in the day"
  srun matlab -nodisplay -r "clust_make_eracloudrtp_sergio_sarta_YYMMDD_loopGG; exit"
elif [[ "$1" -eq "3" ]]; then
  echo "cmd line arg = 3, making uniform clear from allfov"
  srun matlab -nodisplay -r "clustmake_uniformclear_from_clearcalc; exit"
elif [[ "$1" -eq "4" ]]; then
  echo "cmd line arg = 4, making tile ERA5 monthly avg"
  srun matlab -nodisplay -r "clust_loop_make_monthly_tile_center_asc_or_desc; exit"
  #only have about months 007 to 019 of this
  #srun matlab -nodisplay -r "clust_loop_make_monthly2m_tile_center_asc_or_desc; exit"
elif [[ "$1" -eq "5" ]]; then
  echo "cmd line arg = 5, making tile ERA5 monthly individiual days/ERA grid points"
  srun matlab -nodisplay -r "clust_loop_make_monthly_tile_28points; exit"
elif [[ "$1" -eq "6" ]]; then
  echo "cmd line arg = 6, making tile ERA5 monthly averaged 2m data"
  srun matlab -nodisplay -r "clust_make_monthlyavg_2m; exit"
elif [[ "$1" -eq "7" ]]; then
  echo "cmd line arg = 7, making tile clustmake_whole_12month_ERA5.m"
  srun matlab -nodisplay -r "clustmake_whole_12month_ERA5; exit"
elif [[ "$1" -eq "8" ]]; then
  echo "cmd line arg = 8, making tile cluster_driver_put_together_globalavg_profiles"
  srun matlab -nodisplay -r "cluster_driver_put_together_globalavg_profiles; exit"
fi
