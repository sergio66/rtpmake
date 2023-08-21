#!/bin/bash

# run this with sbatch -p high_mem --array=1-X --output='testslurm' sergio_matlab_jobB.sbatch
# run this with sbatch -p high_mem --array=1-Y sergio_matlab_jobB.sbatch 1 for filelist
# run this with sbatch -p high_mem --array=1-Z --output='/dev/null' sergio_matlab_jobB.sbatch 
# run this with sbatch -p high_mem --array=1-240 sergio_matlab_jobB.sbatch 2    for loop over a day
# run this with sbatch -p high_mem --array=1-415 sergio_matlab_jobB.sbatch 4    for tile center avg
# run this with sbatch -p high_mem --array=1-415 sergio_matlab_jobB.sbatch 5    for tile with 7x4=28 ERA grid points/tile, 16 days .. set quue to medium, time to 24 hours
# run this with sbatch -p high_mem --array=1-240 sergio_matlab_jobB.sbatch 6    for loop over a day, interp ERA
# run this with sbatch -p high_mem --array=1-240 sergio_matlab_jobB.sbatch 7    for randomly generate 240 granules spanning specifiled year

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
##SBATCH --qos=short+
##SBATCH --time=0:59:00 

# for command line argument selection 5
#SBATCH --qos=medium+
#SBATCH --time=23:59:00 

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
  echo "cmd line arg = 2, choose day, loop over all granules in the day, no interp"
  srun matlab -nodisplay -r "clust_make_eracloudrtp_sergio_sarta_YYMMDD_loopGG; exit"
elif [[ "$1" -eq "3" ]]; then
  echo "cmd line arg = 3, making uniform clear from allfov"
  srun matlab -nodisplay -r "clustmake_uniformclear_from_clearcalc; exit"
elif [[ "$1" -eq "4" ]]; then
  echo "cmd line arg = 4, making tile 16day avg"
  srun matlab -nodisplay -r "clust_loop_make16day_tile_center; exit"
elif [[ "$1" -eq "5" ]]; then
  echo "cmd line arg = 5, making tile 16day individiual days/ERA grid points"
  srun matlab -nodisplay -r "clust_loop_make16day_tile_28points; exit"
elif [[ "$1" -eq "6" ]]; then
  echo "cmd line arg = 6, making one whole day interp ERA, SYMBOLIC LINK to clustbatch_make_eracloudrtp_sergio_sarta_filelist_interp_YYMMDD_loopGG.m"
  srun matlab -nodisplay -r "clustbatch_eracloudrtp_sarta_filelist_interp_YYMMDD_loopGG; exit"
elif [[ "$1" -eq "7" ]]; then
  echo "cmd line arg = 7, making randowm data through year!!!"
  srun matlab -nodisplay -r "clustbatch_make_eracloudrtp_sarta_filelist_interp_YY_randMMDDGG; exit"
fi
