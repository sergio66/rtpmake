#!/bin/bash

# have 12 years of data (2002/09 to 2014/08) so 12x12 = 144
# run this with sbatch -p high_mem --array=1-144 sergio_matlab_jobB.sbatch 1          for cmip6
# run this with sbatch -p cpu2021  --array=1-144 sergio_matlab_jobB.sbatch 2          for amip6

#  Name of the job:
#SBATCH --job-name=CLUST_MAKE_ACMIP6_RTP
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
  echo "cmd line arg = DNE, cmip6 filelist run"
  srun matlab -nodisplay -r "clust_compute_cmip6_profile_rtpfiles;exit"
elif [[ "$1" -eq "1" ]]; then
  echo "cmd line arg = 1, cmip6 filelist run"
  srun matlab -nodisplay -r "clust_compute_cmip6_profile_rtpfiles;exit"
elif [[ "$1" -eq "2" ]]; then
  echo "cmd line arg = 1, amip6 filelist run"
  srun matlab -nodisplay -r "clust_compute_amip6_profile_rtpfiles;exit"
fi
