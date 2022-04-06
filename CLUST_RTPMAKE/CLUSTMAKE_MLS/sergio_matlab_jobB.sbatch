#!/bin/bash

# have 18 years of data (2004/09 to 2020/08) so 12x18 = 216
# run this with sbatch -p high_mem --array=1-216 sergio_matlab_jobB.sbatch 1          for mls

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
  echo "cmd line arg = DNE, mls filelist run"
  srun matlab -nodisplay -r "clust_compute_mls_profile_rtpfiles;exit"
elif [[ "$1" -eq "1" ]]; then
  echo "cmd line arg = 1, mls filelist run"
  srun matlab -nodisplay -r "clust_compute_mls_profile_rtpfiles;exit"
fi
