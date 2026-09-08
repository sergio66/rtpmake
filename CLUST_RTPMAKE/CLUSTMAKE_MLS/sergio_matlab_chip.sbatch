#!/bin/bash

# have 18 years of data (2002/09 to 2020/08) so 12x18 = 216
# have 20 years of data (2002/09 to 2022/08) so 12x20 = 240
# have 23 years of data (2002/09 to 2025/08) so 12x23 = 276
#
# run this with sbatch -p cpu2024  --array=1-216 sergio_matlab_chip.sbatch 1          for mls 2004-2022
# run this with sbatch -p cpu2024  --array=1-252 sergio_matlab_chip.sbatch 1          for mls 2004-2025
#    or after running sqfail
#               sbatch -p high_mem --exclude=cnode006,cnode007,cnode008,cnode027,cnode028 --array=1-216 sergio_matlab_chip.sbatch 1   for mls; 

########################################################################

#SBATCH --job-name=CLUST_MAKE_ACMIP6_RTP    #  Name of the job:
#SBATCH -N1                                 #  N specifies that 1 job step is to be allocated per matlab instance
#SBATCH --cpus-per-task 1                   #  number of cores per matlab session available for parallel jobs
##SBATCH --partition=2024                   # desired partition
#SBATCH --partition=match                   # desired partition
#SBATCH --cluster=chip-cpu                  # desired cluster
#SBATCH --account=pi_strow
#SBATCH --qos=shared                        # qos to get as many cpu2024 as possible, else put pi_strow

#SBATCH --mem=12000                         #  This is in MB, very aggressive but I have been running outta memory
#SBATCH --time=0:59:00    

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
  srun /usr/ebuild/installs/software/MATLAB/2023b/bin/matlab -nodisplay -r "clust_compute_mls_profile_rtpfiles;exit"
elif [[ "$1" -eq "1" ]]; then
  echo "cmd line arg = 1, mls filelist run"
  srun /usr/ebuild/installs/software/MATLAB/2023b/bin/matlab -nodisplay -r "clust_compute_mls_profile_rtpfiles;exit"
fi
