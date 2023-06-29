#!/bin/bash

# run this with   sbatch --exclude=cnode[225,267] --array=1-48  --output='testslurm' sergio_matlab_jobB.sbatch
# run this with   sbatch --exclude=cnode[225,267] --array=1-48  sergio_matlab_jobB.sbatch 1   for filelist
# run this with   sbatch --exclude=cnode[225,267] --array=1-240 sergio_matlab_jobB.sbatch 2   for loop over a day
# run this with   sbatch -p cpu2021               --array=1-48  sergio_matlab_jobB.sbatch 1   for filelist
# run this with   sbatch -p cpu2021               --array=1-240 sergio_matlab_jobB.sbatch 2   for loop over a day
# run this with   sbatch -p cpu2021               --array=1-240 sergio_matlab_jobB.sbatch 3   for filelist, with interp analysis

#  Name of the job:
#SBATCH --job-name=CLUST_MAKE_ECM_RTP

#  N specifies that 1 job step is to be allocated per instance of
#matlab
#SBATCH -N1

#  This specifies the number of cores per matlab session will be
#available for parallel jobs
#SBATCH --cpus-per-task 1

#  Specify the desired partition develop/batch/prod
##SBATCH --partition=batch
#SBATCH --partition=high_mem

#  Specify the qos and run time (format:  dd-hh:mm:ss)
##SBATCH --qos=medium
##SBATCH --time=3:59:00 

#SBATCH --qos=short+
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
  srun matlab -nodisplay -r "clustbatch_make_ecmcloudrtp_sergio_sarta_filelist;exit"
elif [[ "$1" -eq "1" ]]; then
  echo "cmd line arg = 1, filelist run"
  srun matlab -nodisplay -r "clustbatch_make_ecmcloudrtp_sergio_sarta_filelist;exit"
elif [[ "$1" -eq "2" ]]; then
  echo "cmd line arg = 2, loop over day"
  srun matlab -nodisplay -r "clust_make_ecmcloudrtp_sergio_sarta_YYMMDD_loopGG_2834; exit"
elif [[ "$1" -eq "3" ]]; then
  echo "cmd line arg = 3, filelist run, but interp over EMC 3hour files"
  srun matlab -nodisplay -r "clustbatch_make_ecmcloudrtp_sergio_sarta_filelist_interp; exit"
fi
