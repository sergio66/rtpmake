#!/bin/bash

# run this with   sbatch --exclude=cnode[225,267] --array=1-48  --output='testslurm' sergio_matlab_chip.sbatch
# run this with   sbatch --exclude=cnode[225,267] --array=1-48  sergio_matlab_chip.sbatch 1   for filelist
# run this with   sbatch --exclude=cnode[225,267] --array=1-240 sergio_matlab_chip.sbatch 2   for loop over a day
# run this with   sbatch -p 2021                  --array=1-48  sergio_matlab_chip.sbatch 1   for filelist
# run this with   sbatch -p 2021                  --array=1-240 sergio_matlab_chip.sbatch 2   for loop over a day
# run this with   sbatch -p 2021                  --array=1-240 sergio_matlab_chip.sbatch 3   for filelist, with interp analysis
# run this with   sbatch -p 2021                  --array=1-72  sergio_matlab_chip.sbatch 4   for chirp tilelist

# >>>>>>>>>>>>>>>>>>>>>>>>> ########################################################################

#  Name of the job:
#SBATCH --job-name=CLUST_MAKE_ECM_RTP

# requeue; if held type scontrol release JOBID
#SBATCH --requeue

#  N specifies that 1 job step is to be allocated per instance of matlab
#SBATCH -N1

#  This specifies the number of cores per matlab session will be
#available for parallel jobs
#SBATCH --cpus-per-task 1

#SBATCH --account=pi_strow

#SBATCH --cluster=chip-cpu

########################################################################
#  Specify the qos and run time (format:  dd-hh:mm:ss)

## 100 layer clouds
##SBATCH --qos=medium
##SBATCH --time=3:59:00

# TwoSlab clouds jacs rads, loop 
##SBATCH --qos=medium
##SBATCH --time=3:59:00

# clear, or TwoSlab clouds, or fluxes
##SBATCH --qos=short
##SBATCH --time=0:59:00

# tiles since we loop over 10 or so tiles for each of the 72 lonbins
##SBATCH --time=7:59:00
#SBATCH --time=23:59:00

#########################

#  Specify the desired partition develop/batch/prod
##SBATCH --partition=batch
## can shorthand this as -p high_mem on the command line
##SBATCH --partition=2018
##SBATCH --partition=2021
##SBATCH --partition=2024

#### to make sure we have exclusive privilege
##SBATCH --qos=pi_strow
##SBATCH --partition=pi_strow

#### we run on 2024 but run risk of someone pre-empting
#SBATCH --qos=shared
#SBATCH --partition=2024

########################################################################

#### https://hpc-community.unige.ch/t/srun-fatal-slurm-mem-per-cpu-slurm-mem-per-gpu-and-slurm-mem-per-node-are-mutually-exclusive/722/5
##  This is in MB, less aggressive
#SBATCH --mem=16000

########################################################################

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
elif [[ "$1" -eq "4" ]]; then
  echo "cmd line arg = 4, tile filelist run"
  srun matlab -nodisplay -r "clustbatch_make_ecmcloudrtp_sergio_sarta_filelist_chirptiles; exit"
fi
