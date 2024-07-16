#!/bin/bash

# run this with   sbatch -p cpu2021 --array=1-X --output='testslurm' sergio_matlab_jobB.sbatch
# run this with   sbatch -p cpu2021 --array=1-Z --output='/dev/null' sergio_matlab_jobB.sbatch 
# run this with   sbatch -p cpu2021 --array=1-Y   sergio_matlab_jobB.sbatch 1       for filelist, no interp
# run this with   sbatch -p cpu2021 --array=1-240 sergio_matlab_jobB.sbatch 2       for loop over a day, no interp
# run this with   sbatch -p cpu2021 --array=1-Y   sergio_matlab_jobB.sbatch 3       for filelist, with interp analysis
# run this with   sbatch -p cpu2021 --array=1-240 sergio_matlab_jobB.sbatch 4       for loop over a day, with interp analysis
# run this with   sbatch -p cpu2021 --array=1-Y   sergio_matlab_jobB.sbatch 5       for loop over a day, with interp analysis, but you can give a whole bunch of days thru a filelist
# run this with   sbatch -p cpu2021 --array=1-240 sergio_matlab_jobB.sbatch 4       for making uniform clear from allfov

#  Name of the job:
#SBATCH --job-name=MAKE_CRISRTP
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
##SBATCH --qos=medium+
##SBATCH --time=23:59:00 

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
  echo "cmd line arg = DNE, filelist run, no Interp"
  srun matlab -nodisplay -r "iInterp = 0; clustbatch_make_eraORecm_cloudrtp_sergio_sarta_filelist;exit"
elif [[ "$1" -eq "1" ]]; then
  echo "cmd line arg = 1, filelist run, no Interp"
  srun matlab -nodisplay -r "iInterp = 0; clustbatch_make_eraORecm_cloudrtp_sergio_sarta_filelist;exit"
elif [[ "$1" -eq "2" ]]; then
  echo "cmd line arg = 2, choose day, loop over all granules in the day, noInterp"
  srun matlab -nodisplay -r "iInterp = 0; clustbatch_make_eraORecm_cloudrtp_sergio_sarta_YYMMDD_loopGG; exit"
elif [[ "$1" -eq "3" ]]; then
  echo "cmd line arg = 3, filelist run with interp"
  srun matlab -nodisplay -r "iInterp = 1; clustbatch_make_eraORecm_cloudrtp_sergio_sarta_filelist;exit"
elif [[ "$1" -eq "4" ]]; then
  echo "cmd line arg = 4, choose day, loop over all granules in the day, yes Interp"
  srun matlab -nodisplay -r "iInterp = 1; clustbatch_make_eraORecm_cloudrtp_sergio_sarta_YYMMDD_loopGG; exit"
elif [[ "$1" -eq "5" ]]; then
  echo "cmd line arg = 5, filelist sets the day, loop over all granules in the day, yes Interp"
  #srun matlab -nodisplay -r "iInterp = 1; clustbatch_make_eraORecmcloudrtp_sergio_sarta_filelist_YYMMDD_loopGG; exit"   NAME TOO LONG
  srun matlab -nodisplay -r "iInterp = 1; clustbatch_make_eraORecmcloudrtp_filelist_YYMMDD_loopGG; exit"
elif [[ "$1" -eq "6" ]]; then
  echo "cmd line arg = 6, making uniform clear from allfov"
  srun matlab -nodisplay -r "iInterp = 0; clustmake_uniformclear_from_clearcalc; exit"
fi
