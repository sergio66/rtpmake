#!/bin/bash

# run this with sbatch --array=1-48 --output='testslurm' sergio_matlab_jobB.sbatch 

#  Name of the job:
#SBATCH --job-name=CLUST_MAKERTP
#  N specifies that 1 job step is to be allocated per instance of
#matlab
#SBATCH -N1
#  This specifies the number of cores per matlab session will be
#available for parallel jobs
#SBATCH --cpus-per-task 1
#  Specify the desired partition develop/batch/prod
#SBATCH --partition=batch
#  Specify the qos and run time (format:  dd-hh:mm:ss)
#SBATCH --qos=medium
#SBATCH --time=3:59:00 
#  This is in MB, very aggressive but I have been running outta memory
#SBATCH --mem-per-cpu=12000
#  This forces use of hpcf2013 nodes add this to your sbatch script
#SBATCH --constraint=hpcf2013

#srun matlab -nodisplay -r "clustbatch_make_ecmcloudrtp_sergio_sarta_filelist; exit" 
#srun matlab -nodisplay -r "clustbatch_make_eracloudrtp_sergio_sarta_filelist; exit"
#srun matlab -nodisplay -r "clust_make_eraclearrtp_sergio_sarta_YYMMDD_loopGG; exit"
#srun matlab -nodisplay -r "clust_make_ecmclearrtp_sergio_sarta_YYMMDD_loopGG; exit"

srun matlab -nodisplay -r "clustbatch_make_eracloudrtp_nadir_sarta_filelist; exit"
