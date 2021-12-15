%% have to hardcode  YYMMDD before starting, loops over granules

%% run with
%% sbatch --array=N1-N2 --output='testslurm' sergio_matlab_jobB.sbatch
%% sbatch --array=1-48 sergio_matlab_jobB.sbatch 
%% N1 = 1, N2 = number of files to be processed

addpath /home/sergio/MATLABCODE

JOB = str2num(getenv('SLURM_ARRAY_TASK_ID'));
warning('off', 'MATLAB:imagesci:hdfeos:removalWarningHDFSW');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% March 11, 2011 is my good example
%yymmdd0  = input('Enter [YYYY MM DD] : ');
%iaGlist = input('Enter [GranStart GranList] : ');

yymmdd0  = [2011 03 11];   %% for JGR paper

yymmdd0  = [2011 01 11];   %% for JGR paper
yymmdd0  = [2011 06 12];   %% for JGR paper <<< 2011/06/11 is BAD
yymmdd0  = [2011 07 11];   %% for JGR paper

yymmdd0  = [2012 09 20];   %% for SNO with CrIS, did with OLD landfrac, and with new one (oDec2013)
yymmdd0  = [2013 10 25];   %% for Ruben/John co-located O3 profiles over MD, job = 71,182

yymmdd0  = [2009 06 01];   %% for work with Ali/George
yymmdd0  = [2009 03 13];   %% for work with Ali/George

yymmdd0  = [2006 08 03];   %% for Zhibo, grans 63.64.80

yymmdd0  = [2006 09 03];   %% for paul testing
yymmdd0  = [2006 09 04];   %% for paul testing
yymmdd0  = [2006 09 05];   %% for paul testing

yymmdd0  = [2013 06 22];   %% for Andy/Chris OCO

yymmdd0  = [2008 08 24];   %% for Chenxi I showed two MODIS (Aqua) land granules, 
                           %% which are MYD2008237.1955 and MYD2005091.1315

yymmdd0  = [2011 03 11];   %% for JGR paper
yymmdd0  = [2013 08 27];   %% for CrIS high res, do AIRS for fun
yymmdd0  = [2014 12 05];   %% for CrIS high res, do AIRS for fun

yymmdd0  = [2014 02 08];   %% atmospheric river for March 2016 AIRS STM 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iaGlist  = JOB;

cloud_set_defaults_run_maker

