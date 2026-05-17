%% have to hardcode  YYMMDD before starting, loops over granules
%% gets in few extra vars, such as vertical velocity, divergence, potential vorticity

%% see http://www.ecmwf.int/products/catalogue/I.html

%% sbatch --array=N1-N2 --output='testslurm' sergio_matlab_jobB.sbatch
%% sbatch --array=1-48 sergio_matlab_jobB.sbatch 
%% N1 = 1, N2 = number of files to be processed

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

JOB = str2num(getenv('SLURM_ARRAY_TASK_ID'));
if length(JOB) == 0
  JOB = 2;
  JOB = 3;    
  JOB = 1;
  JOB = 200;
  JOB = 205;
  JOB = 195;  
end

warning('off', 'MATLAB:imagesci:hdfeos:removalWarningHDFSW');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% March 11, 2011 is my good example
%yymmdd0  = input('Enter [YYYY MM DD] : ');
%iaGlist = input('Enter [GranStart GranList] : ');

yymmdd0 = [2012 09 20];   %% for SNO with CrIS

yymmdd0 = [2011 01 11];   %% for JGR paper
yymmdd0 = [2011 06 12];   %% for JGR paper <<< 2011/06/11 is BAD
yymmdd0 = [2011 07 11];   %% for JGR paper
yymmdd0 = [2011 03 11];   %% for JGR paper
yymmdd0 = [2024 11 13];   %% WHYMSIE

iaGlist  = JOB;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iPertTCC = +1;  %% use tcc model 1 (best so far)
iPertTCC = -1;  %% use default tcc in ECMWF  <<<<<<<<<<<<<<<<< DEFAULT >>>>>>>>>>>>>>

iSlabCld_CumSumStrowORGeorge = +1; %% strow,  cumsum 9999, cloud at PEAK of wgt fcn <<<< DEFAULT >>>>>>>
iSlabCld_CumSumStrowORGeorge = -1; %% aumann, cumsum -1,   cloud at mean of cld profile

if iPertTCC > 0
  error('nupe not coded')
else
  cloud_set_defaults_run_maker_uvw
end

