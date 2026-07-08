%% run with
%% sbatch --array=N1-N2 sergio_matlab_jobB.sbatch  2 for no iInterp, no uvw
%% sbatch --array=N1-N2 sergio_matlab_jobB.sbatch -2 for no iInterp, yes uvw

%% sbatch --array=N1-N2 sergio_matlab_chip.sbatch  4 for iInterp, no uvw
%% sbatch --array=N1-N2 sergio_matlab_chip.sbatch -4 for iInterp, yes uvw
%% sbatch --array=N1-N2 sergio_matlab_chip.sbatch -44 for update iInterp, yes uvw

%% N1 = 1, N2 = number of files to be processed

%% specify text file which has YY MM DD GG lst that needs to be processed 
%% see eg https://www.ssec.wisc.edu/datacenter/NOAA20/GLOBAL2021_01_25_025.gif for NOAA20
%%        https://www.ssec.wisc.edu/datacenter/polar_orbit_tracks/
%%        https://www.ssec.wisc.edu/datacenter/polar_orbit_tracks/data/JPSS-1/
%%        https://www.ssec.wisc.edu/datacenter/polar_orbit_tracks/data/NPP/    ------------------->>>>>>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('assumes basic uvw file already made by clustbatch_make_eraORecm_cloudrtp_sergio_sarta_YYMMDD_loopGG.m')
disp('  and then takes in that file and reruns updates to get_richardson_number_levelsX.m')

addpath0

%%  check_all_jobs_done('/asl/s1/sergio/rtp/j1_ccast_hires/allfov/2019/04/25//cloudy_cris_fsr_ecm_sarta_baum_ice.2019.04.25.',240,'.rtp');

system_slurm_stats

if ~exist('iUVW')
  iUVW = -1;   %% do not add in windspeeds u,v,w
  iUVW = +1;   %% do     add in windspeeds u,v,w  
end

JOB = str2num(getenv('SLURM_ARRAY_TASK_ID'));
if length(JOB) == 0

  %%% these are JPSS CRIS 2024/11/13
  JOB = 49;   %% daytime Australia, getting veddy veddy high PBLH
  JOB = 48;   %% daytime Australia, getting veddy veddy high PBLH
  JOB = 213;  %% S.Africa to Antartica, ocean
  JOB = 200;  %% california, ocean
  JOB = 210;  %% indian ocean near india
  JOB = 203;  %% pacific ocean near Mexico and CA
  JOB = 235;  %% pacific oean near samoa? nah???!! middle of nowhere
  JOB = 118;  %% africa
end

gran = JOB;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% see eg /umbc/rs/pi_sergio/WorkDirDec2025/matlabcode/PBL_Retrievals/HALO_BdryLayer/PBL_Hgt_from_poemNew/loop_driver_compute_PBLH_poemNew2.m
%% see eg /umbc/rs/pi_sergio/WorkDirDec2025/matlabcode/PBL_Retrievals/HALO_BdryLayer/PBL_Hgt_from_poemNew/cluster_driver_compute_PBLH_poemNew.m
%% see eg /umbc/rs/pi_sergio/WorkDirDec2025/matlabcode/PBL_Retrievals/HALO_BdryLayer/PBL_Hgt_from_poemNew/set_RiSTR.m

%% these had wrong potential enegy equation (wrong speed term in denom) ugh!!!!
RiSTR_INFILE = 'Ri0.27'; %% orig
RiSTR_INFILE = 'Ri0.25'; %% era5
RiSTR_INFILE = 'Ri0.20'; %% try

%% now fixed potential energy denom
RiSTR_INFILE = 'Ri0.25Fixed'; %% era5
RiSTR_INFILE = 'Ri0.27Fixed'; %% orig

%% using static energy
RiSTR_INFILE = 'Ri0.25StaticE'; %% era5, using static energy
RiSTR_INFILE = 'Ri0.20StaticE'; %% era5, using static energy

RiSTR_INFILE = 'Ri0.25StaticE_NEW';  %% era5, using static energy; forced a restriction to 4 km, else PBLH = salt (and check lapse rate/stability at PBLH)
RiSTR_INFILE = 'Ri0.25StaticE_NEW2'; %% era5, using static energy; forced a restriction to 4 km, else PBLH = salt (and check lapse rate/stability at PBLH); also using the claude code
                                     %% made by get_richardson_number_levels? get_richardson_number_levels_v2?
				     
RiSTR_INFILE = 'T2TRY0'; %% ecmwf, using static energy; forced a restriction to 4 km, else PBLH = salt (and check lapse rate/stability at PBLH); also using the claude code ESP use t2m 
                         %% instead of SKT in numerator of Ri(z) made by get_richardson_number_levels_v3.m. Pretty good. Also has pd0.pblh_nwp created from the daily ERA5 pd0 = quick_get_ERA5_pblh(pd0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RiSTR_OUTFILE = 'XTRY';
RiSTR_OUTFILE = 'YTRY';
RiSTR_OUTFILE = 'ZTRY';
RiSTR_OUTFILE = 'WTRY';   %% same as Ztry,      uses skt and u10,v10
RiSTR_OUTFILE = 'T2TRY';  %% same as Wtry,Ztry, uses T2m and u10,v10  no need to have input RiSTR_INFILE = 'T2TRY'; and output RiSTR_OUTFILE = 'T2TRY';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rtpIN         = ['/home/sergio/nogit/sergio_temp_rtp_files/j1_ccast_hires/allfov/2024/11/13/interp_analysis_cloudy_cris_fsr_ecm_sarta_baum_ice.2024.11.13.' num2str(gran,'%03d') '.rtp'];
rtpIN         = ['/home/sergio/nogit/sergio_temp_rtp_files/j1_ccast_hires/allfov/2024/11/13/interp_analysis_cloudy_jpss_fsr_ecm_sarta_baum_ice.2024.11.13.' num2str(gran,'%03d') '.rtp'];

wspeedfileIN  = ['/home/sergio/nogit/sergio_temp_rtp_files/j1_ccast_hires/allfov/2024/11/13/' RiSTR_INFILE '/interp_analysis_uvw_cloudy_cris_fsr_ecm_sarta_baum_ice.2024.11.13.' num2str(gran,'%03d') '.mat'];

wspeedfileOUT = ['/home/sergio/nogit/sergio_temp_rtp_files/j1_ccast_hires/allfov/2024/11/13/' RiSTR_OUTFILE '/interp_analysis_uvw_cloudy_cris_fsr_ecm_sarta_baum_ice.2024.11.13.' num2str(gran,'%03d') '.mat']; 

if ~exist(wspeedfileIN)
  fprintf(1,'%s DNE \n',wspeedfileIN)
end
if ~exist(rtpIN)
  fprintf(1,'%s DNE \n',rtpIN)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[h,ha,p,pa] = rtpread(rtpIN);
loader = ['load ' wspeedfileIN];
eval(loader);

p.v10 = xpdmat.v10;
p.u10 = xpdmat.u10;
p.u   = xpdmat.u;
p.v   = xpdmat.v;
p.w   = xpdmat.w;
p.t2m = xpdmat.t2m;

if isfield(xpdmat,'pblh_nwp')
  p.pblh_nwp = xpdmat.pblh_nwp;
end

xhd00   = xhd0;
xpdmat0 = xpdmat;
clear xhd0 xpdmat

%[xhd0,xpdmat] = get_richardson_number_levels(h,ha,p,pa);
%[xhd0,xpdmat] = get_richardson_number_levels_v2(h,ha,p,pa);
[xhd0,xpdmat] = get_richardson_number_levels_v3(h,ha,p,pa);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

comment = 'see /umbc/rs/pi_sergio/WorkDirDec2025/matlabcode/PBL_Retrievals/HALO_BdryLayer/PBL_Hgt_from_poemNew/cluster_driver_compute_PBLH_poemNew.m';
set_Ri_critical
saver = ['save ' wspeedfileOUT ' xhd0 xpdmat RiCritical iVers_Ri'];
if ~exist(wspeedfileOUT)
  eval(saver)
  fprintf(1,'saved %s \n',wspeedfileOUT)
else
  fprintf(1,'%s already exists, not saving \n',wspeedfileOUT)
end
