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
addpath0

%%  check_all_jobs_done('/asl/s1/sergio/rtp/j1_ccast_hires/allfov/2019/04/25//cloudy_airs_l1c_ecm_sarta_baum_ice.2019.04.25.',240,'.rtp');

system_slurm_stats

if ~exist('iUVW')
  iUVW = -1;   %% do not add in windspeeds u,v,w
  iUVW = +1;   %% do     add in windspeeds u,v,w  
end

disp('assumes basic uvw file already made by clustbatch_make_eraORecm_cloudrtp_sergio_sarta_YYMMDD_loopGG.m')

JOB = str2num(getenv('SLURM_ARRAY_TASK_ID'));
if length(JOB) == 0
  JOB = 213;  %% S.Africa to Antartica
  JOB = 200;  %% california  
  JOB = 49;   %% daytime Australia, getting veddy veddy high PBLH
  JOB = 48;   %% daytime Australia, getting veddy veddy high PBLH  
end

gran = JOB;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% see eg /umbc/rs/pi_sergio/WorkDirDec2025/matlabcode/PBL_Retrievals/HALO_BdryLayer/PBL_Hgt_from_poemNew/loop_driver_compute_PBLH_poemNew2.m
%% see eg /umbc/rs/pi_sergio/WorkDirDec2025/matlabcode/PBL_Retrievals/HALO_BdryLayer/PBL_Hgt_from_poemNew/cluster_driver_compute_PBLH_poemNew.m

%% these had wrong potential enegy equation (wrong speed term in denom) ugh!!!!
RiSTR = 'Ri0.27'; %% orig
RiSTR = 'Ri0.25'; %% era5
RiSTR = 'Ri0.20'; %% try

%% now fixed potential energy denom
RiSTR = 'Ri0.25Fixed'; %% era5
RiSTR = 'Ri0.27Fixed'; %% orig

%%%%%%%%%%%%%%%%%%%%%%%%%

%% using static energy
RiSTR = 'Ri0.25StaticE'; %% era5, using static energy
RiSTR = 'Ri0.20StaticE'; %% era5, using static energy

RiSTR = 'Ri0.25StaticE_NEW'; %% era5, using static energy; forced a restriction to 4 km, else PBLH = salt (and check lapse rate/stability at PBLH)
RiSTR = 'Ri0.25StaticE_NEW2'; %% era5, using static energy; forced a restriction to 4 km, else PBLH = salt (and check lapse rate/stability at PBLH); also using the claude code

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rtpIN         = ['/home/sergio/nogit/sergio_temp_rtp_files/j1_ccast_hires/allfov/2024/11/13/interp_analysis_cloudy_airs_l1c_ecm_sarta_baum_ice.2024.11.13.' num2str(gran,'%03d') '.rtp'];
wspeedfileIN  = ['/home/sergio/nogit/sergio_temp_rtp_files/j1_ccast_hires/allfov/2024/11/13/' RiSTR '/interp_analysis_uvw_cloudy_airs_l1c_ecm_sarta_baum_ice.2024.11.13.' num2str(gran,'%03d') '.mat'];

wspeedfileOUT = ['/home/sergio/nogit/sergio_temp_rtp_files/j1_ccast_hires/allfov/2024/11/13/XTRY/interp_analysis_uvw_cloudy_airs_l1c_ecm_sarta_baum_ice.2024.11.13.' num2str(gran,'%03d') '.mat'];
wspeedfileOUT = ['/home/sergio/nogit/sergio_temp_rtp_files/j1_ccast_hires/allfov/2024/11/13/YTRY/interp_analysis_uvw_cloudy_airs_l1c_ecm_sarta_baum_ice.2024.11.13.' num2str(gran,'%03d') '.mat'];
wspeedfileOUT = ['/home/sergio/nogit/sergio_temp_rtp_files/j1_ccast_hires/allfov/2024/11/13/ZTRY/interp_analysis_uvw_cloudy_airs_l1c_ecm_sarta_baum_ice.2024.11.13.' num2str(gran,'%03d') '.mat'];
wspeedfileOUT = ['/home/sergio/nogit/sergio_temp_rtp_files/j1_ccast_hires/allfov/2024/11/13/WTRY/interp_analysis_uvw_cloudy_airs_l1c_ecm_sarta_baum_ice.2024.11.13.' num2str(gran,'%03d') '.mat']; %% same as Ztry, uses skt and u10,v10

%% wspeedfileOUT = ['/home/sergio/nogit/sergio_temp_rtp_files/j1_ccast_hires/allfov/2024/11/13/MAXTRY/interp_analysis_uvw_cloudy_airs_l1c_ecm_sarta_baum_ice.2024.11.13.' num2str(gran,'%03d') '.mat']; %% same as Ztry, uses skt and u10,v10


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

xhd00   = xhd0;
xpdmat0 = xpdmat;
clear xhd0 xpdmat

[xhd0,xpdmat] = get_richardson_number_levels(h,ha,p,pa);

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
