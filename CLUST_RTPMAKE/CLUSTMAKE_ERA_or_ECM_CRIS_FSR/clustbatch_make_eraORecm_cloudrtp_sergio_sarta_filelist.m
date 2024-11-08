%% run with
%%   sbatch --array=N1-N2 --output='testslurm' sergio_matlab_jobB.sbatch
%%   sbatch -p cpu2021 --array=1-Y   sergio_matlab_jobB.sbatch 1
%% N1 = 1, N2 = number of files to be processed

%% specify text file which has YY MM DD GG lst that needs to be processed 
%% see eg https://www.ssec.wisc.edu/datacenter/NOAA20/GLOBAL2021_01_25_025.gif for NOAA20
%%        https://www.ssec.wisc.edu/datacenter/polar_orbit_tracks/
%%        https://www.ssec.wisc.edu/datacenter/polar_orbit_tracks/data/JPSS-1/
%%        https://www.ssec.wisc.edu/datacenter/polar_orbit_tracks/data/NPP/    ------------------->>>>>>

addpath /asl/matlib/h4tools
addpath /asl/matlib/rtptools
addpath /asl/matlib/aslutil
addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/TIME
addpath /home/sergio/MATLABCODE/PLOTTER
addpath /home/sergio/MATLABCODE/matlib/clouds/sarta/

system_slurm_stats

set_filelist

if ~exist('iInterp')
  iInterp = +1;
end

if ~exist('iERAorECM')
  iERAorECM = +1; %% till June 2019
  iERAorECM = -1; %% after June 2019
end

if ~exist('iSNPPorJ1orJ2')
  iSNPPorJ1orJ2 = +0; %% SuomiNPP
  iSNPPorJ1orJ2 = +1; %% J1
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

JOB = str2num(getenv('SLURM_ARRAY_TASK_ID'));
%JOB = 14  %% NPP H. Dorian 2019/09/03 at 6.45 am
%JOB = 15  %% J-1 H. Dorian 2019/09/03 at 7.30 am
%JOB = 4   %% NPP HALO lidar 2019/04/25 at 22.10 pm --> JPSS 2019/04/25 at 21.20 pm, kinda misses HALO
%JOB = 3   %% NPP HALO lidar 2019/04/25 at 22.10 pm --> JPSS 2019/04/25 at 21.25 pm, gets most of HALO /asl/s1/sergio/rtp/j1_ccast_hires/allfov/2019/04/25//cloudy_airs_l1c_ecm_sarta_baum_ice.2019.04.25.214.rtp

if length(JOB) == 0
  JOB = 1;
  JOB = 6;
end

thefilelist = load(filelist);
thefilelist = thefilelist(JOB,:);

yymmdd0      = thefilelist(1:3); %% YY MM DD
iaGlist      = thefilelist(4);   %% granule
yymmdddggstr = ['.' num2str(thefilelist(1),'%04d') '.' num2str(thefilelist(2),'%02d') '.' num2str(thefilelist(3),'%02d') '.'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iPertTCC = +1;  %% use tcc model 1 (best so far)
iPertTCC = -1;  %% use default tcc in ECMWF  <<<<<<<<<<<<<<<<< DEFAULT >>>>>>>>>>>>>>

iSlabCld_CumSumStrowORGeorge = +1; %% strow,  cumsum 9999, cloud at PEAK of wgt fcn <<<< DEFAULT >>>>>>>
iSlabCld_CumSumStrowORGeorge = -1; %% aumann, cumsum -1,   cloud at mean of cld profile

yy = yymmdd0(1); mm = yymmdd0(2); dd = yymmdd0(3); gg = iaGlist;

if iSNPPorJ1orJ2 == 0
  NONONOdout = ['/asl/rtp/cris/npp_ccast_hires/allfov/' num2str(yy,'%04d') '/' num2str(mm,'%02d') '/' num2str(dd,'%02d') '/'];
  dout = ['/asl/s1/sergio/rtp/npp_ccast_hires/allfov/' num2str(yy,'%04d') '/' num2str(mm,'%02d') '/' num2str(dd,'%02d') '/'];
elseif iSNPPorJ1orJ2 == 1
  NONONOdout = ['/asl/rtp/cris/j1_ccast_hires/allfov/' num2str(yy,'%04d') '/' num2str(mm,'%02d') '/' num2str(dd,'%02d') '/'];
  dout = ['/asl/s1/sergio/rtp/j1_ccast_hires/allfov/' num2str(yy,'%04d') '/' num2str(mm,'%02d') '/' num2str(dd,'%02d') '/'];
else
  error('unknow SNPP, J1 or ... ?')
end
if ~exist(dout)
  mker = ['!mkdir -p ' dout];
  eval(mker)
end

if iInterp <= 0
  if iERAorECM == 1
    fout = ['fsr_allfov_era_' num2str(gg,'%03d') '*.rtp'];
    fout = ['cloudy_airs_l1c_era_sarta_baum_ice' yymmdddggstr '*.rtp'];
    fout = ['cloudy_airs_l1c_era_sarta_baum_ice' yymmdddggstr num2str(gg,'%03d') '.rtp'];
  elseif iERAorECM == -1
    fout = ['fsr_allfov_ecm_' num2str(gg,'%03d') '*.rtp'];
    fout = ['cloudy_airs_l1c_ecm_sarta_baum_ice' yymmdddggstr '*.rtp'];
    fout = ['cloudy_airs_l1c_ecm_sarta_baum_ice' yymmdddggstr num2str(gg,'%03d') '.rtp'];
  end
else
  if iERAorECM == 1
    fout = ['interp_analysis_fsr_allfov_era_' num2str(gg,'%03d') '*.rtp'];
    fout = ['interp_analysis_cloudy_airs_l1c_era_sarta_baum_ice' yymmdddggstr '*.rtp'];
    fout = ['interp_analysis_cloudy_airs_l1c_era_sarta_baum_ice' yymmdddggstr num2str(gg,'%03d') '.rtp'];
  elseif iERAorECM == -1
    fout = ['interp_analysis_fsr_allfov_ecm_' num2str(gg,'%03d') '*.rtp'];
    fout = ['interp_analysis_cloudy_airs_l1c_ecm_sarta_baum_ice' yymmdddggstr '*.rtp'];
    fout = ['interp_analysis_cloudy_airs_l1c_ecm_sarta_baum_ice' yymmdddggstr num2str(gg,'%03d') '.rtp'];
  end
end
  
ee = dir([dout '/' fout]);
if length(ee) == 0
  if iERAorECM == 1
    [hd0, ha0, pd0, pa0, tstr] = cris_l1c_to_rtp_sergio(yy,mm,dd,gg,'era',iSNPPorJ1orJ2,iInterp);
  elseif iERAorECM == -1  
    [hd0, ha0, pd0, pa0, tstr] = cris_l1c_to_rtp_sergio(yy,mm,dd,gg,'ecmwf',iSNPPorJ1orJ2,iInterp);
  end

  rtpwrite([dout '/' fout],hd0, ha0, pd0, pa0);
  fprintf(1,'wrote out %s \n',[dout '/' fout]);

  i900 = find(hd0.vchan >= 900,1);
  tobs = rad2bt(900,pd0.robs1(i900,:));
  tclr = rad2bt(900,pd0.sarta_rclearcalc(i900,:));
  tcld = rad2bt(900,pd0.rcalc(i900,:));
  addpath /home/sergio/MATLABCODE/PLOTTER
  figure(1); clf; scatter_coast(pd0.rlon,pd0.rlat,25,tobs); title('BT 900 obs FSR'); cx1 = caxis; colormap jet
  figure(2); clf; scatter_coast(pd0.rlon,pd0.rlat,25,tclr); title('BT 900 clr');     cx2 = caxis; colormap jet
  figure(3); clf; scatter_coast(pd0.rlon,pd0.rlat,25,tcld); title('BT 900 cld');     cx3 = caxis; colormap jet

  cx(1) = min([cx1(1)  cx2(1) cx3(1)]);
  cx(2) = max([cx1(2)  cx2(2) cx3(2)]);
  figure(1); caxis(cx); 
  figure(2); caxis(cx); 
  figure(3); caxis(cx); 
  fprintf(1,'DONE : %s written out  \n',[dout '/' fout])  
else
  fprintf(1,'%s already exists \n',[dout '/' fout])
end
