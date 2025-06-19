%% run with
%% sbatch --array=N1-N2 --output='testslurm' sergio_matlab_jobB.sbatch
%% N1 = 1, N2 = number of files to be processed

%% specify text file which has YY MM DD GG lst that needs to be processed 
%% see eg https://www.ssec.wisc.edu/datacenter/NOAA20/GLOBAL2021_01_25_025.gif for NOAA20
%%        https://www.ssec.wisc.edu/datacenter/polar_orbit_tracks/
%%        https://www.ssec.wisc.edu/datacenter/polar_orbit_tracks/data/JPSS-1/
%%        https://www.ssec.wisc.edu/datacenter/polar_orbit_tracks/data/NPP/    ------------------->>>>>>

addpath /asl/matlib/h4tools
addpath /asl/matlib/rtptools
addpath /asl/matlib/aslutil
addpath /home/sergio/MATLABCODE/TIME
addpath /home/sergio/MATLABCODE/PLOTTER
addpath /home/sergio/MATLABCODE/matlib/clouds/sarta/
addpath /home/sergio/git/rtp_prod2/util

addpath /home/sergio/MATLABCODE

%%  check_all_jobs_done('/asl/s1/sergio/rtp/j1_ccast_hires/allfov/2019/04/25//cloudy_airs_l1c_ecm_sarta_baum_ice.2019.04.25.',240,'.rtp');

system_slurm_stats

%set_filelist

if ~exist('iInterp')
  iInterp = +1;
  iInterp = -1;
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
if length(JOB) == 0
  JOB = 120;
  JOB = 214;
  JOB = 212; JOB = 099; %% LA fires from D. Tobin
end

%JOB = 53

yymmdd0  = [2022 01 13]; ddLoop = [];  %% ECMWF says ATMS shows gravity waves from Tonga
yymmdd0  = [2022 01 15]; ddLoop = [14 : 22];  %% ECMWF says ATMS shows gravity waves from Tonga
yymmdd0  = [2019 04 25]; ddLoop = [];  %% HALO day that Eric processed, see ~/MATLABCODE/CRODGERS_FAST_CLOUD/HALO_BdryLayer/Proposal2024/driver_compare_AI.m
yymmdd0  = [2019 04 26]; ddLoop = [];  %% rather surprisingly, he did this day???
yymmdd0  = [2025 01 08]; ddLoop = [];  %% FIres over LA, from Dave Tobin
yymmdd0  = [2024 11 13]; ddLoop = [];  %% FIres over LA, from Dave Tobin

if length(ddLoop) == 0
  ddLoop = yymmdd0(3);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iiddloop = 1 : length(ddLoop)
  yymmdd0 = [yymmdd0(1:2) ddLoop(iiddloop)]

  thefilelist = yymmdd0;
  iaGlist  = 001 : 240;
  iaGlist = iaGlist(JOB);
  gg = iaGlist;
  
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

    addpath /home/sergio/MATLABCODE/PLOTTER

    i900 = find(hd0.vchan >= 900,1);
    tobs = rad2bt(900,pd0.robs1(i900,:));
    tclr = rad2bt(900,pd0.sarta_rclearcalc(i900,:));
    tcld = rad2bt(900,pd0.rcalc(i900,:));
    figure(1); clf; scatter_coast(pd0.rlon,pd0.rlat,25,tobs); title('BT 900 obs FSR'); cx1 = caxis; colormap jet
    figure(2); clf; scatter_coast(pd0.rlon,pd0.rlat,25,tclr); title('BT 900 clr');     cx2 = caxis; colormap jet
    figure(3); clf; scatter_coast(pd0.rlon,pd0.rlat,25,tcld); title('BT 900 cld');     cx3 = caxis; colormap jet

%{
%% CO
    i2160 = find(hd0.vchan >= 2160.00,1);
    i2162 = find(hd0.vchan >= 2161.75,1);
    tobs = rad2bt(2160,pd0.robs1(i2160,:))-rad2bt(2162,pd0.robs1(i2162,:));
    tclr = rad2bt(2160,pd0.sarta_rclearcalc(i2160,:)) - rad2bt(2162,pd0.sarta_rclearcalc(i2162,:));
    tcld = rad2bt(2160,pd0.rcalc(i2160,:)) - rad2bt(2162,pd0.rcalc(i2162,:));

    figure(1); clf; scatter_coast(pd0.rlon,pd0.rlat,25,tobs); title('BT 2161 obs FSR'); cx1 = caxis; colormap jet
    figure(2); clf; scatter_coast(pd0.rlon,pd0.rlat,25,tclr); title('BT 2161 clr');     cx2 = caxis; colormap jet
    figure(3); clf; scatter_coast(pd0.rlon,pd0.rlat,25,tcld); title('BT 2161 cld');     cx3 = caxis; colormap jet
%}
  
    cx(1) = min([cx1(1)  cx2(1) cx3(1)]);
    cx(2) = max([cx1(2)  cx2(2) cx3(2)]);
    figure(1); caxis(cx); 
    figure(2); caxis(cx); 
    figure(3); caxis(cx); 
    fprintf(1,'DONE : %s written out  \n',[dout '/' fout])  
  else
    fprintf(1,'%s already exists \n',[dout '/' fout])
  end

end
