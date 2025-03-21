%% run with
%% sbatch --array=N1-N2 --output='testslurm' sergio_matlab_jobB.sbatch
%% sbatch --array=1-48 sergio_matlab_jobB.sbatch 
%% N1 = 1, N2 = number of files to be processed

%% note sergio_matlab_jobB.sbatch reminds you that
%% echo cmd line arg = 6, making one whole day interp ERA, SYMBOLIC LINK to clustbatch_make_eracloudrtp_sergio_sarta_filelist_interp_YYMMDD_loopGG.m
%% ie 
%% ls -lt clustbatch_eracloudrtp_sarta_filelist_interp_YYMMDD_loopGG.m
%%     clustbatch_eracloudrtp_sarta_filelist_interp_YYMMDD_loopGG.m -> clustbatch_make_eracloudrtp_sergio_sarta_filelist_interp_YYMMDD_loopGG.m

addpath /asl/matlib/h4tools
addpath /asl/matlib/rtptools
addpath /asl/matlib/aslutil
addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/TIME
addpath /home/sergio/MATLABCODE/PLOTTER
addpath /home/sergio/MATLABCODE/matlib/clouds/sarta/

system_slurm_stats

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
if length(JOB) == 0
  JOB = 1;
end
%JOB = 20

warning('off', 'MATLAB:imagesci:hdfeos:removalWarningHDFSW');

%% specify text file which has YY MM DD GG lst that needs to be processed
set_filelist

thefilelist = load(filelist);
thefilelist = thefilelist(JOB,1:3);

%%thefilelist = [2019 04 25];
iaGlist = 001 : 240;

%%%%%%%%%%%%%%%%%%%%%%%%%
%% this is for /home/sergio/MATLABCODE/CRODGERS_FAST_CLOUD/Various/tonga_volcano_jan2022_jpss.txt : 
%%  if JOB = 1:6 then this will do HungaTonga days 14 + (1:6) = 15-20, those 6 missing granules
%%  iInterp = -1;
%%  thefilelist = [2022 01 JOB+14];
%%  iaGlist = [52 53 54 55 56 240];
%%  thefilelist = [2022 01 JOB+14];
%%  iaGlist = [52 240];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ggx = 1 : length(iaGlist)
  gg = iaGlist(ggx);
  
  yymmdd0 = thefilelist;
  yymmdddggstr = ['.' num2str(thefilelist(1),'%04d') '.' num2str(thefilelist(2),'%02d') '.' num2str(thefilelist(3),'%02d') '.'];
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  iPertTCC = +1;  %% use tcc model 1 (best so far)
  iPertTCC = -1;  %% use default tcc in ECMWF  <<<<<<<<<<<<<<<<< DEFAULT >>>>>>>>>>>>>>
  
  iSlabCld_CumSumStrowORGeorge = +1; %% strow,  cumsum 9999, cloud at PEAK of wgt fcn <<<< DEFAULT >>>>>>>
  iSlabCld_CumSumStrowORGeorge = -1; %% aumann, cumsum -1,   cloud at mean of cld profile
  
  yy = yymmdd0(1); mm = yymmdd0(2); dd = yymmdd0(3); gg = iaGlist(ggx);
  
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
      [hd0, ha0, pd0, pa0, tstr] = cris_l1c_to_rtp_sergio(yy,mm,dd,gg,'era',  iSNPPorJ1orJ2,iInterp);
    elseif iERAorECM == -1  
      [hd0, ha0, pd0, pa0, tstr] = cris_l1c_to_rtp_sergio(yy,mm,dd,gg,'ecmwf',iSNPPorJ1orJ2,iInterp);
    end
  
    rtpwrite([dout '/' fout],hd0, ha0, pd0, pa0);
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

end
