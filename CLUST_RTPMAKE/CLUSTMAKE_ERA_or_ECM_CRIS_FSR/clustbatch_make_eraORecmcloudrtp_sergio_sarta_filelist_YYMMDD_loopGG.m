%% run with
%% sbatch --array=N1-N2 --output='testslurm' sergio_matlab_jobB.sbatch 
%% sbatch --array=1-48 sergio_matlab_jobB.sbatch 
%% N1 = 1, N2 = number of files to be processed

%% note sergio_matlab_jobB.sbatch reminds you that
%% echo cmd line arg = 6, making one whole day interp ERA, SYMBOLIC LINK to clustbatch_make_eracloudrtp_sergio_sarta_filelist_interp_YYMMDD_loopGG.m
%% ie 
%% ls -lt clustbatch_eracloudrtp_sarta_filelist_interp_YYMMDD_loopGG.m
%%     clustbatch_eracloudrtp_sarta_filelist_interp_YYMMDD_loopGG.m -> clustbatch_make_eracloudrtp_sergio_sarta_filelist_interp_YYMMDD_loopGG.m

%% [sergio@chip-login1 CLUSTMAKE_ERA_or_ECM_CRIS_FSR]$ ls -lt /home/sergio/git/matlabcode/L2Readers/GetL1_JPSS_CrIS
%% lrwxrwxrwx 1 sergio pi_sergio  67 Dec 19 05:44 Readme_Get_JPSS1_CriS_L1 -> ../../GET_NWP_ERA5_MERRA2_ECMWF_data_NOTES/Readme_Get_JPSS1_CriS_L1
%% so this is /home/sergio/git/matlabcode/GET_NWP_ERA5_MERRA2_ECMWF_data_NOTES/Readme_Get_JPSS1_CriS_L1
%% which says get L1 data from    //sounder.gesdisc.eosdis.nasa.gov/data/JPSS1_Sounder_Level1/SNDRJ1CrISL1B.2/

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath0

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

iaGlist = 234;   %% testing
iaGlist = 236;   %% testing
iaGlist = 182;   %% testing 
iaGlist = 001 : 240;
iaGlist = 180 : 240;
iaGlist = 220;   %% testing

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

  dout_set
  fout_set
    
  if ~exist(dout)
    mker = ['!mkdir -p ' dout];
    eval(mker)
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
