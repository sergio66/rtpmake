%% run with
%% sbatch --array=N1-N2 sergio_matlab_jobB.sbatch  2 for no iInterp, no uvw
%% sbatch --array=N1-N2 sergio_matlab_jobB.sbatch -2 for no iInterp, yes uvw

%% sbatch --array=N1-N2 sergio_matlab_chip.sbatch  4 for iInterp, no uvw
%% sbatch --array=N1-N2 sergio_matlab_chip.sbatch -4 for iInterp, yes uvw

%% N1 = 1, N2 = number of files to be processed

%% specify text file which has YY MM DD GG lst that needs to be processed 
%% see eg https://www.ssec.wisc.edu/datacenter/NOAA20/GLOBAL2021_01_25_025.gif for NOAA20
%%        https://www.ssec.wisc.edu/datacenter/polar_orbit_tracks/
%%        https://www.ssec.wisc.edu/datacenter/polar_orbit_tracks/data/JPSS-1/
%%        https://www.ssec.wisc.edu/datacenter/polar_orbit_tracks/data/NPP/    ------------------->>>>>>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath0

%% airs_l1c ---> cris_fsr
%%  check_all_jobs_done('/asl/s1/sergio/rtp/j1_ccast_hires/allfov/2019/04/25//cloudy_cris_fsr_ecm_sarta_baum_ice.2019.04.25.',240,'.rtp');

system_slurm_stats

%set_filelist

if ~exist('iInterp')
  iInterp = -1;
  iInterp = +1;
end

if ~exist('iERAorECM')
  iERAorECM = +1; %% till June 2019
  iERAorECM = -1; %% after June 2019

  iERAorECM = +1; %% this is ERA5
end

if ~exist('iSNPPorJ1orJ2')
  iSNPPorJ1orJ2 = +0; %% SuomiNPP
  iSNPPorJ1orJ2 = +1; %% J1
end

if ~exist('iUVW')
  iUVW = -1;   %% do not add in windspeeds u,v,w
  iUVW = +1;   %% do     add in windspeeds u,v,w  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

JOB = str2num(getenv('SLURM_ARRAY_TASK_ID'));
if length(JOB) == 0
  JOB = 120;
  JOB = 214;
  JOB = 099; %% LA fires from D. Tobin
  JOB = 180;
  JOB = 21;
  
  JOB = 49;   %% daytime Australia, getting veddy veddy high PBLH
  JOB = 118;  %% africa
  JOB = 48;   %% daytime Australia, getting veddy veddy high PBLH
  JOB = 213;  %% southern ocean
  JOB = 203;  %% pacific ocean near Mexico and CA
  
end

warning('off', 'MATLAB:imagesci:hdfeos:removalWarningHDFSW');

%JOB = 53

yymmdd0  = [2022 01 13]; ddLoop = [];  %% ECMWF says ATMS shows gravity waves from Tonga
yymmdd0  = [2022 01 15]; ddLoop = [14 : 22];  %% ECMWF says ATMS shows gravity waves from Tonga
yymmdd0  = [2019 04 25]; ddLoop = [];  %% HALO day that Eric processed, see ~/MATLABCODE/CRODGERS_FAST_CLOUD/HALO_BdryLayer/Proposal2024/driver_compare_AI.m
yymmdd0  = [2019 04 26]; ddLoop = [];  %% rather surprisingly, he did this day???
yymmdd0  = [2025 01 08]; ddLoop = [];  %% FIres over LA, from Dave Tobin
yymmdd0  = [2024 11 13]; ddLoop = [];  %% WHYMSIE

if length(ddLoop) == 0
  ddLoop = yymmdd0(3);
end

iERAorECM = +1;
if iERAorECM == -1
  cfg.model = 'ecmwf';
elseif iERAorECM == +1
  cfg.model = 'era5';
else
  error('gsk;jksjlksjslkjhs')
end

if (strcmp(cfg.model,'era') | strcmp(cfg.model,'merra2')) & iUVW > 0
  cfg.model
  iUVW
  error('hmm have not coded this up')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iiddloop = 1 : length(ddLoop)
  yymmdd0 = [yymmdd0(1:2) ddLoop(iiddloop)];

  thefilelist = yymmdd0;
  iaGlist  = 001 : 240;
  iaGlist = iaGlist(JOB);
  gg = iaGlist;
  
  yymmdddggstr = ['.' num2str(thefilelist(1),'%04d') '.' num2str(thefilelist(2),'%02d') '.' num2str(thefilelist(3),'%02d') '.'];

  fprintf(1,'loop %3i of %3i  [yymmdd0 gg] = %4i/%2i/%2i %3i \n',iiddloop,length(ddLoop),[yymmdd0 gg])
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  iPertTCC = +1;  %% use tcc model 1 (best so far)
  iPertTCC = -1;  %% use default tcc in ECMWF  <<<<<<<<<<<<<<<<< DEFAULT >>>>>>>>>>>>>>
  
  iSlabCld_CumSumStrowORGeorge = +1; %% strow,  cumsum 9999, cloud at PEAK of wgt fcn <<<< DEFAULT >>>>>>>
  iSlabCld_CumSumStrowORGeorge = -1; %% aumann, cumsum -1,   cloud at mean of cld profile
  
  yy = yymmdd0(1); mm = yymmdd0(2); dd = yymmdd0(3); gg = iaGlist;

  dout_set
  if iUVW < 0
    fout_set
  else
    fout_set_uvw
  end
  
  if ~exist(dout)
    mker = ['!mkdir -p ' dout];
    eval(mker)
  end

  if iUVW < 0
    ee = dir([dout '/' fout]);
  else
    ee = dir([dout '/' fmatout]);
  end
  if length(ee) == 0
    if iERAorECM == 1
      [hd0, ha0, pd0, pa0, tstr] = cris_l1c_to_rtp_sergio_ecm_or_era5(yy,mm,dd,gg,'era5',iSNPPorJ1orJ2,iInterp,iUVW);
    elseif iERAorECM == -1  
      [hd0, ha0, pd0, pa0, tstr] = cris_l1c_to_rtp_sergio_ecm_or_era5(yy,mm,dd,gg,'ecmwf',iSNPPorJ1orJ2,iInterp,iUVW);
    end

    if iUVW < 0
      rtpwrite([dout '/' fout],hd0, ha0, pd0, pa0);
    else
      rtpwrite([dout '/' fout],hd0, ha0, pd0, pa0);    
      %[xhd0,xpdmat] = get_richardson_number_levels(hd0,ha0,pd0,pa0);
      %[xhd0,xpdmat] = get_richardson_number_levels_v2(hd0,ha0,pd0,pa0);
      [xhd0,xpdmat] = get_richardson_number_levels_v3(hd0,ha0,pd0,pa0);      
      set_Ri_critical

  %%%%%%%%%%%%%%%%%%%%%%%%%
  boo = find(hd0.vchan >= 1231,1); %% 754
  figure(5); scatter_coast(pd0.rlon,pd0.rlat,20,rad2bt(1231,pd0.robs1(754,:)));                          title('BT1231 obs')
  figure(6); scatter_coast(pd0.rlon,pd0.rlat,20,rad2bt(1231,pd0.rcalc(754,:)));                          title('BT1231 calc')
										                         
  figure(1); scatter_coast(pd0.rlon,pd0.rlat,20,pd0.stemp);                                              title('ERA5 SKT [K]')
  figure(2); scatter_coast(pd0.rlon,pd0.rlat,20,pd0.pblh_nwp);                                           title('ERA5 OFFICIAL PBLH [km]'); caxis([0 4])
  %%% p.salti is in meters but xpdmat.salti o is in km				                         
  figure(3); scatter_coast(pd0.rlon,pd0.rlat,20,xpdmat.zPBLH_Ri-xpdmat.salti);                           title('SERGIO PBLH [km] '); caxis([0 4])
  figure(4); scatter_coast(pd0.rlon,pd0.rlat,20,xpdmat.PBLH_Ri_flag);                                    title('SERGIO PBLH flag');
  figure(9); scatter_coast(pd0.rlon,pd0.rlat,20,xpdmat.zPBLH_Ri_from_compute_pblh_Ri/1000-xpdmat.salti); title('SERGIO from compute\_pblh\_Ri')
  figure(9); caxis([0 4]);
  keyboard_nowindow
  %% saverx = ['save cris_2024_11_18_g' num2str(JOB,'%03d') '.mat hd0 pd0 xhd0 xpdmat']; eval(saverx)
  %%%%%%%%%%%%%%%%%%%%%%%%%

      saver = ['save ' dout '/' fmatout  ' xhd0 xpdmat RiCritical iVers_Ri'];
      eval(saver)      
      %% plot_richardson_PBLH_levels

      %%
      %% see /umbc/rs/pi_sergio/WorkDirDec2025/matlabcode/PBL_Retrievals/HALO_BdryLayer/PBL_Hgt_from_poemNew/cluster_driver_compute_PBLH_poemNew.m
      %%   for use of [yhd0,ypdmat] = get_richardson_number_layers(hoemNew,poemNew,xhd0,xpdmat,iPlot);
      %%   saved into eg fnameOUT = ['/home/sergio/nogit/sergio_temp_rtp_files/j1_ccast_hires/allfov/2024/11/13/retr_cloudy_cris_fsr_ecm_sarta_baum_ice.2024.11.13.' num2str(gran,'%03d') '_layers_PBLH_Ri.mat'];

      figure(4); clf; scatter_coast(xpdmat.rlon,xpdmat.rlat,50,xpdmat.zPBLH_Ri/1000); title('Ri PBLH [km]')
    end
    
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
    if iUVW < 0    
      fprintf(1,'DONE : %s written out  \n',[dout '/' fout])
    else
      fprintf(1,'DONE : %s written out  \n',[dout '/' fmatout])
    end
    
  else
    if iUVW < 0
      fprintf(1,'%s already exists \n',[dout '/' fout])
    elseif iUVW > 0
      fprintf(1,'%s already exists \n',[dout '/' fmatout])
    end
    
  end

end
