%% run with
%% sbatch --array=N1-N2 --output='testslurm' sergio_matlab_jobB.sbatch
%% N1 = 1, N2 = number of files to be processed

%% specify text file which has YY MM DD GG lst that needs to be processed 
filelist = '/home/sergio/MATLABCODE/MBL_CLOUDS/MAGIC/magicsonde_matchup_airsfile_list.txt';
filelist = '/home/sergio/MATLABCODE/MBL_CLOUDS/TWP_quicksonde/twp_stuff.txt';
filelist = '/home/sergio/MATLABCODE/MBL_CLOUDS/TWP_quicksonde/twp_stuff_NEW.txt';
filelist = '/home/sergio/MATLABCODE/MBL_CLOUDS/TWP_quicksonde/twp_stuff_FIXED.txt';
filelist = '/home/sergio/MATLABCODE/RTPMAKE/RESULTS_for_PAPERS/MFILES_AtmosphericRiver/atmospheric_river_2014_02_08.txt';
filelist = '/home/sergio/MATLABCODE/RTPMAKE/CLUST_RTPMAKE/CLUSTMAKE_ECM/hurricane_2002_09_06_g044.txt';
filelist = '/home/sergio/MATLABCODE/CRODGERS_FAST_CLOUD/SONDE_VALIDATION/Strow_LIN/lin_filelist.txt';
filelist = '/home/sergio/MATLABCODE/RTPMAKE/CLUST_RTPMAKE/CLUSTMAKE_ERA/sullivan_ozone_intrusin_aug_2014.txt';
filelist = '/home/sergio/MATLABCODE/RTPMAKE/CLUST_RTPMAKE/CLUSTMAKE_ERA/dcc_2011_03_11_g039.txt';

JOB = str2num(getenv('SLURM_ARRAY_TASK_ID'));
JOB = 1

thefilelist = load(filelist);
thefilelist = thefilelist(JOB,:);

yymmdd0  = thefilelist(1:3); %% YY MM DD
iaGlist  = thefilelist(4);   %% granule
     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% creates an rtp file for ONE granule
%% can be modified for more!

klayers = '/asl/packages/klayers/Bin/klayers_airs';
sarta   = '/asl/packages/sartaV108/Bin/sarta_apr08_m140_wcon_nte';

klayers = '/asl/packages/klayers/Bin/klayers_airs';
sarta   = '/asl/packages/sartaV108/Bin/sarta_apr08_m140_wcon_nte';

klayers = '/asl/packages/klayersV205/BinV201/klayers_airs';
sarta   = '/asl/packages/sartaV108_PGEv6/Bin/sarta_airs_PGEv6_postNov2003';

addpath /asl/matlab2012/airs/readers
addpath /asl/matlib/aslutil
addpath /asl/matlib/science
addpath /asl/matlib/rtptools
addpath /asl/matlib/h4tools/
addpath /asl/matlib/rtptools/
addpath /asl/matlib/gribtools/
addpath /home/sergio/MATLABCODE/matlib/clouds/sarta

%addpath /home/strow/cress/Work/Rtp
% addpath /home/strow/Matlab/Grib     WARNING /home/strow/Matlab/Grib/rtpadd_grib_data.m DIFFERENT than /asl/matlib/gribtools/rtpadd_era_data.m
% addpath /home/sergio/MATLABCODE/CRIS_HiRes             %% for sergio_fill_ecmwf
% addpath /home/strow/Git/rtp_prod2/grib                  %% for fill_ecm
addpath /asl/packages/rtp_prod2/grib

theinds = (1 : 2378)';

%%%%%%%%%%%%%%%%%%%%%%%%%
%%% set sarta exec

run_sarta.clear = +1;
run_sarta.cloud = +1;
run_sarta.cumsum = -1;    %% this is "closer" to MRO but since cliuds are at centroid, does not do too well with DCC
run_sarta.cumsum = 9999;  %% larrabee likes this, puts clouds high so does well for DCC

codeX = 0; %% use default with A. Baran params
codeX = 1; %% use new     with B. Baum, P. Yang params

code0 = '/asl/packages/sartaV108/BinV201/sarta_apr08_m140_iceaggr_waterdrop_desertdust_slabcloud_hg3_wcon_nte';
code1 = '/home/sergio/SARTA_CLOUDY/BinV201/sarta_apr08_m140x_iceGHMbaum_waterdrop_desertdust_slabcloud_hg3';

if codeX == 0
  icestr = '_sarta_baran_ice';
  run_sarta.sartacloud_code = code0;
elseif codeX == 1
  icestr = '_sarta_baum_ice';
  run_sarta.sartacloud_code = code1;
else
  error('codeX???')
end

%icestr = ['NEWLANDFRAC/cloudy_airs_l1b_era' icestr '.'];
%icestr = ['cloudy_airs_l1b_era' icestr '.'];
icestr = ['xcloudy_airs_l1b_era' icestr '.'];

%%%%%%%%%%%%%%%%%%%%%%%%%

for ixx = 1 : length(iaGlist)
  ix = iaGlist(ixx);
  clear p h hattr pattr prof yymmddgg

  yymmddgg = [yymmdd0 ix];

  ystr = num2str(yymmddgg(1));
  mstr = num2str(yymmddgg(2),'%02d');
  dstr = num2str(yymmddgg(3),'%02d');
  gstr = num2str(yymmddgg(4),'%03d');

  fdirOUT = ['/asl/data/rtprod_airs/' ystr '/' mstr '/' dstr '/'];
  fdirOUT = ['/asl/rtp/rtprod_airs/' ystr '/' mstr '/' dstr '/'];
  fdirOUT = ['/asl/s1/sergio/rtp/rtp_airibrad_v5/' ystr '/' mstr '/' dstr '/'];

  if ~exist(fdirOUT)
    mker = ['!/bin/mkdir -p ' fdirOUT];
    fprintf(1,'mker = %s \n',mker);
    eval(mker)
  end
  
  fnameOUT= [fdirOUT icestr ystr '.' mstr '.' dstr '.' gstr '.rtp'];

  eeP = exist(fnameOUT);

  if eeP == 0
    fprintf(1,' making %s \n',fnameOUT);
    toucher = ['!touch ' fnameOUT];
    eval(toucher)
  
    year  = yymmddgg(1);
    month = yymmddgg(2);
    day   = yymmddgg(3);
    gran  = yymmddgg(4);
    if mod(year,4) == 0
      mos = [31 29 31 30 31 30 31 31 30 31 30 31];  %% leap year
    else
      mos = [31 28 31 30 31 30 31 31 30 31 30 31];  %% normal year
      end
    days_so_far = 0;
    if month > 1
      days_so_far = sum(mos(1:month-1));
    end
    days_so_far = days_so_far + day;
    filename = ['/strowdataN/data/airs/Aqua_AIRS_Level1/AIRIBRAD.005/' ystr '/'];
    filename = ['/asl/data/airs/AIRIBRAD/' ystr '/'];
    filename = [filename num2str(days_so_far,'%03d') '/'];
    dir0 = filename;
    filename = [filename 'AIRS.' ystr '.' mstr '.' dstr '.' gstr];
    filename = [filename '.L1B.AIRS_Rad.v5*.hdf'];

    thedir = dir(filename);
    if length(thedir) == 1
      fname = [dir0 thedir.name];
    else
      fprintf(1,'%s \n',filename);
      disp('file does not exist');
      return
    end

    [a,b,c] = sdload(fname);
    p.rlat = [a.Latitude(:)'];
    p.rlon = [a.Longitude(:)'];
    p.rtime = [a.Time(:)'];

    %[meantime, f, prof] = readl1b_all(fname);  %% has the old rtime 1993
    [meantime, f, prof] = xreadl1b_all(fname);  %% has the new rtime 1958
    p = prof;

    p.pobs = zeros(size(p.solazi));
    p.upwell = ones(size(p.solazi));
    %p.irinst = AIRSinst*ones(1,nobs);
    %p.findex = grannum*ones(1,nobs);

    plot(p.rlon,p.rlat,'.')

    pa = {{'profiles','rtime','seconds since 1993'}};
    ha = {{'header','hdf file',filename}};

    h.pfields=5; % (1=prof + 4=IRobs);

    h.nchan = length(theinds);
    h.ichan = theinds;;
    h.vchan = f(h.ichan);;

    %%% this is NEW
    p.landfrac_fromL1B = p.landfrac;
    p.salti_fromL1B = p.salti;
    [salti, landfrac] = usgs_deg10_dem(p.rlat, p.rlon);
    p.landfrac = landfrac;
    p.salti    = salti;

    clrfields = {'SP','SKT','10U','10V','TCC','CI','T','Q','O3'};
    cldfields = {'SP','SKT','10U','10V','TCC','CI','T','Q','O3',...
                 'CC','CIWC','CLWC'};

    %[h,ha,p,pa] = rtpadd_era_data(h,ha,p,pa,cldfields); %%% add on era
    [p,h] = fill_era(p,h);
    
    p0 = p;

    %[h,ha,p,pa] = rtpadd_emis_DanZhou2(h,ha,p,pa);
    %p = Prof_add_emis(p,yymmddgg(1),yymmddgg(2),yymmddgg(3));  %% broken crap by whoever
    %p = rtpadd_emis_DanZhou(h,ha,p,pa);   %% lso totally broken crap
    %[h,ha,p,pa] = rtpadd_emis_wis(h,ha,p,pa);
    %addpath /asl/rtp_prod2/emis/
    %addpath /asl/rtp_prod2/util/
    addpath /asl/packages/rtp_prod2/emis/
    addpath /asl/packages/rtp_prod2/util/    
    p.rlon = wrapTo180(p.rlon);
    [p,pa] = rtp_add_emis(p,pa);
    
    %figure(1)
    %scatter_coast(p.rlon,p.rlat,10,p.nemis); 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% default stuff, what I usually do
    p0 = p;
    run_sarta.clear = -1;
    [p0] = driver_sarta_cloud_rtp(h,ha,p,pa,run_sarta);

    %%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% no sarta calcs needed at this point, just do the first cut at cloud slabs
    xrun_sarta = run_sarta;
    xrun_sarta.clear = -1;
    xrun_sarta.cloud = -1;
    xrun_sarta.klayers_code = klayers;
    [xp] = driver_sarta_cloud_rtp(h,ha,p,pa,xrun_sarta);

    tempx.sarta = run_sarta.sartacloud_code;
    tempx.klayers = klayers;

    %% now optimize the cloud slabs
    addpath /home/sergio/MATLABCODE/RTPMAKE/CLUST_RTPMAKE/SUBSET_BEST_CLOUD
    addpath /home/sergio/MATLABCODE/PLOTTER    
    [pjunk,g,FOV_ECM,fxx,alltOBS,alltECM,closestBT] = subset_best_cloud(h,ha,xp,pa,tempx)
    
    %% finally subset the random profiles and run klayers/sarta again
    p.cprtop = pjunk.cprtop;
    p.cprbot = pjunk.cprbot;
    p.cngwat = pjunk.cngwat;
    p.cpsize = pjunk.cpsize;
    p.cfrac = pjunk.cfrac;
    p.ctype = pjunk.ctype;
    p.cprtop2 = pjunk.cprtop2;
    p.cprbot2 = pjunk.cprbot2;
    p.cngwat2 = pjunk.cngwat2;
    p.cpsize2 = pjunk.cpsize2;
    p.cfrac2 = pjunk.cfrac2;
    p.cfrac12 = pjunk.cfrac12;    
    p.ctype2 = pjunk.ctype2;

    p = sanity_check(p);
    p = sanity_check(p);

    %% AT THIS POINT CAN CHOOSE THE RANDOM PROFILES and subset .... or if not, then do all
    [h,p] = subset_rtp_allcloudfields(h,p,[],[],random_selection_of_profiles);
    
    %% the final she-bang .... does clear and cloud calcs etc
    run_sarta.ForceNewSlabs = -1;   %% do not alter the slabs now!!!    
    run_sarta.clear = +1;           
    [p2] = driver_sarta_cloud_rtp(h,ha,p,pa,run_sarta);

tobs = rad2bt(1231,p2.robs1(1291,:));
tcal0 = rad2bt(1231,p0.rcalc(1291,:));
tcalf = rad2bt(1231,p2.rcalc(1291,:));
figure(1); scatter_coast(p0.rlon,p0.rlat,10,real(tobs)); title('BT1231 obs'); colormap jet
figure(2); scatter_coast(p0.rlon,p0.rlat,10,real(tcal0)); title('BT1231 cal0'); colormap jet
figure(3); scatter_coast(p0.rlon,p0.rlat,10,real(tcalf)); title('BT1231 calf'); colormap jet

tobs = rad2bt(h.vchan,p2.robs1);
tcal0 = rad2bt(h.vchan,p0.rcalc);
tcalf = rad2bt(h.vchan,p2.rcalc);
figure(4);
  plot(h.vchan(g),nanmean(tobs(g,:)'-tcal0(g,:)'),'b',h.vchan(g),nanstd(tobs(g,:)'-tcal0(g,:)'),'c',...
       h.vchan(g),nanmean(tobs(g,:)'-tcalf(g,:)'),'r',h.vchan(g),nanstd(tobs(g,:)'-tcalf(g,:)'),'m')
hl = legend('mean(obs-cal0)','std(obs-cal0)','mean(obs-calf)','std(obs-calf)','location','best'); set(hl,'fontsize',10);
grid;

    fnamex = fnameOUT;
    [h,ha,p2x,pa] = rtptrim_sartacloud(h,ha,p2,pa);
    rtpwrite(fnamex,h,ha,p2x,pa)
  else
    fprintf(1,' %s already exists \n',fnameOUT)
  end
end
