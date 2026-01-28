%{
https://www.mathworks.com/help/matlab/matlab_prog/suppress-warnings.html?requestedDomain=www.mathworks.com

Warning: HDFSD will be removed in a future release. Use MATLAB.IO.HDF4.SD instead.
> In hdfsd (line 261)
In sdload (line 116)
In cloud_set_defaults_run_maker (line 121)
In clustbatch_make_ecmcloudrtp_sergio_sarta_filelist (line 32)

>> w = warning('query','last')
w = identifier: 'MATLAB:imagesci:hdf:removalWarningHDFSD'
         state: 'on'
>> id = w.identifier;
>> warning('off',id)
>> lastwarn
ans = 
HDFSW will be removed in a future release. Use MATLAB.IO.HDFEOS.SW instead.
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iv5or6 = 5;   %% AIRS L1B
iv5or6 = 6;   %% AIRS L1C
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% creates an rtp file for ONE granule
%% can be modified for more!

klayers = '/asl/packages/klayers/Bin/klayers_airs';
sarta   = '/asl/packages/sartaV108/Bin/sarta_apr08_m140_wcon_nte';

klayers = '/asl/packages/klayers/Bin/klayers_airs';
sarta   = '/asl/packages/sartaV108/Bin/sarta_apr08_m140_wcon_nte';

klayers = '/asl/packages/klayersV205/BinV201/klayers_airs';
sarta   = '/asl/packages/sartaV108_PGEv6/Bin/sarta_airs_PGEv6_postNov2003';

klayers  = '/home/sergio/git/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/KLAYERS_RTPv221_150levs_80km/klayersV205_0_80km/BinV221/klayers_airs';
sartaCld = '/home/sergio/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/JACvers/bin/jac_airs_l1c_2834_cloudy_may19_prod';
sartaCld = '/home/sergio/git/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/JACvers/bin/jac_airs_l1c_2834_cloudy_jan25_H2020';

%% addpath /home/sergio/MATLABCODE
%% addpath /asl/matlab2012/airs/readers
%% addpath /asl/matlib/aslutil
%% %addpath /asl/matlib/science
%% addpath /home/sergio/MATLABCODE/matlib/science/
%% addpath /asl/matlib/rtptools
%% addpath /asl/matlib/h4tools/
%% addpath /asl/matlib/rtptools/
%% addpath /asl/matlib/gribtools/
%% addpath /asl/matlib/time
%% addpath /home/sergio/MATLABCODE/matlib/clouds/sarta
%% addpath /home/sergio/MATLABCODE
%% addpath /home/sergio/MATLABCODE/matlib/rtp_prod2/emis
%% 
%% % addpath /home/strow/cress/Work/Rtp
%% % addpath /home/strow/Matlab/Grib     WARNING /home/strow/Matlab/Grib/rtpadd_grib_data.m DIFFERENT than /asl/matlib/gribtools/rtpadd_era_data.m
%% % addpath /home/sergio/MATLABCODE/CRIS_HiRes             %% for sergio_fill_ecmwf
%% % addpath /home/strow/Git/rtp_prod2/grib                  %% for fill_ecm
%% % addpath /asl/packages/rtp_prod2/grib
%% addpath /home/sergio/MATLABCODE/matlib/rtp_prod2/grib
%% addpath /home/sergio/MATLABCODE/matlib/rtp_prod2/util
%% addpath /home/sergio/MATLABCODE/NANROUTINES

%% addpath /home/sergio/MATLABCODE/TIME
%% addpath /home/sergio/MATLABCODE/matlib/rtp_prod2/emis/
%% addpath /home/sergio/MATLABCODE/matlib/rtp_prod2/util/
%% 
%% %addpath /asl/rtp_prod2/emis/
%% %addpath /asl/rtp_prod2/util/    
%% %addpath /asl/packages/rtp_prod2/emis/
%% %addpath /asl/packages/rtp_prod2/util/

%%%%%%%%%%%%%%%%%%%%%%%%%

addpath /home/sergio/git/matlabcode
addpath /home/sergio/git/matlabcode/matlibSergio/matlab2012/airs/readers
addpath /home/sergio/git/matlabcode/matlibSergio/matlib/science/
%% addpath /home/sergio/git/matlabcode/matlibSergio/matlib/clouds/sarta
addpath /home/sergio/git/sergio_matlib/matlib/clouds/sarta/

%% addpath /home/strow/cress/Work/Rtp
%% addpath /home/strow/Matlab/Grib     WARNING /home/strow/Matlab/Grib/rtpadd_grib_data.m DIFFERENT than /asl/matlib/gribtools/rtpadd_era_data.m
%% addpath /home/sergio/git/matlabcode/CRIS_HiRes             %% for sergio_fill_ecmwf
%% addpath /home/strow/Git/rtp_prod2/grib                  %% for fill_ecm
%% addpath /asl/packages/rtp_prod2/grib

%addpath /asl/rtp_prod2/emis/
%addpath /asl/rtp_prod2/util/    
%addpath /asl/packages/rtp_prod2/emis/
%addpath /asl/packages/rtp_prod2/util/

addpath /home/sergio/git/rtp_prod2/emis

%addpath /home/sergio/git/matlabcode/matlibSergio/matlib/rtp_prod2/emis    %%%% >>>>>>>>
%addpath /home/sergio/git/matlabcode/matlibSergio/matlib/rtp_prod2/grib    %%%% >>>>>>>>
%addpath /home/sergio/git/matlabcode/matlibSergio/matlib/rtp_prod2/util    %%%% >>>>>>>>
%addpath /home/sergio/git/matlabcode/matlibSergio/matlib/rtp_prod2/emis/   %%%% >>>>>>>>

addpath /home/sergio/git/matlabcode/TIME
addpath /home/sergio/git/matlabcode/matlibSergio/matlib/aslutil
addpath /home/sergio/git/matlabcode/matlibSergio/matlib/science
addpath /home/sergio/git/matlabcode/matlibSergio/matlib/h4tools
addpath /home/sergio/git/matlabcode/matlibSergio/matlib/rtptools
addpath /home/sergio/git/matlabcode/matlibSergio/matlib/gribtools
%% addpath /home/sergio/git/matlabcode/matlibSergio/matlib/TIME

addpath ../GRIB
addpath /home/sergio/git/matlabcode/PLOTTER
    
%%%%%%%%%%%%%%%%%%%%%%%%%

if iv5or6 == 5
  theinds = (1 : 2378)';
else
  theinds = (1 : 2645)';
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%% set sarta exec

run_sarta.clear = +1;
run_sarta.cloud = +1;

if iSlabCld_CumSumStrowORGeorge > 0
  run_sarta.cumsum = 9999;  %% strow pick, cloud at PEAK of wgt fcn
else
  run_sarta.cumsum = -1;  %% aumann pick, cloud at wgt mean of profile
end
run_sarta.klayers_code = klayers;

codeX = 0; %% use default with A. Baran params
codeX = 1; %% use new     with B. Baum, P. Yang params

code0 = '/asl/packages/sartaV108/BinV201/sarta_apr08_m140_iceaggr_waterdrop_desertdust_slabcloud_hg3_wcon_nte';
code1 = '/home/sergio/SARTA_CLOUDY/BinV201/sarta_apr08_m140x_iceGHMbaum_waterdrop_desertdust_slabcloud_hg3';
code1 = sartaCld;

if codeX == 0
  icestr = '_sarta_baran_ice';
  run_sarta.sartacloud_code = code0;
elseif codeX == 1
  icestr = '_sarta_baum_ice';
  run_sarta.sartacloud_code = code1;
else
  error('codeX???')
end

run_sarta.sartaclear_code = run_sarta.sartacloud_code;

%icestr = ['NEWLANDFRAC/cloudy_airs_l1b_ecm' icestr '.'];
if iv5or6 == 5
  icestr = ['cloudy_airs_l1b_ecm' icestr '.'];
elseif iv5or6 == 6
  icestr = ['cloudy_airs_l1c_ecm' icestr '.'];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% typically iaGList == JOB so this loop runs ONCE .. so you control multiple files via the cluster

%% for ixx = 1 : length(iaGlist)
%%   fprintf(1,'main loop : index ixx = %2i of 72 lonbins = %2i \n',ixx)
%%   thetilename = [dir0{JOB} iaGlist(ixx).name]; %% orig ... choose a dir0{JOB},  loop over 1:72 lonbins per JOB
%%   moostr = iaGlist(ixx).name;
%%   moostr = moostr(1:end-3);
%%   fdirOUT = ['/home/sergio/git/matlabcode/QUICKTASKS_TELECON/SuddenStratWarming_SSW/TestRTP/' num2str(JOB,'%03i') '/'];
%%
%% then do " .... now the rest .... "

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% <<< for ixx = 1 : len_dir0  %%% DEFAULT >>>
%
% lonbin 36 has errors right at the beginning
%for ixx = len_dir0 : -1 : 1
%for ixx = len_dir0/2 : -1 : 1
%for ixx = len_dir0/2+1 : len_dir0
%for ixx = 2 : len_dir0

%% for ixx = 1 : len_dir0
%%   fprintf(1,'main loop : index ixx = %2i of len(dir0) = %2i \n',ixx,len_dir0)
%%   lonbin_list = dir([dir0{ixx} '/*.nc']);
%%   thetilename = [dir0{ixx} lonbin_list(JOB).name];  %% new ...  chose a lonbin_JOB), loop over thedir0 list
%%   moostr = lonbin_list(JOB).name;
%%   moostr = moostr(1:end-3);
%%   fdirOUT = ['/home/sergio/git/matlabcode/QUICKTASKS_TELECON/SuddenStratWarming_SSW/TestRTP/' num2str(ixx,'%03i') '/']; 
%%
%% then do " .... now the rest .... "

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ixx = 1 : len_dir0
%for ixx = len_dir0/2 : -1 : 1
  disp(' ')
  disp('new process new process new process new process new process new process')
  disp(' ')
  fprintf(1,'main loop : index ixx = %2i of len(dir0) = %2i \n',ixx,len_dir0)
  lonbin_list = dir([dir0{ixx} '/*.nc']);
  thetilename = [dir0{ixx} lonbin_list(JOB).name];  %% new ...  chose a lonbin_JOB, loop over thedir0 list
  fprintf(1,'ixx = %2i JOB = %3i dir0(ixx) = %s lonbin_list(JOB).name = %s \n',ixx,JOB,dir0{ixx},lonbin_list(JOB).name)
  moostr = lonbin_list(JOB).name;
  moostr = moostr(1:end-3);
  fdirOUT = ['/home/sergio/git/matlabcode/QUICKTASKS_TELECON/SuddenStratWarming_SSW/TestRTP/' num2str(ixx,'%03i') '/'];

  %%% .... now the rest ....
  %%% .... now the rest ....
  %%% .... now the rest ....
  
  fprintf(1,'thetile = %s \n',thetilename)
  fprintf(1,'will save to %s \n',fdirOUT)

  clear p h hattr pattr prof yymmddgg
  if ~exist(fdirOUT)
    mker = ['!mkdir -p ' fdirOUT];
    eval(mker);
    fprintf(1,'made %s \n',fdirOUT)
  end
  
  if iSlabCld_CumSumStrowORGeorge == 1
    fnameOUT= [fdirOUT icestr moostr .rtp'];
  else
    fnameOUT= [fdirOUT icestr moostr  '_cumsum_-1.rtp'];
  end

  eeP = exist(fnameOUT);
  if eeP > 1
    blah = dir(fnameOUT);
    if blah.bytes < 1e8
      fprintf(1,'%s is of size %10i so saying it DNE \n',fnameOUT,blah.bytes)
      eeP = 0;
    end
  end
  fprintf(1,'%s exists 0/1 %2i  \n',fnameOUT,eeP)

%%%%%%%%%%%%%%%%%%%%%%%%%
%% comment this out for running to completion WHICH IS WHAT YOU WANT
%% uncomment if you want to see which files are processed
%end
%return
%if iamanidiot
%%%%%%%%%%%%%%%%%%%%%%%%%
  if eeP == 0
    %!ls -lt /asl/models/ecmwf
    %!ls -lt /asl/models/ecmwf/2020

    fprintf(1,' making %s \n',fnameOUT);
    toucher = ['!touch ' fnameOUT];
    eval(toucher)
  
    a = read_netcdf_lls(thetilename);
    filename =  thetilename;
    iaValidIndices = 1 : a.total_obs;

    iV5or6 = 6;
    p.rlat = a.lat(iaValidIndices)';
    p.rlon = a.lon(iaValidIndices)';
    p.rtime = a.tai93(iaValidIndices)' + offset1958_to_1993;
    p.landfrac = a.land_frac(iaValidIndices)';
    p.solzen   = a.sol_zen(iaValidIndices)';
    p.satzen   = a.sat_zen(iaValidIndices)';
    p.ascflag  = a.asc_flag(iaValidIndices)';            
    p.robs1    = a.rad(:,iaValidIndices);
    p.radqc    = a.rad_qc(iaValidIndices)';

    p.pobs = zeros(size(p.solzen));
    p.upwell = ones(size(p.solzen));
    %p.irinst = AIRSinst*ones(1,nobs);
    %p.findex = grannum*ones(1,nobs);

    scatter_coast(p.rlon,p.rlat,10,rad2bt(1231,p.robs1(1520,:))); colormap jet; title('BT1231(lon,lat)')
    pause(0.1);
    
    pa = {{'profiles','rtime','seconds since 1993'}};
    ha = {{'header','hdf file',filename}};

    h.pfields = 5; % (1=prof + 4=IRobs);
    h.ptype   = 0;

    if iv5or6 == 5
      h.nchan = length(theinds);
      h.ichan = 1:2378;
      h.vchan = f(h.ichan);
    else
      wavoo = load('/home/sergio/git/matlabcode/h2645structure.mat');    
      h.nchan = length(wavoo.h.ichan);
      h.ichan = wavoo.h.ichan;;
      h.vchan = wavoo.h.vchan;
    end

    %%% this is NEW
    %p.landfrac_fromL1B = p.landfrac;
    %p.salti_fromL1B = p.salti;
    %[salti, landfrac] = usgs_deg10_dem(p.rlat, p.rlon);
    %p.landfrac = landfrac;
    %p.salti    = salti;

    clrfields = {'SP','SKT','10U','10V','TCC','CI','T','Q','O3'};
    cldfields = {'SP','SKT','10U','10V','TCC','CI','T','Q','O3',...
                 'CC','CIWC','CLWC'};

    %[h,ha,p,pa] = rtpadd_ecmwf_data(h,ha,p,pa,cldfields); %%% add on ecm
    %[p,h] = sergio_fill_ecmwf(p,h,'/asl/data/ecmwf/',-1);
    %save test_2002_09_08_g044_sergio.mat p h

    which fill_ecmwf
    disp('calling fill_ecmwf')

   %p00 = p

    [p,h] = fill_ecmwf(p,h);

    [xyy,xmm,xdd,xhh] = tai2utcSergio(p.rtime);        %%% <<<<<<<<<<<<<<<<<<<<<<<<<<<<< for SdSM old time
    time_so_far = (xyy-2000) + ((xmm-1)+1)/12;
    co2ppm = 368 + 2.077*time_so_far;  %% 395.6933
    p.co2ppm = co2ppm;
    run_sarta.co2ppm = p.co2ppm;
    fprintf(1,'CLIMATOLOGY co2ppm for FIRST %4i/%2i/%2i = %8.6f ppmv\n',xyy(1),xmm(1),xdd(1),p.co2ppm(1));
    fprintf(1,'CLIMATOLOGY co2ppm for LAST  %4i/%2i/%2i = %8.6f ppmv\n',xyy(end),xmm(end),xdd(end),p.co2ppm(end));

    p0 = p;

    %[h,ha,p,pa] = rtpadd_emis_DanZhou2(h,ha,p,pa);
    %p = Prof_add_emis(p,yymmddgg(1),yymmddgg(2),yymmddgg(3));  %% broken crap by whoever
    %p = rtpadd_emis_DanZhou(h,ha,p,pa);   %% lso totally broken crap
    %[h,ha,p,pa] = rtpadd_emis_wis(h,ha,p,pa);
    p.rlon = wrapTo180(p.rlon);

    if exist('/asl/data/iremis/danz/danz_interpolant.mat')
      [p,pa] = rtp_add_emis(p,pa);
    else  
      disp('no /asl/data/iremis/danz/danz_interpolant.mat so skip rtp_add_emis  .... use constant emis')
      disp('no /asl/data/iremis/danz/danz_interpolant.mat so skip rtp_add_emis  .... use constant emis')
      disp('no /asl/data/iremis/danz/danz_interpolant.mat so skip rtp_add_emis  .... use constant emis')
      p.nemis = ones(size(p.stemp)) * 2;
      p.efreq = [600 3000]' * ones(1,length(p.stemp));
      p.emis  = [0.98 0.98]' * ones(1,length(p.stemp));
      p.rho = (1-p.emis)/pi;
    end
    
    %figure(1)
    %scatter_coast(p.rlon,p.rlat,10,p.nemis); 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [p2] = driver_sarta_cloud_rtp(h,ha,p,pa,run_sarta);

    fnamex = fnameOUT;
    [h,ha,p2x,pa] = rtptrim_sartacloud(h,ha,p2,pa);
    rtpwrite(fnamex,h,ha,p2x,pa)
    fprintf(1,'saved %s \n',fnamex)

    %{
    i1231 = find(h.vchan >= 1231,1);
    tobs1231 = real(rad2bt(h.vchan(i1231),p2x.robs1(i1231,:)));
    tclr1231 = rad2bt(h.vchan(i1231),p2x.sarta_rclearcalc(i1231,:));
    tcld1231 = rad2bt(h.vchan(i1231),p2x.rcalc(i1231,:));

    addpath /home/sergio/MATLABCODE/PLOTTER
    scatter_coast(p2x.rlon,p2x.rlat,10,tobs1231)
    scatter_coast(p2x.rlon,p2x.rlat,10,tobs1231-tcld1231); caxis([-10 +10]); colorbar
    %}   

    %{
    addpath /home/sergio/MATLABCODE/matlib/clouds/sarta/
    addpath /home/sergio/MATLABCODE/matlib/clouds/TCC/
    tcc1 = tcc_method1(p);
    tcc2 = tcc_method2(p);
    [tcc3A,tcc3B] = tcc_method3(p);

    figure(1); scatter_coast(p2x.rlon,p2x.rlat,10,tcc1); colormap jet; title('METHOD 1')
    figure(2); scatter_coast(p2x.rlon,p2x.rlat,10,tcc2); colormap jet; title('METHOD 2')
    figure(3); scatter_coast(p2x.rlon,p2x.rlat,10,tcc3A); colormap jet; title('METHOD 3A')
    figure(4); scatter_coast(p2x.rlon,p2x.rlat,10,tcc3B); colormap jet; title('METHOD 3B')
    figure(5); scatter_coast(p2x.rlon,p2x.rlat,10,p.tcc); colormap jet; title('ORIG')
    
    [Y,I] = sort(p.tcc);
    figure(6); plot(p.tcc(I),tcc1(I),'b',p.tcc(I),tcc2(I),'g',p.tcc(I),tcc3A(I),'r',p.tcc(I),tcc3B(I),'m')
    corr1 = linearcorrelation(p.tcc,tcc1);
    corr2 = linearcorrelation(p.tcc,tcc2);
    corr3A = linearcorrelation(p.tcc,tcc3A);
    corr3B = nanlinearcorrelation(p.tcc,tcc3B);
    [corr1 corr2 corr3A corr3B]    
    %}
    
  else
    fprintf(1,' %s already exists \n',fnameOUT)
  end
end
