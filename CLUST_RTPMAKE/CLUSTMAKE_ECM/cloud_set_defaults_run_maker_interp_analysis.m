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

addpath /home/sergio/MATLABCODE

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

addpath /asl/matlab2012/airs/readers
addpath /asl/matlib/aslutil
addpath /asl/matlib/science
addpath /asl/matlib/rtptools
addpath /asl/matlib/h4tools/
addpath /asl/matlib/rtptools/
addpath /asl/matlib/gribtools/
addpath /home/sergio/MATLABCODE/matlib/clouds/sarta
addpath /home/sergio/MATLABCODE/TIME
addpath /home/sergio/MATLABCODE/PLOTTER
addpath /home/sergio/MATLABCODE/COLORMAP

%addpath /home/strow/cress/Work/Rtp
%addpath /home/strow/Matlab/Grib

%addpath /home/sergio/MATLABCODE/CRIS_HiRes             %% for sergio_fill_ecmwf
%addpath /home/strow/Git/rtp_prod2/grib                  %% for fill_ecm 
%addpath /asl/rtp_prod2/grib/                           %% for fill_ecmwf

%addpath  /asl/packages/rtp_prod2/grib
%addpath  /home/sbuczko1/git/rtp_prod2/grib

addpath /home/strow/git/rtp_prod2/grib
addpath /home/sergio/MATLABCODE/RTPMAKE/CLUST_RTPMAKE/GRIB
addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS/Strow_humidity/convert_humidity/

if iv5or6 == 5
  theinds = (1 : 2378)';
else
  theinds = (1 : 2645)';
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%% set sarta exec

run_sarta.clear = +1;
run_sarta.cloud = +1;
run_sarta.cumsum = 9999;

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

%icestr = ['NEWLANDFRAC/cloudy_airs_l1b_ecm' icestr '.'];
%icestr = ['interp_analysis_cloudy_airs_l1b_ecm' icestr '.'];
if iv5or6 == 5
  icestr = ['interp_analysis_cloudy_airs_l1b_ecm' icestr ];
elseif iv5or6 == 6
  icestr = ['interp_analysis_cloudy_airs_l1c_ecm' icestr ];
end
icestr = [icestr '_timeoffset_' num2str(iTimeOffset,'%04d') '_'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% typically iaGList == JOB so this loop runs ONCE .. so you control muktiple files via the cluster

for ixx = 1 : length(iaGlist)
  ix = iaGlist(ixx);
  clear p h hattr pattr prof yymmddgg

  yymmddgg = [yymmdd0 ix];

  ystr = num2str(yymmddgg(1));
  mstr = num2str(yymmddgg(2),'%02d');
  dstr = num2str(yymmddgg(3),'%02d');
  gstr = num2str(yymmddgg(4),'%03d');

  % fdirOUT = ['/asl/data/rtprod_airs/' ystr '/' mstr '/' dstr '/'];
  fdirOUT = ['/asl/rtp/rtprod_airs/' ystr '/' mstr '/' dstr '/'];
  if iv5or6 == 5
    fdirOUT = ['/asl/s1/sergio/rtp/rtp_airibrad_v5/' ystr '/' mstr '/' dstr '/'];  %% till 2018
  elseif iv5or6 == 6
    fdirOUT = ['/asl/s1/sergio/rtp/rtp_airicrad_v6/' ystr '/' mstr '/' dstr '/'];  %% after 2018
  end

  if ~exist(fdirOUT)
    mker = ['!mkdir -p ' fdirOUT];
    eval(mker);
    fprintf(1,'made %s \n',fdirOUT)
  end

  if iSlabCld_CumSumStrowORGeorge == 1
    fnameOUT= [fdirOUT icestr ystr '.' mstr '.' dstr '.' gstr '.rtp'];
  else
    fnameOUT= [fdirOUT icestr ystr '.' mstr '.' dstr '.' gstr '_cumsum_-1.rtp'];
  end

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

    if iv5or6 == 5
      %% L1B
      % filename = ['/strowdataN/data/airs/Aqua_AIRS_Level1/AIRIBRAD.005/' ystr '/'];
      filename = ['/asl/data/airs/AIRIBRAD/' ystr '/'];
      filename = [filename num2str(days_so_far,'%03d') '/'];
      dir0 = filename;
      filename = [filename 'AIRS.' ystr '.' mstr '.' dstr '.' gstr];
      filename = [filename '.L1B.AIRS_Rad.v5*.hdf'];
    elseif iv5or6 == 6
      %% L1C
      filename = ['/asl/data/airs/L1C/' ystr '/'];
      filename = ['/asl/data/airs/L1C_v672/' ystr '/'];
      if yymmddgg(1) <= 2021
        if days_so_far <= 266
          filename = ['/asl/airs/l1c_v672/' ystr '/'];
        else
          filename = ['/asl/airs/l1c_v674/' ystr '/'];
        end
      elseif yymmddgg(1) > 2021
        filename = ['/asl/airs/l1c_v674/' ystr '/'];
      end
      filename = [filename num2str(days_so_far,'%03d') '/'];

      dir0 = filename;
      filename = [filename 'AIRS.' ystr '.' mstr '.' dstr '.' gstr];
      filename = [filename '.L1C.AIRS_Rad.v6*.hdf'];
    end

    thedir = dir(filename);
    iSimulateData = -1;
    if length(thedir) == 1
      fname = [dir0 thedir.name];
    else
      fprintf(1,'%s \n',filename);
      disp('file does not exist');

      %% excess wet bulb, 2019_06_23
      xjunk = '/home/sergio/MATLABCODE/WetBulbTemperatures/usa_2019_06_23.mat';  

      %% MLS by Werner say a big system and I am STRETCHING it
      xjunk = '/home/sergio/MATLABCODE/WetBulbTemperatures/usa_2019_08_27.mat';

      %% https://earthsky.org/earth/study-predicts-deadly-heat-in-persian-gulf/ hot day in the Persian Gulf
      xjunk = '/home/sergio/MATLABCODE/WetBulbTemperatures/middle_east_2015_07_31.mat';
      xjunk = '/home/sergio/MATLABCODE/WetBulbTemperatures/middle_east_2020_08_23.mat';
      xjunk = '/home/sergio/MATLABCODE/WetBulbTemperatures/middle_east_2020_07_29.mat';
      xjunk = '/home/sergio/MATLABCODE/WetBulbTemperatures/middle_east_2020_07_29pm2.mat';

      %% excess wet bulb, 2020_08_23
      xjunk = '/home/sergio/MATLABCODE/WetBulbTemperatures/usa_2020_08_21.mat';
      xjunk = '/home/sergio/MATLABCODE/WetBulbTemperatures/usa_2020_08_23.mat';    %% hmm, not as excessive as I thought!
      xjunk = '/home/sergio/MATLABCODE/WetBulbTemperatures/usa_2020_08_20.mat';    %% this one looked better
      xjunk = '/home/sergio/MATLABCODE/WetBulbTemperatures/usa_2020_08_23pm2.mat'; %% hmm, not as excessive as I thought!

      xjunk = '/home/sergio/MATLABCODE/WetBulbTemperatures/india_2019_05_11.mat';
      xjunk = '/home/sergio/MATLABCODE/WetBulbTemperatures/india_2019_05_11pm2.mat';


      %%%%%%%%%%%%%%%%%%%%%%%%%

      strjunk = ['load ' xjunk '  as needed ? ']; fprintf(1,'xjunk = %s \n',xjunk);
      %iSimulateData = input(strjunk);
      iSimulateData = 1

      if iSimulateData < 0
        return
      end
    end

    if iv5or6 == 5
      [a,b,c] = sdload_quiet(fname);
    elseif iv5or6 == 6
      %% find v6_readl2cc.m  /asl/*/rtp_prod2/airs/readers
      addpath /home/sergio/MATLABCODE/matlib/rtp_prod2/airs/readers
      addpath /asl/matlib/time

      if iSimulateData > 0
        fprintf(1,'loading %s \n',xjunk);
        pjunk = load(xjunk); 

        gdata.rtime = pjunk.p0.rtime;
        gdata.rlat = pjunk.p0.rlat;
        gdata.rlon = pjunk.p0.rlon;
        gdata.solazi = ones(size(pjunk.p0.rlon)) * 00;
        gdata.solzen = ones(size(pjunk.p0.rlon)) * 150;
          wonk = length(gdata.solzen);
          gdata.solzen(1:wonk/2) = 40; 
        gdata.satzen = ones(size(pjunk.p0.rlon)) * 22;
        gdata.satazi = ones(size(pjunk.p0.rlon)) * 0;
        gdata.robs1 = zeros(2645,length(gdata.rtime));
       
        [gdata.salti, gdata.landfrac] = usgs_deg10_dem(gdata.rlat, gdata.rlon);
      else
        [eq_x_tai, f, gdata, attr, opt] = read_airicrad(fname);  % Steve
        f0_2645 = f;
      end

hdffile = '/home/sergio/MATLABCODE/airs_l1c_srf_tables_lls_20181205.hdf';   % what he gave in Dec 2018
vchan2834 = hdfread(hdffile,'freq');
f = vchan2834;
load sarta_chans_for_l1c.mat
theinds2645 = ichan;
f2645 = f(ichan);

%      a = read_airs_l1c(fname);   %% Chris Hepplewhite
%      theinds2645 = cell2mat(a.chanID); theinds2645 = theinds2645';
%      f2645   = cell2mat(a.freq);   f2645   = f2645';
% plot(f2645)
% plot(chanID)

    end

    if iv5or6 == 5
      p.rlat = [a.Latitude(:)'];
      p.rlon = [a.Longitude(:)'];
      p.rtime = [a.Time(:)'];

      %[meantime, f, prof] = readl1b_all(fname);  %% has the old rtime 1993
      [meantime, f, prof] = readl1b_all(fname);  %% has the new rtime 1958    
      p = prof;
    elseif iv5or6 == 6
      p = gdata;
    end

    p.rtime = p.rtime + iTimeOffset*60;  %% <<<<<<< add on timeoffset , convert to seconds >>>>>>>>>>

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

    %[h,ha,p,pa] = rtpadd_ecmwf_data(h,ha,p,pa,cldfields); %%% add on ecm
    %[p,h] = sergio_fill_ecmwf(p,h,'/asl/data/ecmwf/',-1);
    %save test_2002_09_08_g044_sergio.mat p h

    %%%%%%%%%%%%%%%%%%%%%%%%%
    %orig_8_ECMfiles_interp_analysis
    new_8_ECMfiles_interp_analysis
    %%%%%%%%%%%%%%%%%%%%%%%%%

    p0 = p;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %[h,ha,p,pa] = rtpadd_emis_DanZhou2(h,ha,p,pa);
    %p = Prof_add_emis(p,yymmddgg(1),yymmddgg(2),yymmddgg(3));  %% broken crap by whoever
    %p = rtpadd_emis_DanZhou(h,ha,p,pa);   %% lso totally broken crap
    %[h,ha,p,pa] = rtpadd_emis_wis(h,ha,p,pa);
    %addpath /asl/rtp_prod2/emis/
    %addpath /asl/rtp_prod2/util/    
    %addpath /asl/packages/rtp_prod2/emis/
    %addpath /asl/packages/rtp_prod2/util/
    %addpath /asl/rtp_prod2/emis/
    %addpath /asl/rtp_prod2/util/
 
    addpath /home/sergio/MATLABCODE/matlib/rtp_prod2/emis/
    addpath /home/sergio/MATLABCODE/matlib/rtp_prod2/util/

    p.rlon = wrapTo180(p.rlon);
    [p,pa] = rtp_add_emis(p,pa);

    %figure(1)
    %scatter_coast(p.rlon,p.rlat,10,p.nemis); 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [p2] = driver_sarta_cloud_rtp(hB2,ha,p,pa,run_sarta);

    fnamex = fnameOUT;
    [h,ha,p2x,pa] = rtptrim_sartacloud(h,ha,p2,pa);
    rtpwrite(fnamex,h,ha,p2x,pa)
    fprintf(1,'saved to %s \n',fnamex)
  else
    fprintf(1,' %s already exists \n',fnameOUT)
  end
end
