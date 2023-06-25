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
%addpath /asl/matlib/science
addpath /home/sergio/MATLABCODE/matlib/science/
addpath /asl/matlib/rtptools
addpath /asl/matlib/h4tools/
addpath /asl/matlib/rtptools/
addpath /asl/matlib/gribtools/
addpath /asl/matlib/time
addpath /home/sergio/MATLABCODE/matlib/clouds/sarta
addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/matlib/rtp_prod2/emis

% addpath /home/strow/cress/Work/Rtp
% addpath /home/strow/Matlab/Grib     WARNING /home/strow/Matlab/Grib/rtpadd_grib_data.m DIFFERENT than /asl/matlib/gribtools/rtpadd_era_data.m
% addpath /home/sergio/MATLABCODE/CRIS_HiRes             %% for sergio_fill_ecmwf
% addpath /home/strow/Git/rtp_prod2/grib                  %% for fill_ecm
% addpath /asl/packages/rtp_prod2/grib
addpath /home/sergio/MATLABCODE/matlib/rtp_prod2/grib
addpath /home/sergio/MATLABCODE/matlib/rtp_prod2/util

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
if iv5or6 == 5
  icestr = ['cloudy_airs_l1b_ecm' icestr '.'];
elseif iv5or6 == 6
  icestr = ['cloudy_airs_l1c_ecm' icestr '.'];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% typically iaGList == JOB so this loop runs ONCE .. so you control multiple files via the cluster

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
    !ls -lt /asl/models/ecmwf
    !ls -lt /asl/models/ecmwf/2020

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

    read_in_L1B_or_L1C

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

    p.pobs = zeros(size(p.solzen));
    p.upwell = ones(size(p.solzen));
    %p.irinst = AIRSinst*ones(1,nobs);
    %p.findex = grannum*ones(1,nobs);

    plot(p.rlon,p.rlat,'.')

    pa = {{'profiles','rtime','seconds since 1993'}};
    ha = {{'header','hdf file',filename}};

    h.pfields = 5; % (1=prof + 4=IRobs);
    h.ptype   = 0;

    if iv5or6 == 5
      h.nchan = length(theinds);
      h.ichan = 1:2378;
      h.vchan = f(h.ichan);
    else
      h.nchan = length(theinds2645);
      h.ichan = theinds2645;;
      h.vchan = f2645;
    end

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

which fill_ecmwf
disp('calling fill_ecmwf')

%p00 = p

    [p,h] = fill_ecmwf(p,h);

addpath /home/sergio/MATLABCODE/TIME
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
    %addpath /asl/rtp_prod2/emis/
    %addpath /asl/rtp_prod2/util/    
    %addpath /asl/packages/rtp_prod2/emis/
    %addpath /asl/packages/rtp_prod2/util/

    addpath /home/sergio/MATLABCODE/matlib/rtp_prod2/emis/
    addpath /home/sergio/MATLABCODE/matlib/rtp_prod2/util/

    p.rlon = wrapTo180(p.rlon);
    [p,pa] = rtp_add_emis(p,pa);

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
addpath /home/sergio/MATLABCODE/NANROUTINES
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
