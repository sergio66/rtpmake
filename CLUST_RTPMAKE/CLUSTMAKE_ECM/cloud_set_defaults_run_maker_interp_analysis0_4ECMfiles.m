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
    orig_4_ECMfiles_interp_analysis
    %%%%%%%%%%%%%%%%%%%%%%%%%    
    p0 = p;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
p321 = find(p.plevs(:,1) > 321,1);

figure(1); scatter_coast(pClosest.rlon,pClosest.rlat,30,p.gas_1(p321,:));
figure(2); scatter_coast(pClosest.rlon,pClosest.rlat,30,pClosest.gas_1(p321,:));
figure(3); scatter_coast(pClosest.rlon,pClosest.rlat,30,pB1.gas_1(p321,:));
figure(4); scatter_coast(pClosest.rlon,pClosest.rlat,30,pB2.gas_1(p321,:));

rh321      = mixr2rh(p.gas_1(p321,:)*1000,p.plevs(p321,:),p.ptemp(p321,:),1);  %% from IDL routines
irionrh321 = mixr2rh(p.gas_1(p321,:)*1000,p.plevs(p321,:),p.ptemp(p321,:),0);  %% from IDL routines, with SVP a linear mix of I/W (Irion)
rah = p.gas_1(p321,:); %% SH in g/g
rah = rah./(1-rah);       %% mixR in g/g
irionrh321 = mixr2rh(rah*1000,p.plevs(p321,:),p.ptemp(p321,:),0);  %% from IDL routines, with SVP a linear mix of I/W (Irion)
xrh321 = convert_humidity(p.plevs(p321,:),p.ptemp(p321,:),p.gas_1(p321,:),'mixing ratio','relative humidity');  %% from Strow
xrh321 = xrh321*100;
figure(1); clf; scatter_coast(pClosest.rlon,pClosest.rlat,30,irionrh321); caxis([0 120]); colormap jet; title('analysis interp')
axis([115 145 10 40]); colormap(irion_idl_colormap);
    
rh321      = mixr2rh(pClosest.gas_1(p321,:)*1000,pClosest.plevs(p321,:),pClosest.ptemp(p321,:),1);  %% from IDL routines
irionrh321 = mixr2rh(pClosest.gas_1(p321,:)*1000,pClosest.plevs(p321,:),pClosest.ptemp(p321,:),0);  %% from IDL routines, with SVP a linear mix of I/W (Irion)
rah = pClosest.gas_1(p321,:); %% SH in g/g
rah = rah./(1-rah);       %% mixR in g/g
irionrh321 = mixr2rh(rah*1000,pClosest.plevs(p321,:),pClosest.ptemp(p321,:),0);  %% from IDL routines, with SVP a linear mix of I/W (Irion)
xrh321 = convert_humidity(pClosest.plevs(p321,:),pClosest.ptemp(p321,:),pClosest.gas_1(p321,:),'mixing ratio','relative humidity');  %% from Strow
xrh321 = xrh321*100;
figure(2); clf; scatter_coast(pClosest.rlon,pClosest.rlat,30,irionrh321); caxis([0 120]); colormap jet; title('pClosest in time')
axis([115 145 10 40]); colormap(irion_idl_colormap);

[hx,hax,px,pax] = rtpread('/asl/rtp/rtprod_airs/2002/09/06/cloudy_airs_l1b_ecm_sarta_baum_ice.2002.09.06.044.rtp');
p321 = find(px.plevs(:,1) > 321,1)
rh321      = mixr2rh(px.gas_1(p321,:)*1000,px.plevs(p321,:),px.ptemp(p321,:),1);  %% from IDL routines
irionrh321 = mixr2rh(px.gas_1(p321,:)*1000,px.plevs(p321,:),px.ptemp(p321,:),0);  %% from IDL routines, with SVP a linear mix of I/W (Irion)
rah = px.gas_1(p321,:); %% SH in g/g
rah = rah./(1-rah);       %% mixR in g/g
irionrh321 = mixr2rh(rah*1000,px.plevs(p321,:),px.ptemp(p321,:),0);  %% from IDL routines, with SVP a linear mix of I/W (Irion)
xrh321 = convert_humidity(px.plevs(p321,:),px.ptemp(p321,:),px.gas_1(p321,:),'mixing ratio','relative humidity');  %% from Strow
xrh321 = xrh321*100;
figure(3); clf; scatter_coast(pClosest.rlon,pClosest.rlat,30,irionrh321); caxis([0 120]); colormap jet; title('px')
axis([115 145 10 40]); colormap(irion_idl_colormap);

for ip = 1 : 3
  figure(ip); colormap jet
  axis([115 145 10 40]); colormap(irion_idl_colormap);  
end
%}

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
