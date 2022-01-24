addpath /asl/matlab2012/airs/readers
addpath /asl/matlib/aslutil
%addpath /asl/matlib/science
addpath  addpath /home/sergio/MATLABCODE/matlib/science/
addpath /asl/matlib/rtptools
addpath /asl/matlib/h4tools/
addpath /asl/matlib/rtptools/
addpath /asl/matlib/gribtools/
addpath /asl/matlib/time
addpath /home/sergio/MATLABCODE/matlib/clouds/sarta
addpath /home/sergio/MATLABCODE

% addpath /home/strow/cress/Work/Rtp
% addpath /home/strow/Matlab/Grib     WARNING /home/strow/Matlab/Grib/rtpadd_grib_data.m DIFFERENT than /asl/matlib/gribtools/rtpadd_era_data.m
% addpath /home/sergio/MATLABCODE/CRIS_HiRes             %% for sergio_fill_ecmwf
% addpath /home/strow/Git/rtp_prod2/grib                  %% for fill_ecm
% addpath /asl/packages/rtp_prod2/grib
addpath /home/sergio/MATLABCODE/matlib/rtp_prod2/grib

load  /home/motteler/shome/obs_stats/airs_tiling/obs_2019_s6.mat

%%%%%%%%%%%%%%%%%%%%%%%%%
%%% set sarta exec

run_sarta.iWhichInterp = 0;
run_sarta.clear = -1;
run_sarta.cloud = -1;
run_sarta.cumsum = -1;    %% this is "closer" to MRO but since cliuds are at centroid, does not do too well with DCC
run_sarta.cumsum = 9999;  %% larrabee likes this, puts clouds high so does well for DCC

code0 = '/asl/packages/sartaV108/BinV201/sarta_apr08_m140_iceaggr_waterdrop_desertdust_slabcloud_hg3_wcon_nte';
code1 = '/home/sergio/SARTA_CLOUDY/BinV201/sarta_apr08_m140x_iceGHMbaum_waterdrop_desertdust_slabcloud_hg3';
code1 = '/home/chepplew/gitLib/sarta/bin/airs_l1c_2834_cloudy_may19_prod_v3';

%%%%%%%%%%%%%%%%%%%%%%%%%

clear p h hattr pattr prof yymmddgg

hdffile = '/home/sergio/MATLABCODE/airs_l1c_srf_tables_lls_20181205.hdf';   % what he gave in Dec 2018
vchan2834 = hdfread(hdffile,'freq');
f = vchan2834;
load sarta_chans_for_l1c.mat
theinds2645 = ichan;
f2645 = f(ichan);

p.rlat = lat_list';
p.rlon = lon_list';
addpath /home/sergio/MATLABCODE/TIME
p.rtime = tai93_list' + offset1958_to_1993;

iJump = 10000;
iJump = 1000;
p1.rlat = p.rlat(1:iJump:length(lat_list));
p1.rlon = p.rlon(1:iJump:length(lat_list));
p1.rtime = p.rtime(1:iJump:length(lat_list));

pa = {{'profiles','rtime','seconds since 1993'}};
ha = {{'header','hdf file','miaow'}};

h.pfields=1; % (1=prof + 4=IRobs);

h.nchan = length(theinds2645);
h.ichan = theinds2645;
h.vchan = f2645;

[salti, landfrac] = usgs_deg10_dem(p.rlat, p.rlon);
  p.landfrac = landfrac; 
  p.salti    = salti;
[salti, landfrac] = usgs_deg10_dem(p1.rlat, p1.rlon);
  p1.landfrac = landfrac; 
  p1.salti    = salti;

clrfields = {'SP','SKT','10U','10V','TCC','CI','T','Q','O3'};
cldfields = {'SP','SKT','10U','10V','TCC','CI','T','Q','O3',...
                 'CC','CIWC','CLWC'};

%[h,ha,p,pa] = rtpadd_era_data(h,ha,p,pa,cldfields); %%% add on era
tic
[p1,h] = fill_era(p1,h);
toc
[p,h] = fill_era(p,h);
toc
error('nana')

addpath /home/sergio/MATLABCODE/TIME
[xyy,xmm,xdd,xhh] = tai2utcSergio(p.rtime);        %%% <<<<<<<<<<<<<<<<<<<<<<<<<<<<< for SdSM old time
time_so_far = (xyy-2000) + ((xmm-1)+1)/12;
co2ppm = 368 + 2.077*time_so_far;  %% 395.6933
p.co2ppm = co2ppm;
fprintf(1,'CLIMATOLOGY co2ppm for FIRST %4i/%2i/%2i = %8.6f ppmv\n',xyy(1),xmm(1),xdd(1),p.co2ppm(1));
fprintf(1,'CLIMATOLOGY co2ppm for LAST  %4i/%2i/%2i = %8.6f ppmv\n',xyy(end),xmm(end),xdd(end),p.co2ppm(end));
    
addpath /home/sergio/MATLABCODE/matlib/rtp_prod2/emis/
addpath /home/sergio/MATLABCODE/matlib/rtp_prod2/util/    
p1.satzen = zeros(size(p1.stemp));
p1.wspeed = ones(size(p1.stemp)) * 10;

p1.rlon = wrapTo180(p1.rlon);
[p1,pa] = rtp_add_emis(p1,pa);

p.rlon = wrapTo180(p.rlon);
[p,pa] = rtp_add_emis(p,pa);
    
[pX] = driver_sarta_cloud_rtp(h,ha,p,pa,run_sarta);
[p1X] = driver_sarta_cloud_rtp(h,ha,p1,pa,run_sarta);
error(';lskg;lkgs')
