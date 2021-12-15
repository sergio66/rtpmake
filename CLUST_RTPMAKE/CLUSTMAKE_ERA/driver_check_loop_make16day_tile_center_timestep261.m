addpath /home/sergio/MATLABCODE
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

addpath /home/sergio/MATLABCODE/RTPMAKE/CLUST_RTPMAKE/GRIB

disp('checking really high col water from clust_loop_make16day_tile_center.m for 2014/01/last week')
disp('see ../CLUSTMAKE_ERA5/compare_ERA5_daily137_ERA_daily60.m')

%%%%%%%%%%%%%%%%%%%%%%%%%
%%% set sarta exec

code1 = '/home/sergio/SARTA_CLOUDY/BinV201/sarta_apr08_m140x_iceGHMbaum_waterdrop_desertdust_slabcloud_hg3';
code1 = '/home/sergio/SARTA_CLOUDY/BinV201/xsarta_apr08_m140_iceGHMbaum_waterdrop_desertdust_slabcloud_hg3';
code1 = '/home/chepplew/gitLib/sarta/bin/airs_l1c_2834_cloudy_may19_prod_v3';
run_sarta.sartaclear_code = code1;
run_sarta.sartacloud_code = code1;

run_sarta.clear = +1;
run_sarta.cloud = +1;
run_sarta.cumsum = -1;    %% this is "closer" to MRO but since cliuds are at centroid, does not do too well with DCC
run_sarta.cumsum = 9999;  %% larrabee likes this, puts clouds high so does well for DCC

era0  = load('/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA/Tile_Center/ASC/era_tile_center_timestep_261_closestINtime.mat');
h = era0.hnew_ip;
p = era0.pnew_ip;
  p = rmfield(p,'stemp');
  p = rmfield(p,'ptemp');
  p = rmfield(p,'gas_1');
  p = rmfield(p,'gas_3');
  p = rmfield(p,'cc');
  p = rmfield(p,'ciwc');
  p = rmfield(p,'clwc');
  p = rmfield(p,'tcc');
  p = rmfield(p,'cpsize');
  p = rmfield(p,'cngwat');
  p = rmfield(p,'cfrac');
  p = rmfield(p,'cprtop');
  p = rmfield(p,'cprbot');
  p = rmfield(p,'cpsize2');
  p = rmfield(p,'cngwat2');
  p = rmfield(p,'cfrac2');
  p = rmfield(p,'cprtop2');
  p = rmfield(p,'cprbot2');
  p = rmfield(p,'cfrac12');

[p,h] = fill_era_interp(p,h);                           %%% add on era the NEW WAY

addpath /home/sergio/MATLABCODE/TIME
[xyy,xmm,xdd,xhh] = tai2utcSergio(p.rtime);        %%% <<<<<<<<<<<<<<<<<<<<<<<<<<<<< for SdSM old time
time_so_far = (xyy-2000) + ((xmm-1)+1)/12;
co2ppm = 368 + 2.077*time_so_far;  %% 395.6933
p.co2ppm = co2ppm;
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

   pa = {{'profiles','rtime','seconds since 1993'}};
   pa = {{'profiles','rtime','seconds since 1958'}};
   ha = {{'header','hdf file','check 2014/01/28'}};

[p,pa] = rtp_add_emis(p,pa);

%figure(1)
%scatter_coast(p.rlon,p.rlat,10,p.nemis); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[p2] = driver_sarta_cloud_rtp(h,ha,p,pa,run_sarta);

[h,ha,p2x,pa] = rtptrim_sartacloud(h,ha,p2,pa);

klayers = '/asl/packages/klayersV205/BinV201/klayers_airs';
rtpwrite('junk.ip.rtp',h,ha,p2x,pa);
sss = ['!' klayers ' fin=junk.ip.rtp fout=junk.op.rtp >& ugh']; eval(sss)
[hjunk,~,pjunk,~] = rtpread('junk.op.rtp');
mmw = mmwater_rtp(hjunk,pjunk);
scatter_coast(pjunk.rlon,pjunk.rlat,10,mmw)

%%% find where mmw = outrageously huge!!!!!
boo = find(mmw == max(mmw),1); mmw(boo);
[yyyy,mmmm,dddd,hhhh] = tai2utcSergio(p2.rtime(boo));
fprintf(1,'mmw huge (%8.4f mm) at pixel %5i %4i/%2i/%2i %8.6f at rlon/lat %8.4f %8.4f \n',mmw(boo),boo,yyyy,mmmm,dddd,hhhh,p2.rlon(boo),p2.rlat(boo))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
check_loop_make16day_tile_center_timestep261_subset1939
