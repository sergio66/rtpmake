disp('checking really high col water from clust_loop_make16day_tile_center.m for2014/01/last week')

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
[h,p] = subset_rtp_allcloudfields(h,p,[],[],1939);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f1 = '/asl/models/era/2007/11/20071125_sfc.nc';
%f2 = '/asl/models/era/2007/11/20071125_lev.nc';

s1 = read_netcdf_lls(f1);
%s2 = read_netcdf_lls(f2);

yres = 0.75;
[mm,nn,oo] = size(s1.skt); fprintf(1,'ERA-I size(skt) = %4i %4i %4i \n',mm,nn,oo);  %% 480 x 241 x 4

yres = 0.75;
  [(360-0.75)/(mm-1) 180/(nn-1)]  %% remember N/S is unique hence (180)/(nn-1), while E-W wraps around hence (360-0.75)/(mm-1)
  era_rlon = (1:480); era_rlon = era_rlon-1; era_rlon = 0 + era_rlon*0.75;
  era_rlat = (1:241); era_rlat = era_rlat-1; era_rlat = -90 + era_rlat*0.75;
  [era_Y,era_X] = meshgrid(era_rlat,era_rlon); 

junkx = find(era_rlon <= p.rlon); junkx = junkx(end); [junkx era_rlon(junkx) p.rlon era_rlon(junkx+1)]
junky = find(era_rlat <= p.rlat); junky = junky(end); [junky era_rlat(junky) p.rlat era_rlat(junky+1)]
[h,p] = replicate_rtp_headprof(h,p,1,5);
%
%  1         2
%    orig=3
%  4         5
%
p.rlon(1) = era_rlon(junkx); p.rlon(2) = era_rlon(junkx+1); p.rlon(4) = era_rlon(junkx);   p.rlon(5) = era_rlon(junkx+1); 
p.rlat(1) = era_rlat(junky); p.rlat(2) = era_rlat(junky);   p.rlat(4) = era_rlat(junky+1); p.rlat(5) = era_rlat(junky+1); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(1,'p.rtime = %25.22e \n',p.rtime(1))
[yyyJ,mmmJ,dddJ,hhhJ] = tai2utcSergio(p.rtime(1))

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
iRemovePalts = input('remove palts (-1,no) (+1,yes BEST BEST BEST) : ');
if iRemovePalts > 0
  p = rmfield(p,'palts');  %% new
end

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

addpath /home/sergio/MATLABCODE/matlib/rtp_prod2/emis
addpath /home/sergio/MATLABCODE/matlib/rtp_prod2/util
addpath /home/sergio/MATLABCODE/matlib/rtp_prod2/util/time
  p.rlon = wrapTo180(p.rlon);
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
if exist('/asl/s1/sergio/rtp/rtp_airicrad_v6/2014/01/29/cloudy_airs_l1c_era_sarta_baum_ice.2014.01.29.031_WRONG_RTIME.rtp')
  [hfake,~,pfake,~] = rtpread('/asl/s1/sergio/rtp/rtp_airicrad_v6/2014/01/29/cloudy_airs_l1c_era_sarta_baum_ice.2014.01.29.031_WRONG_RTIME.rtp');
  rtpwrite('junk.ip.rtp',hfake,ha,pfake,pa);
  sss = ['!' klayers ' fin=junk.ip.rtp fout=junk.op.rtp >& ugh']; eval(sss)
  [hjunk2,~,pjunk2,~] = rtpread('junk.op.rtp');
  mmw2 = mmwater_rtp(hjunk2,pjunk2);
  boo2 = find(mmw2 == max(mmw2),1); mmw2(boo2);
  fprintf(1,'mmw timy (%8.4f mm) at pixel %5i %4i/%2i/%2i %8.6f at rlon/lat %8.4f %8.4f \n',mmw2(boo2),boo2,yyyy,mmmm,dddd,hhhh,pjunk2.rlon(boo),pjunk2.rlat(boo))
  scatter_coast(pjunk2.rlon,pjunk2.rlat,10,mmw2); hold on; plot(p2x.rlon(3),p2x.rlat(3),'kx','linewidth',3,'markersize',5);

pN = pjunk.plevs(1:100,:)-pjunk.plevs(2:101,:);
pD = log(pjunk.plevs(1:100,:)./pjunk.plevs(2:101,:));
pjunk.plays = pN./pD;

pN = pjunk2.plevs(1:100,:)-pjunk2.plevs(2:101,:);
pD = log(pjunk2.plevs(1:100,:)./pjunk2.plevs(2:101,:));
pjunk2.plays = pN./pD;

[pjunk.spres(boo) pjunk2.spres(boo2)]
lays = pjunk.nlevs(boo)-1; lays = 1:lays;
lays2 = pjunk2.nlevs(boo2)-1; lays2 = 1:lays2;

figure(1); subplot(121); loglog(p2x.gas_1(:,boo),p2x.plevs(:,boo),'b.-',pfake.gas_1(:,boo2),pfake.plevs(:,boo2),'r.-'); set(gca,'ydir','reverse'); ylim([0.1 1010])
  hl = legend('wierd','normal','location','best','fontsize',8); title('WV levels'); 
figure(1); subplot(122); semilogy(p2x.ptemp(:,boo),p2x.plevs(:,boo),'b.-',pfake.ptemp(:,boo2),pfake.plevs(:,boo2),'r.-'); set(gca,'ydir','reverse');ylim([0.1 1010]); xlim([180 300])
  hl = legend('wierd','normal','location','best','fontsize',8); title('T levels'); 

figure(2); subplot(121); loglog(pjunk.gas_1(lays,boo),pjunk.plays(lays,boo),'b.-',pjunk2.gas_1(lays2,boo2),pjunk2.plays(lays2,boo2),'r.-'); set(gca,'ydir','reverse'); ylim([0.005 1010])
  hl = legend('wierd','normal','location','best','fontsize',8); title('WV layers'); 
figure(2); subplot(122); semilogy(pjunk.ptemp(lays,boo),pjunk.plays(lays,boo),'b.-',pjunk2.ptemp(lays2,boo2),pjunk2.plays(lays2,boo2),'r.-'); set(gca,'ydir','reverse');ylim([0.005 1010]); xlim([180 300])
  hl = legend('wierd','normal','location','best','fontsize',8); title('T layers'); 



figure(1); subplot(121); loglog(p2x.gas_1(:,boo),p2x.plevs(:,boo),'b.-',pfake.gas_1(:,boo2),pfake.plevs(:,boo2),'r.-'); set(gca,'ydir','reverse'); ylim([700 1010])
  hl = legend('wierd','normal','location','best','fontsize',8); title('WV levels'); 
figure(1); subplot(122); semilogy(p2x.ptemp(:,boo),p2x.plevs(:,boo),'b.-',pfake.ptemp(:,boo2),pfake.plevs(:,boo2),'r.-'); set(gca,'ydir','reverse');ylim([700 1010]); 
  hl = legend('wierd','normal','location','best','fontsize',8); title('T levels'); 

figure(2); subplot(121); loglog(pjunk.gas_1(lays,boo),pjunk.plays(lays,boo),'b.-',pjunk2.gas_1(lays2,boo2),pjunk2.plays(lays2,boo2),'r.-'); set(gca,'ydir','reverse'); ylim([700 1010])
  hl = legend('wierd','normal','location','best','fontsize',8); title('WV layers'); 
figure(2); subplot(122); semilogy(pjunk.ptemp(lays,boo),pjunk.plays(lays,boo),'b.-',pjunk2.ptemp(lays2,boo2),pjunk2.plays(lays2,boo2),'r.-'); set(gca,'ydir','reverse'); ylim([700 1010]); 
  hl = legend('wierd','normal','location','best','fontsize',8); title('T layers'); 


end


