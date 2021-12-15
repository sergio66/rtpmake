addpath /asl/matlib/aslutil
addpath /asl/matlib/science
addpath /asl/matlib/rtptools
addpath /asl/matlib/h4tools/
addpath /asl/matlib/rtptools/
addpath /asl/matlib/gribtools/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% see /home/sergio/klayersV205/Src_rtpV105_100layercloudamountsize
%% klayers = '/home/sergio/klayersV205/Bin/klayers_airs_x_testme';

%% see /home/sergio/klayersV205/Src_rtpV201_100layercloudamountsize
klayers = '/home/sergio/klayersV205/BinV201/klayers_airs_x_testmeCLOUDversion';

%% see /home/sergio/SARTA_CLOUDY/v107_Src_pclsam_ctype_100layer
sarta = '/home/sergio/SARTA_CLOUDY/Bin/sarta_dec05_iceaggr_waterdrop_volz_1_May07_100layer_testme';

%% see /home/sergio/SARTA_CLOUDY/v108_Src_pclsam_ctype_100layer_rtpV201
sarta = '/home/sergio/SARTA_CLOUDY/Bin/sarta_dec05_iceaggr_waterdrop_volz_1_May07_100layer_testme';

%% see /home/sergio/SARTA_CLOUDY/v108_Src_rtpV201_pclsam_slabcloud_hg3_100layerNEW
sarta = '/home/sergio/SARTA_CLOUDY/BinV201/sarta_apr08_m140x_iceGHMbaum_waterdrop_desertdust_slabcloud_hg3_100layerNEW';

%[h0,ha,p0,pa] = rtpread('/asl/s1/sergio/LES/RTP/rico.50446.0.nc_perturbWV_10.rtp');
%[h0,ha,p0,pa] = rtpread('/asl/data/rtprod_airs/2011/03/11/cloudy_airs_l1b_era_sarta_baum_ice.2011.03.11.240.rtp');
[h0,ha,p0,pa] = rtpread('/asl/data/rtprod_airs/2011/03/11/cloudy_airs_l1b_era_sarta_baum_ice.2011.03.11.230.rtp');

h0.ngas = h0.ngas + 2;
h0.glist = [h0.glist;  [201 202]'];
h0.gunit = [h0.gunit; h0.gunit];

p0.ctype   = ones(size(p0.stemp)) * 201;
p0.ctype2  = ones(size(p0.stemp)) * 101;
p0.gas_201 = p0.ciwc;  sum201 = nansum(p0.ciwc);
p0.gas_202 = p0.clwc;  sum202 = nansum(p0.clwc);
p0.cfrac  = zeros(size(p0.cfrac)); oo = find(sum201 > 0); p0.cfrac(oo) = 1;
p0.cfrac2 = zeros(size(p0.cfrac2)); oo = find(sum202 > 0); p0.cfrac2(oo) = 1;
p0.cfrac12 = max(p0.cfrac,p0.cfrac2);

p0.cprtop = -9999 * ones(size(p0.cfrac));
p0.cprbot = -9999 * ones(size(p0.cfrac));
p0.cprtop2 = -9999 * ones(size(p0.cfrac));
p0.cprbot2 = -9999 * ones(size(p0.cfrac));

p0.cprtop = 0 * ones(size(p0.cfrac));
p0.cprbot = p0.spres;
p0.cprtop2 = 0 * ones(size(p0.cfrac));
p0.cprbot2 = p0.spres;

[h0,ha,p0,pa] = rtpadd_emis_wis(h0,ha,p0,pa);
p0 = rmfield(p0,'cprtop');
p0 = rmfield(p0,'cprtop2');
p0 = rmfield(p0,'cprbot');
p0 = rmfield(p0,'cprbot2');

iJump = 10;
iJump = 2;
iJump = 100;
[h0,p0] = subset_rtp_clouds(h0,p0,[],[],1:iJump:length(p0.stemp));
p0.ctype   = ones(size(p0.stemp)) * 201;
p0.ctype2  = ones(size(p0.stemp)) * 101;
p0.gas_201 = p0.ciwc;  sum201 = nansum(p0.ciwc);
p0.gas_202 = p0.clwc;  sum202 = nansum(p0.clwc);
p0.cfrac  = zeros(size(p0.cfrac)); oo = find(sum201 > 0); p0.cfrac(oo) = 1;
p0.cfrac2 = zeros(size(p0.cfrac2)); oo = find(sum202 > 0); p0.cfrac2(oo) = 1;
p0.cfrac12 = max(p0.cfrac,p0.cfrac2)  
  oo = find(sum201 <= 0); p0.cfrac12(oo) = p0.cfrac(oo);
  oo = find(sum202 <= 0); p0.cfrac12(oo) = p0.cfrac2(oo);

rtpwrite('xtest.ip.rtp',h0,ha,p0,pa);

klayerser = ['!' klayers ' fin=xtest.ip.rtp fout=xtest.op.rtp nwant=10 listg=1,2,3,4,5,6,9,12,201,202'];
sartaer   = ['!' sarta   ' fin=xtest.op.rtp fout=xtest.rp.rtp'];

disp('klayers')
eval(klayerser)
[h1,ha,p1,pa] = rtpread('xtest.op.rtp');
plot(p1.cngwat2+p1.cngwat,sum(p1.gas_201)+sum(p1.gas_202),'ko',p1.cngwat2,sum(p1.gas_202),'bx',p1.cngwat,sum(p1.gas_201),'rs')
xlabel('cngwat1+cngwat2'); ylabel('sum(p1.gas_201)+sum(p1.gas_202)')
title('comparing CNGWATS')

disp('sarta')
eval(sartaer)
[h2,ha,p2,pa] = rtpread('xtest.rp.rtp');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% plot(h0.vchan,rad2bt(h0.vchan,p0.rcalc),'b',h0.vchan,rad2bt(h0.vchan,p2.rcalc),'r')

addpath /home/sergio/MATLABCODE/CLOUD/
g = dogoodchan;
tobs = real(rad2bt(h0.vchan,p0.robs1));
t0   = rad2bt(h0.vchan,p0.rcalc);
t2   = rad2bt(h0.vchan,p2.rcalc);

if length(p0.stemp) < 100
  plot(h0.vchan(g),tobs(g,:)-t0(g,:),'b',h0.vchan(g),tobs(g,:)-t2(g,:),'r')
  axis([650 2700 -25 +25]); grid on
end

figure(2);
dbt = -25 : 1 : +25;
n0 = hist(tobs(1291,:)-t0(1291,:),dbt);
n2 = hist(tobs(1291,:)-t2(1291,:),dbt);
  plot(dbt,n0,dbt,n2,'r');

figure(1);
dbt = 200 : 1 : 300;
nobs = hist(tobs(1291,:),dbt);
n0 = hist(t0(1291,:),dbt);
n2 = hist(t2(1291,:),dbt);
  plot(dbt,nobs,'k',dbt,n0,'b',dbt,n2,'r');
