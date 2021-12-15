function [p, h, pattr, pClosest] = fill_era_interp(p, h, pattr)

%% based on /home/sergio/MATLABCODE/RTPMAKE/CLUST_RTPMAKE/CLUSTMAKE_ECM/cloud_set_defaults_run_maker_interp_analysis.m

addpath /home/sergio/MATLABCODE/TIME
addpath /home/sergio/MATLABCODE/PLOTTER
addpath /home/sergio/MATLABCODE/COLORMAP
addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS
addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS/IDL_WV_ROUTINES/atmos_phys/MATLAB/
addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS/Strow_humidity/convert_humidity/

fprintf(1,'mean p.rtime = %8.6e ',mean(p.rtime));
[xyyh,xmmh,xddh,xhhh] = tai2utcSergio(mean(p.rtime));
fprintf(1,' corresponds to following : %4i %2i %2i %6.3f \n',xyyh,xmmh,xddh,xhhh);
%fprintf(1,'when hhh = %6.3f is same as %2i : %2i \n',hhh,floor(xhhh),floor((xhhh-floor(xhhh))*60));

frac1 = ones(size(p.rtime));
frac2 = zeros(size(p.rtime));
timeslot = zeros(size(p.rtime));

if nargin == 2
  pattr = [];
end

[xyyh,xmmh,xddh,xhhh] = tai2utcSergio(p.rtime);
oo_0_6 = find(xhhh >= 0 & xhhh <= 6);
if length(oo_0_6) > 0
  tB1_0_6 = utc2taiSergio(xyyh,xmmh,xddh,0.0);
  tB2_0_6 = utc2taiSergio(xyyh,xmmh,xddh,6.0);
  frac1(oo_0_6) = 1 - (p.rtime(oo_0_6)-tB1_0_6(oo_0_6))/(6*60*60);
  donk = find(frac1(oo_0_6) < 0); frac1(oo_0_6(donk)) = 0;    
  donk = find(frac1(oo_0_6) > 1); frac1(oo_0_6(donk)) = 1;    
  frac2(oo_0_6) = 1-frac1(oo_0_6);
  timeslot(oo_0_6) = 1;
end
oo_6_12 = find(xhhh > 6 & xhhh <= 12);
if length(oo_6_12) > 0
  tB1_6_12 = utc2taiSergio(xyyh,xmmh,xddh,6.0);
  tB2_6_12 = utc2taiSergio(xyyh,xmmh,xddh,12.0);      
  frac1(oo_6_12) = 1 - (p.rtime(oo_6_12)-tB1_6_12(oo_6_12))/(6*60*60);
  donk = find(frac1(oo_6_12) < 0); frac1(oo_6_12(donk)) = 0;    
  donk = find(frac1(oo_6_12) > 1); frac1(oo_6_12(donk)) = 1;    
  frac2(oo_6_12) = 1-frac1(oo_6_12);
  timeslot(oo_6_12) = 2;
end
oo_12_18 = find(xhhh > 12 & xhhh <= 18);
if length(oo_12_18) > 0
  tB1_12_18 = utc2taiSergio(xyyh,xmmh,xddh,12.0);
  tB2_12_18 = utc2taiSergio(xyyh,xmmh,xddh,18.0);      
  frac1(oo_12_18) = 1 - (p.rtime(oo_12_18)-tB1_12_18(oo_12_18))/(6*60*60);
  donk = find(frac1(oo_12_18) < 0); frac1(oo_12_18(donk)) = 0;    
  donk = find(frac1(oo_12_18) > 1); frac1(oo_12_18(donk)) = 1;    
  frac2(oo_12_18) = 1-frac1(oo_12_18);
  timeslot(oo_12_18) = 3;
end
oo_18_24 = find(xhhh > 18 & xhhh <= 24);
if length(oo_18_24) > 0
  tB1_18_24 = utc2taiSergio(xyyh,xmmh,xddh,18.0);
  tB2_18_24 = utc2taiSergio(xyyh,xmmh,xddh,24.0);
  frac1(oo_18_24) = 1 - (p.rtime(oo_18_24)-tB1_18_24(oo_18_24))/(6*60*60);
  donk = find(frac1(oo_18_24) < 0); frac1(oo_18_24(donk)) = 0;    
  donk = find(frac1(oo_18_24) > 1); frac1(oo_18_24(donk)) = 1;    
  frac2(oo_18_24) = 1-frac1(oo_18_24);
  timeslot(oo_18_24) = 4;
end

if length(oo_0_6)+length(oo_6_12)+length(oo_12_18)+length(oo_18_24) ~= length(p.rtime)
  error('oops something wrong in partitioning the times')
end
fprintf(1,'partitioning into 0-6-12-18-24 = %4i %4i %4i %4i \n',length(oo_0_6),length(oo_6_12),length(oo_12_18),length(oo_18_24))

figure(7); plot(oo_0_6,frac1(oo_0_6),'b.',oo_6_12,frac1(oo_6_12),'b.',oo_12_18,frac1(oo_12_18),'b.',oo_18_24,frac1(oo_18_24),'k.')
ylim([0 1])
pause(0.1)

disp('doing basic usual')
  [pClosest,hClosest] = fill_era(p,h);
  p = pClosest;
  h = hClosest;

disp('now subbing in the time intervals as needed')
if length(oo_0_6) > 0
  p = fracB1B2(h,p,frac1,frac2,oo_0_6,tB1_0_6,tB2_0_6);        p1 = p;
end
if length(oo_6_12) > 0
  p = fracB1B2(h,p,frac1,frac2,oo_6_12,tB1_6_12,tB2_6_12);     p2 = p;
end
if length(oo_12_18) > 0
  p = fracB1B2(h,p,frac1,frac2,oo_12_18,tB1_12_18,tB2_12_18);  p3 = p;
end
if length(oo_18_24) > 0
  p = fracB1B2(h,p,frac1,frac2,oo_18_24,tB1_18_24,tB2_18_24);  p4 = p;
end

figure(8); scatter_coast(pClosest.rlon+rand(size(p.stemp))*0.1,pClosest.rlat+rand(size(p.stemp))*0.1,30,p.stemp-pClosest.stemp);
if length(p.rlon) > 10
  xlim([min(p.rlon)-10 max(p.rlon)-10]);  ylim([min(p.rlat)-10 max(p.rlat)-10]);
end
figure(9); scatter(xhhh,p.stemp-pClosest.stemp,10,p.landfrac); colorbar; xlabel('xhhh (UTC)'); ylabel('\delta(ST)'); title('colorbar = landfrac')

p321 = find(p.plevs(:,1) > 321,1);

figure(1); scatter_coast(pClosest.rlon,pClosest.rlat,30,p.gas_1(p321,:));
figure(2); scatter_coast(pClosest.rlon,pClosest.rlat,30,pClosest.gas_1(p321,:));
%figure(3); scatter_coast(pClosest.rlon,pClosest.rlat,30,pB1.gas_1(p321,:));
%figure(4); scatter_coast(pClosest.rlon,pClosest.rlat,30,pB2.gas_1(p321,:));
%figure(5); semilogy(pClosest.ptemp,pClosest.plevs,'k',pB1.ptemp,pB1.plevs,'c',pB2.ptemp,pB2.plevs,'b',p.ptemp,p.plevs,'r');
%  set(gca,'ydir','reverse')
%figure(6); loglog(pClosest.gas_1,pClosest.plevs,'k',pB1.gas_1,pB1.plevs,'c',pB2.gas_1,pB2.plevs,'b',p.gas_1,p.plevs,'r');
%  set(gca,'ydir','reverse')

rh321      = mixr2rh(p.gas_1(p321,:)*1000,p.plevs(p321,:),p.ptemp(p321,:),1);  %% from IDL routines
irionrh321 = mixr2rh(p.gas_1(p321,:)*1000,p.plevs(p321,:),p.ptemp(p321,:),0);  %% from IDL routines, with SVP a linear mix of I/W (Irion)
rah = p.gas_1(p321,:); %% SH in g/g
rah = rah./(1-rah);       %% mixR in g/g
irionrh321 = mixr2rh(rah*1000,p.plevs(p321,:),p.ptemp(p321,:),0);  %% from IDL routines, with SVP a linear mix of I/W (Irion)
xrh321 = convert_humidity(p.plevs(p321,:),p.ptemp(p321,:),p.gas_1(p321,:),'mixing ratio','relative humidity');  %% from Strow
xrh321 = xrh321*100;
figure(1); clf; scatter_coast(pClosest.rlon,pClosest.rlat,30,irionrh321); caxis([0 120]); colormap jet; title('analysis interp RH')
colormap(irion_idl_colormap);
    
rh321      = mixr2rh(pClosest.gas_1(p321,:)*1000,pClosest.plevs(p321,:),pClosest.ptemp(p321,:),1);  %% from IDL routines
irionrh321 = mixr2rh(pClosest.gas_1(p321,:)*1000,pClosest.plevs(p321,:),pClosest.ptemp(p321,:),0);  %% from IDL routines, with SVP a linear mix of I/W (Irion)
rah = pClosest.gas_1(p321,:); %% SH in g/g
rah = rah./(1-rah);       %% mixR in g/g
irionrh321 = mixr2rh(rah*1000,pClosest.plevs(p321,:),pClosest.ptemp(p321,:),0);  %% from IDL routines, with SVP a linear mix of I/W (Irion)
xrh321 = convert_humidity(pClosest.plevs(p321,:),pClosest.ptemp(p321,:),pClosest.gas_1(p321,:),'mixing ratio','relative humidity');  %% from Strow
xrh321 = xrh321*100;
figure(2); clf; scatter_coast(pClosest.rlon,pClosest.rlat,30,irionrh321); caxis([0 120]); colormap jet; title('pClosest in time RH')
colormap(irion_idl_colormap);

figure(3); clf; scatter_coast(pClosest.rlon,pClosest.rlat,30,xrh321-irionrh321); colormap jet; title('difference RH :  closest-intep')

for ip = 1 : 3
  figure(ip); colormap jet
  colormap(irion_idl_colormap);  
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
