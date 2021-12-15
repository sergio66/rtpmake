%see clust_loop_make_monthly2m_tile_center_asc.m or clust_loop_make_monthly2m_tile_center_desc.m
%!ls /asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Center/ASC/era5_tile_center_monthly2m*
%!ls /asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Center/DESC/era5_tile_center_monthly2m*

%clear all
for ii = 1 : 5; figure(ii); clf; colormap jet; end; 

%rlat64 = 1:64;

%{
addpath /asl/matlib/h4tools
[hx,ha,px,pa] = rtpread('/asl/s1/sergio/rtp/rtp_airicrad_v6/2019/04/25/cloudy_airs_l1c_ecm_sarta_baum_ice.2019.04.25.222_cumsum_-1.rtp');
scatter_coast(px.rlon,px.rlat,50,px.spres-px.plevs(91,:)); title('EMCWF 91 spres - pnew.plevs(91,:)')
%}

load /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_For_HowardObs_TimeSeries/latB64.mat
rlon = -180 : +180;          rlat = -90 : +90;
drlon = 5;                   drlat = 3;

rlon = -180 : drlon : +180;  rlat = -90 : drlat : +90;
rlon = rlon;                 rlat = latB2;

rlat64 = 0.5*(rlat(1:end-1)+rlat(2:end));

iDorA = input('desc,default (+1) or asc (-1) : ');
if length(iDorA) == 0
  iDorA = +1;
end
if iDorA > 0
  load /asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Center/DESC/era5_tile_center_monthly2m_009.mat
else
  load /asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Center/ASC/era5_tile_center_monthly2m_009.mat
end

scatter_coast(pnew_ip.rlon,pnew_ip.rlat,50,pnew_ip.spres-pnew_ip.plevs(37,:)); title('ERA5 37 spres - pnew.plevs(37,:)')
pnew_ip.spres(50)
subplot(121);
  semilogy(pnew_ip.ptemp(:,50),pnew_ip.plevs(:,50),'bo-',pnew_op.ptemp(:,50),pnew_op.plevs(:,50),'rx-'); set(gca,'ydir','reverse')
  ylim([1e-3 1000])
  ax = axis; axis([180 ax(2)+5 1e-3 1000])
  title('(b)ERA5 LEVELS','fontsize',8)
subplot(122);
  semilogy(pnew_ip.ptemp(:,50),pnew_ip.plevs(:,50),'bo-',pnew_op.ptemp(:,50),pnew_op.plevs(:,50),'rx-'); set(gca,'ydir','reverse')
  ylim([500 1000])
  ax = axis; axis([180 ax(2)+5 500 1000])
  title('(r) KLAYERS LEVELS','fontsize',8)

if ~isfield(pnew_op,'plays')
  playsN = (pnew_op.plevs(1:100,:)-pnew_op.plevs(2:101,:));
  playsD = (pnew_op.plevs(1:100,:) ./ pnew_op.plevs(2:101,:));
  pnew_op.plays = playsN ./ log(playsD);
end

[~,~,~,lowest_layer_info,RH2m_test,T2m_test] = layeramt2RH_v2(hnew_op,pnew_op);
[single(lowest_layer_info(:,1743))' pnew_op.RHSurf(1743)-pnew_op.rh2m(1743)]
figure(1); scatter_coast(pnew_ip.rlon,pnew_ip.rlat,50,lowest_layer_info(2,:)); colormap(jet); title('lowest layer frac')

figure(1); scatter_coast(pnew_ip.rlon,pnew_ip.rlat,50,T2m_test-pnew_op.t2m); colormap(usa2);   title('T2m : UMBC-ERA'); caxis([-1 +1])
figure(2); scatter_coast(pnew_ip.rlon,pnew_ip.rlat,50,RH2m_test-pnew_op.rh2m); colormap(usa2); title('RH 2m : UMBC-ERA'); caxis([-10 +10])

figure(2); plot(lowest_layer_info(2,:),pnew_op.RHSurf-pnew_op.rh2m,'.',lowest_layer_info(2,:),RH2m_test-pnew_op.rh2m,'.')
addpath /home/sergio/MATLABCODE/SHOWSTATS
[nnn,nnx,nny,nmean,nstd] = myhist2d(lowest_layer_info(2,:),pnew_op.RHSurf-pnew_op.rh2m,-0.05:0.05:1,-20:0.5:+20);
[nnn,nnx,nny,nmeanTest,nstdTest] = myhist2d(lowest_layer_info(2,:),RH2m_test-pnew_op.rh2m,-0.05:0.05:1,-20:0.5:+20);
  figure(2); errorbar(-0.05:0.05:1,nmean,nstd); hold on; errorbar((-0.05:0.05:1)+0.01,nmeanTest,nstdTest,'color','r'); hold on; 
       plot(-0.05:0.05:1,nmean,'bo-',(-0.05:0.05:1)+0.01,nmeanTest,'rx-'); hold off; xlabel('lowest layer frac'); ylabel('UMBC-ERA 2m RH'); 
  line([0 1],[0 0],'linewidth',2,'color','k'); xlim([-0.05 1.05])
%boo = find(lowest_layer_info(2,:) < 0.1); hist(pnew_op.RHSurf(boo)-pnew_op.rh2m(boo),100)
%boo = find(lowest_layer_info(2,:) > 0.1); hist(pnew_op.RHSurf(boo)-pnew_op.rh2m(boo),100)
figure(3); scatter_coast(pnew_ip.rlon,pnew_ip.rlat,50,RH2m_test-pnew_op.rh2m); caxis([-10 +10]); colormap(usa2); title('RH : mycalc-ERA5sfc = RHsurf-rh2m')

boo0 = reshape(pnew_op.RHSurf-pnew_op.rh2m,72,64); boo0 = nanmean(boo0,1);
boo1 = reshape(RH2m_test-pnew_op.rh2m,72,64);    boo1 = nanmean(boo1,1);
figure(4); clf; plot(rlat64,boo0,rlat64,boo1,'r'); xlabel('latbin'); ylabel('RHsergio-RHERA'); plotaxis2;
  hl = legend('Orig RH 2m','New RH T2m','location','best','fontsize',10);

boo1 = reshape(T2m_test-pnew_op.t2m,72,64);       boo1 = nanmean(boo1,1);
boo0 = reshape(pnew_op.stemp-pnew_op.t2m,72,64);  boo0 = nanmean(boo0,1);
figure(5); clf; plot(rlat64,boo0,rlat64,boo1,'r'); xlabel('latbin'); ylabel('T2m Xsergio-ERA'); plotaxis2;
  hl = legend('STemp','New T2m','location','best','fontsize',10);

if isfield(pnew_ip,'skt_2m')
  figure(6); clf; scatter_coast(pnew_ip.rlon,pnew_ip.rlat,50,pnew_op.stemp-pnew_op.skt_2m); caxis([-0.1 +0.1]); colormap(usa2); title('Seeing effects of daily avg on stemp')
end

%% see https://earthscience.stackexchange.com/questions/5076/how-to-calculate-specific-humidity-with-relative-humidity-temperature-and-pres
addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS/Strow_humidity/convert_humidity/
%pnew_ip.rh = convert_humidity(pnew_ip.plevs*100,pnew_ip.ptemp,pnew_ip.gas_1,'mixing ratio','relative humidity'); 
pnew_ip.rh = convert_humidity (pnew_ip.plevs*100,pnew_ip.ptemp,pnew_ip.gas_1,'specific humidity','relative humidity');
pnew_ip.rh2m = airtemp_dewpointtemp_2_RH(pnew_ip.t2m,pnew_ip.d2m);

ii = input('enter profile number to look at (1:4608, outside to stop) : ');
while ii >= 1 & ii <= 4608
  figure(6); clf; scatter_coast(pnew_ip.rlon,pnew_ip.rlat,50,pnew_op.stemp); colormap(usa2);
    hold on; plot(pnew_ip.rlon(ii),pnew_ip.rlat(ii),'kx','markersize',10,'linewidth',5); hold off

  figure(7); clf; 
     layx = pnew_op.nlevs(ii)-1;
     semilogy(pnew_ip.ptemp(:,ii),pnew_ip.plevs(:,ii),'b+-',pnew_op.ptemp(1:layx,ii),pnew_op.plays(1:layx,ii),'rx-',...
                   pnew_ip.stemp(ii),pnew_ip.spres(ii),'go',pnew_ip.t2m(ii),pnew_ip.spres(ii),'ko',T2m_test(ii),pnew_ip.spres(ii),'co',...
                   pnew_op.stemp(ii),pnew_op.spres(ii),'g.',pnew_op.t2m(ii),pnew_op.spres(ii),'k.',T2m_test(ii),pnew_ip.spres(ii),'c.',...
                   pnew_op.ptemp(1:layx+1,ii),pnew_op.plays(1:layx+1,ii),'rx--',...
                   'linewidth',2)
  set(gca,'ydir','reverse'); 
   ylim([pnew_op.spres(ii)-100 pnew_op.spres(ii)+50])
   ylim([pnew_op.spres(ii)-100 1050])

  ax = axis; 
  line([ax(1) ax(2)],[pnew_op.plevs(layx+1,ii) pnew_op.plevs(layx+1,ii)],'linewidth',2,'color','k');
  line([ax(1) ax(2)],[pnew_op.plevs(layx+0,ii) pnew_op.plevs(layx+0,ii)],'linewidth',2,'color','k');
  line([ax(1) ax(2)],[pnew_op.plevs(layx-1,ii) pnew_op.plevs(layx-1,ii)],'linewidth',2,'color','k');
  line([ax(1) ax(2)],[pnew_op.plevs(layx-2,ii) pnew_op.plevs(layx-2,ii)],'linewidth',2,'color','k');

  line([ax(1)+0.5 ax(2)-0.5],[pnew_op.plays(layx+0,ii) pnew_op.plays(layx+0,ii)],'linewidth',1,'color','r');
  line([ax(1)+0.5 ax(2)-0.5],[pnew_op.plays(layx-1,ii) pnew_op.plays(layx-1,ii)],'linewidth',1,'color','r');
  line([ax(1)+0.5 ax(2)-0.5],[pnew_op.plays(layx-2,ii) pnew_op.plays(layx-2,ii)],'linewidth',1,'color','r');  

  legend('IP','OP','go = stemp','ko = 2m','co = new guess 2m','location','best','fontsize',8)

  title('T(p)')

  %%%%%%%%%%%%%%%%%%%%%%%%%

  figure(8); clf; 
     layx = pnew_op.nlevs(ii)-1;
     semilogy(pnew_ip.rh(:,ii),pnew_ip.plevs(:,ii),'b+-',pnew_op.RH(1:layx,ii),pnew_op.plays(1:layx,ii),'rx-',...
                   pnew_ip.rh2m(ii),pnew_ip.spres(ii),'go',pnew_op.RHSurf(ii),pnew_ip.spres(ii),'ko',RH2m_test(ii),pnew_ip.spres(ii),'co',...
                   pnew_op.rh2m(ii),pnew_op.spres(ii),'g.',pnew_op.RHSurf(ii),pnew_ip.spres(ii),'k.',RH2m_test(ii),pnew_ip.spres(ii),'c.',...
                   pnew_op.RH(1:layx+1,ii),pnew_op.plays(1:layx+1,ii),'rx--',...
                   'linewidth',2)
  set(gca,'ydir','reverse'); 
   ylim([pnew_op.spres(ii)-100 pnew_op.spres(ii)+50])
   ylim([pnew_op.spres(ii)-100 1050])

  ax = axis; 
  line([ax(1) ax(2)],[pnew_op.plevs(layx+1,ii) pnew_op.plevs(layx+1,ii)],'linewidth',2,'color','k');
  line([ax(1) ax(2)],[pnew_op.plevs(layx+0,ii) pnew_op.plevs(layx+0,ii)],'linewidth',2,'color','k');
  line([ax(1) ax(2)],[pnew_op.plevs(layx-1,ii) pnew_op.plevs(layx-1,ii)],'linewidth',2,'color','k');
  line([ax(1) ax(2)],[pnew_op.plevs(layx-2,ii) pnew_op.plevs(layx-2,ii)],'linewidth',2,'color','k');

  line([ax(1)+0.5 ax(2)-0.5],[pnew_op.plays(layx+0,ii) pnew_op.plays(layx+0,ii)],'linewidth',1,'color','r');
  line([ax(1)+0.5 ax(2)-0.5],[pnew_op.plays(layx-1,ii) pnew_op.plays(layx-1,ii)],'linewidth',1,'color','r');
  line([ax(1)+0.5 ax(2)-0.5],[pnew_op.plays(layx-2,ii) pnew_op.plays(layx-2,ii)],'linewidth',1,'color','r');  

  legend('IP','OP','go = ERA 2m','ko = orig guess 2m','co = new guess 2m','location','best','fontsize',8)

  title('RH(p)')

  %%%%%%%%%%%%%%%%%%%%%%%%%

  ii = input('enter profile number to look at (1:4608, outside to stop) : ');
end

