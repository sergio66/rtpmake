%for ii=1:13; figure(ii); colormap jet; end

qqq0 = (1:length(quants));
qqm1 = meanvaluebin(1:length(quants));

qqq0 = quants(1:length(quants));
qqm1 = meanvaluebin(quants(1:length(quants)));

figure(1); clf
  plot(qqq0,pav.stemp,'bx',qqq0,globalavg.stemp,'ro')
  plot(qqm1,pQav.stemp,'gx',qqm1,globalQavg.stemp,'mo')
  plot(qqq0,pav.stemp,'bx',qqq0,globalavg.stemp,'ro',qqm1,pQav.stemp,'gx',qqm1,globalQavg.stemp,'mo','linewidth',2,'markersize',5)

figure(1); dbt = 200:1:320; semilogy(dbt,hist(p2.stemp,dbt),'b',dbt,hist(rad2bt(1231,p2.rcalc(1520,:)),dbt),'r','linewidth',2); title('hist(stemp,BT1231 cld)')
figure(1); dbt = 200:1:320; semilogy(dbt,hist(p2.stemp,dbt),'b',dbt,hist(rad2bt(1231,p2.rcalc(1520,:)),dbt),'r',dbt,hist(rad2bt(1231,p2.sarta_rclearcalc(1520,:)),dbt),'m','linewidth',2); title('hist(stemp,BT1231 cld,clr')
    hl = legend('stemp','BT1231 cld','BT1231 clr','location','best');

figure(2); clf
  plot(qqq0,globalavg.cfrac,'bx-',qqm1,globalQavg.cfrac,'rx-','linewidth',2,'markersize',5)
    xlabel('Quantile'); ylabel('ice cldfrac'); hl = legend('\int Q(x-->1)','\int Q(x)-->Q(x+1)','location','best');
figure(3); clf
  plot(qqq0,globalavg.iceOD,'bx-',qqm1,globalQavg.iceOD,'rx-','linewidth',2,'markersize',5)
    xlabel('Quantile'); ylabel('ice OD'); hl = legend('\int Q(x-->1)','\int Q(x)-->Q(x+1)','location','best');
figure(4); clf
  plot(qqq0,globalavg.icetop,'bx-',qqm1,globalQavg.icetop,'rx-','linewidth',2,'markersize',5)
    xlabel('Quantile'); ylabel('ice cldtop'); hl = legend('\int Q(x-->1)','\int Q(x)-->Q(x+1)','location','best');
    set(gca,'ydir','reverse')

figure(5); clf
  plot(qqq0,globalavg.cfrac2,'bx-',qqm1,globalQavg.cfrac2,'rx-','linewidth',2,'markersize',5)
    xlabel('Quantile'); ylabel('water cldfrac'); hl = legend('\int Q(x-->1)','\int Q(x)-->Q(x+1)','location','best');
figure(6); clf
  plot(qqq0,globalavg.waterOD,'bx-',qqm1,globalQavg.waterOD,'rx-','linewidth',2,'markersize',5)
    xlabel('Quantile'); ylabel('water OD'); hl = legend('\int Q(x-->1)','\int Q(x)-->Q(x+1)','location','best');
figure(7); clf
  plot(qqq0,globalavg.watertop,'bx-',qqm1,globalQavg.watertop,'rx-','linewidth',2,'markersize',5)
    xlabel('Quantile'); ylabel('water cldtop'); hl = legend('\int Q(x-->1)','\int Q(x)-->Q(x+1)','location','best');
    set(gca,'ydir','reverse')

%%%%%%%%%%%%%%%%%%%%%%%%%
figure(8); clf; colormap jet

subplot(211)
hold on
plot(qqq0,pav.stemp,'cx-','markersize',10,'linewidth',2)
plot(qqm1,pQav.stemp,'mx-','markersize',10,'linewidth',2)
scatter(qqq0,pav.stemp,50,quants,'filled'); colorbar; 
scatter(qqm1,pQav.stemp,150,quants(1:end-1),'filled');
hold off
xlabel('Quantile'); ylabel('stemp'); title('colorbar = quants X-1') 
hl = legend('\int Q(x-->1)','\int Q(x)-->Q(x+1)','location','best');

subplot(212)
hold on
plot(qqq0,rad2bt(1231,pav.cldcalc),'cx-','markersize',10,'linewidth',2)
plot(qqm1,rad2bt(1231,pQav.cldcalc),'mx-','markersize',10,'linewidth',2)
scatter(qqq0,rad2bt(1231,pav.cldcalc),50,quants,'filled'); colorbar; 
scatter(qqm1,rad2bt(1231,pQav.cldcalc),150,quants(1:end-1),'filled');
hold off
xlabel('Quantile'); ylabel('cldcalc'); title('colorbar = quants X-1') 
hl = legend('\int Q(x-->1)','\int Q(x)-->Q(x+1)','location','best');

%%%%%%%%%%%%%%%%%%%%%%%%%
figure(9)
plot(globalavg.ptemp,1:101); xlim([180 300]); set(gca,'ydir','reverse')
  plot(globalavg.ptemp-globalavg.ptemp(:,length(quants)),1:101); set(gca,'ydir','reverse')
semilogy(globalavg.ptemp,globalavg.plevs); xlim([180 300]); set(gca,'ydir','reverse')
  semilogy(globalavg.ptemp-globalavg.ptemp(:,length(quants)),globalavg.plevs); set(gca,'ydir','reverse'); ylim([10 1000])

figure(10)
semilogx(globalavg.gas_1,1:101); set(gca,'ydir','reverse')
  plot(globalavg.gas_1 ./ globalavg.gas_1(:,length(quants)),1:101); set(gca,'ydir','reverse')
loglog(globalavg.gas_1,globalavg.plevs); set(gca,'ydir','reverse')
  semilogy(globalavg.gas_1 ./ globalavg.gas_1(:,length(quants)),globalavg.plevs); set(gca,'ydir','reverse'); ylim([10 1000])

%%%%%

figure(9)
pcolor(qqq0,globalavg.plevs(1:97,:),globalavg.ptemp(1:97,:));  set(gca,'ydir','reverse'); ylim([10 1000]); shading flat; colorbar;
pcolor(qqq0,globalavg.plevs(1:97,:),globalavg.ptemp(1:97,:)-globalavg.ptemp(1:97,length(quants)));  set(gca,'ydir','reverse'); ylim([10 1000]); shading flat; colorbar;
  caxis([-1 +1]*5)
title('\int Q(x)->Q(1) T - T(hottest BT1231)')
colormap(usa2);
shading interp;

figure(10)
pcolor(qqq0,globalavg.plevs(1:97,:),log10(globalavg.gas_1(1:97,:)));  set(gca,'ydir','reverse'); ylim([10 1000]); shading flat; colorbar;
pcolor(qqq0,globalavg.plevs(1:97,:),globalavg.gas_1(1:97,:) ./ globalavg.gas_1(1:97,length(quants)));  set(gca,'ydir','reverse'); ylim([10 1000]); shading flat; colorbar;
title('\int Q(x)->Q(1) Q / Q(hottest BT1231)')
colormap(jet); 
shading interp;

%%%%%%%%%%%%%%%%%%%%%%%%%

iPlotOther = 1;
%iPlotOther = input('plot raw old quantiles? (+1 default)');
if length(iPlotOther) == 0
  iPlotOther = 1;
end

if iPlotOther > 0
  figure(11)
  pcolor(qqm1,globalQavg.plevs(1:97,:),globalQavg.ptemp(1:97,:));  set(gca,'ydir','reverse'); ylim([10 1000]); shading flat; colorbar;
  pcolor(qqm1,globalQavg.plevs(1:97,:),globalQavg.ptemp(1:97,:)-globalQavg.ptemp(1:97,length(quants)-1));  set(gca,'ydir','reverse'); ylim([10 1000]); shading flat; colorbar;
    caxis([-1 +1]*5)
  title('Q(x)->Q(x+1) T - T(hottest BT1231)')
  colormap(usa2);  
  shading interp;

  figure(12)
  pcolor(qqm1,globalQavg.plevs(1:97,:),log10(globalQavg.gas_1(1:97,:)));  set(gca,'ydir','reverse'); ylim([10 1000]); shading flat; colorbar;
  pcolor(qqm1,globalQavg.plevs(1:97,:),globalQavg.gas_1(1:97,:) ./ globalQavg.gas_1(1:97,length(quants)-1));  
    set(gca,'ydir','reverse'); ylim([10 1000]); shading flat; colorbar;
  title('Q(x)->Q(x+1) WV / WV(hottest BT1231)')
  colormap(jet); 
  shading interp;
end
  
figure(13);
scatter_coast(p2.rlon,p2.rlat,50,bt1231); title('BT1231 calculated')
