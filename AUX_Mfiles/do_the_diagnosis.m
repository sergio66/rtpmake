
n0 = hist(tobs(oo), bt);
nP = hist(tcldP(oo),bt);
nS  = hist(tcldS(oo),bt); 
xnS = hist(xtcldS(oo),bt); 
ynS = hist(ytcldS(oo),bt); 
figure(1); plot(bt,[n0; nP; nS; xnS; ynS]); hl = legend('OBS','PCRTM','SARTA','SARTA 2'); set(hl,'fontsize',10)
figure(1); plot(bt,n0,'k','linewidth',3); hold on
           plot(bt,[nP; nS; xnS; ynS]); hold off
hl = legend('OBS','PCRTM','SARTA','SARTA 2','SARTA 3'); set(hl,'fontsize',10)
[nanmean(tobs(oo)-tcldP(oo)) nanstd(tobs(oo)-tcldP(oo)) ...
 nanmean(tobs(oo)-tcldS(oo)) nanstd(tobs(oo)-tcldS(oo)) ...
 nanmean(tobs(oo)-xtcldS(oo)) nanstd(tobs(oo)-xtcldS(oo)) ...
 nanmean(tobs(oo)-ytcldS(oo)) nanstd(tobs(oo)-ytcldS(oo))]

figure(2); plot(bt,n0-nP,'b',bt,n0-nS,'g',bt,n0-xnS,'r',bt,n0-ynS,'c')
hl = legend('PCRTM','SARTA','SARTA 2','SARTA 3'); set(hl,'fontsize',10); grid

figure(2); plot(bt,abs(n0-nP),'b',bt,abs(n0-nS),'g',bt,abs(n0-xnS),'r',bt,abs(n0-ynS),'c')
hl = legend('PCRTM','SARTA','SARTA 2','SARTA 3'); set(hl,'fontsize',10); grid

figure(2)
dd = 50 : 50 : 1000; 
  nS1 = hist(cprtopS1(oo),dd); xnS1 = hist(xcprtopS1(oo),dd); nP1 = hist(cprtopP1(oo),dd); plot(dd,nS1,dd,nP1,dd,xnS1)
  nS2 = hist(cprtopS2(oo),dd); xnS2 = hist(xcprtopS2(oo),dd); nP2 = hist(cprtopP2(oo),dd); plot(dd,nS2,dd,nP2,dd,xnS2)
  %nS1(1) = NaN; nP1(1) = NaN;   nS2(1) = NaN; nP2(1) = NaN;
  plot(dd,nS1,'b',dd,nS2,'r',dd,xnS1,'b--',dd,xnS2,'r--',dd,nP1,'c',dd,nP2,'m')
  hl = legend('S ice','S water','S2 ice','S2 water','P ice','P water','location','northwest');
  set(hl,'fontsize',10)

figure(3)
plot(cprtopP1(oo),tobs(oo)-tcldP(oo),'.',cprtopS1(oo)+25,tobs(oo)-tcldS(oo),'ro',...
                                         xcprtopS1(oo)+25,tobs(oo)-xtcldS(oo),'kx')
  xlabel('ICE TOP (mb)'); ylabel('Obs-Cal'); title('ICE (b) PCRTM (r) SARTA (k) SARTA 2')
  axis([0 1000 -60 +60])

figure(4)
plot(cprtopP2(oo),tobs(oo)-tcldP(oo),'.',cprtopS2(oo)+25,tobs(oo)-tcldS(oo),'ro',...
                                         xcprtopS2(oo)+25,tobs(oo)-xtcldS(oo),'kx')
  xlabel('WATER TOP (mb)'); ylabel('Obs-Cal'); title('WATER (b) PCRTM (r) SARTA (k) SARTA 2')
  axis([0 1000 -60 +60])

biasBins = -60:2.5:+60;
biasBins = -30:2.5:+30;

%{
% [n x y nmean nstd] = myhist2d(fx,fy,limitsx, limitsy,LinearOrLog)
[xnwS xjunk yjunk] = myhist2d(xcprtopS2(oo),tobs(oo)-xtcldS(oo),50 : 50 : 1000,biasBins);
[nwS xjunk yjunk] = myhist2d(cprtopS2(oo),tobs(oo)-tcldS(oo),50 : 50 : 1000,biasBins);
[nwP xjunk yjunk]= myhist2d(cprtopP2(oo),tobs(oo)-tcldP(oo),50 : 50 : 1000,biasBins);

[xniS xjunk yjunk] = myhist2d(xcprtopS1(oo),tobs(oo)-xtcldS(oo),50 : 50 : 1000,biasBins);
[niS xjunk yjunk] = myhist2d(cprtopS1(oo),tobs(oo)-tcldS(oo),50 : 50 : 1000,biasBins);
[niP xjunk yjunk]= myhist2d(cprtopP1(oo),tobs(oo)-tcldP(oo),50 : 50 : 1000,biasBins);

figure(4); pcolor(xjunk,yjunk,log10(niP)); caxis([0 4]); colorbar;
  title('ICE PCRTM'); xlabel('ICE TOP (mb)'); ylabel('Obs-Cal'); 
figure(5); pcolor(xjunk,yjunk,log10(niS)); caxis([0 4]); colorbar; 
  title('ICE SARTA'); xlabel('ICE TOP (mb)'); ylabel('Obs-Cal');
figure(6); pcolor(xjunk,yjunk,log10(xniS)); caxis([0 4]); colorbar; 
  title('ICE SARTA 2'); xlabel('ICE TOP (mb)'); ylabel('Obs-Cal');

figure(7); pcolor(xjunk,yjunk,log10(nwP)); caxis([0 4]); colorbar; 
  title('WATER PCRTM'); xlabel('WATER TOP (mb)'); ylabel('Obs-Cal');
figure(8); pcolor(xjunk,yjunk,log10(nwS)); caxis([0 4]); colorbar; 
  title('WATER SARTA'); xlabel('WATER TOP (mb)'); ylabel('Obs-Cal');
figure(9); pcolor(xjunk,yjunk,log10(xnwS)); caxis([0 4]); colorbar; 
  title('WATER SARTA 2'); xlabel('WATER TOP (mb)'); ylabel('Obs-Cal');
%}

[niP xjunk yjunk]  = myhist2d_plot(cprtopP1(oo),tobs(oo)-tcldP(oo),50 : 50 : 1000, biasBins, -1,4,[0 4]);
  title('ICE PCRTM'); xlabel('ICE TOP (mb)'); ylabel('Obs-Cal');
[niS xjunk yjunk]  = myhist2d_plot(cprtopS1(oo),tobs(oo)-tcldS(oo),50 : 50 : 1000, biasBins, -1,5,[0 4]);
  title('ICE SARTA'); xlabel('ICE TOP (mb)'); ylabel('Obs-Cal');
[xniS xjunk yjunk] = myhist2d_plot(xcprtopS1(oo),tobs(oo)-xtcldS(oo),50 : 50 : 1000, biasBins, -1,6,[0 4]);
  title('ICE SARTA 2'); xlabel('ICE TOP (mb)'); ylabel('Obs-Cal');

[nwP xjunk yjunk]  = myhist2d_plot(cprtopP2(oo),tobs(oo)-tcldP(oo),50 : 50 : 1000, biasBins, -1,7,[0 4]);
  title('WATER PCRTM'); xlabel('WATER TOP (mb)'); ylabel('Obs-Cal');
[nwS xjunk yjunk]  = myhist2d_plot(cprtopS2(oo),tobs(oo)-tcldS(oo),50 : 50 : 1000, biasBins, -1,8,[0 4]);
  title('WATER SARTA'); xlabel('WATER TOP (mb)'); ylabel('Obs-Cal');
[xnwS xjunk yjunk] = myhist2d_plot(xcprtopS2(oo),tobs(oo)-xtcldS(oo),50 : 50 : 1000, biasBins, -1,9,[0 4]);
  title('WATER SARTA 2'); xlabel('WATER TOP (mb)'); ylabel('Obs-Cal');

[sum(nwS(:)) sum(nwP(:)) sum(niS(:)) sum(niP(:)) sum(xnwS(:)) sum(xniS(:))]
ret
%%%%%%%%%%%%%%%%%%%%%%%%%

[niP xjunk yjunk]  = myhist2d_sergio_plot(cprtopP1(oo),tobs(oo)-tcldP(oo),50 : 50 : 1000, biasBins, -1,4,[0 4]);
  title('ICE PCRTM'); xlabel('ICE TOP (mb)'); ylabel('Obs-Cal');
[niS xjunk yjunk]  = myhist2d_sergio_plot(cprtopS1(oo),tobs(oo)-tcldS(oo),50 : 50 : 1000, biasBins, -1,5,[0 4]);
  title('ICE SARTA'); xlabel('ICE TOP (mb)'); ylabel('Obs-Cal');
[xniS xjunk yjunk] = myhist2d_sergio_plot(xcprtopS1(oo),tobs(oo)-xtcldS(oo),50 : 50 : 1000, biasBins, -1,6,[0 4]);
  title('ICE SARTA 2'); xlabel('ICE TOP (mb)'); ylabel('Obs-Cal');

[nwP xjunk yjunk]  = myhist2d_sergio_plot(cprtopP2(oo),tobs(oo)-tcldP(oo),50 : 50 : 1000, biasBins, -1,7,[0 4]);
  title('WATER PCRTM'); xlabel('WATER TOP (mb)'); ylabel('Obs-Cal');
[nwS xjunk yjunk]  = myhist2d_sergio_plot(cprtopS2(oo),tobs(oo)-tcldS(oo),50 : 50 : 1000, biasBins, -1,8,[0 4]);
  title('WATER SARTA'); xlabel('WATER TOP (mb)'); ylabel('Obs-Cal');
[xnwS xjunk yjunk] = myhist2d_sergio_plot(xcprtopS2(oo),tobs(oo)-xtcldS(oo),50 : 50 : 1000, biasBins, -1,9,[0 4]);
  title('WATER SARTA 2'); xlabel('WATER TOP (mb)'); ylabel('Obs-Cal');

[sum(nwS(:)) sum(nwP(:)) sum(niS(:)) sum(niP(:)) sum(xnwS(:)) sum(xniS(:))]



