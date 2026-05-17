ocean = find(xpdmat.landfrac == 0);
land = find(xpdmat.landfrac > 0);

if length(ocean)/length(xpdmat.landfrac) < 0.05
  disp('oops too few ocean points, not separating')
  ocean = land;
elseif length(land)/length(xpdmat.landfrac) < 0.05
  disp('oops too few land points, not separating')
  land = ocean;
end  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1); clf
  pcolor(1:length(xpdmat.stemp),nanmean(xpdmat.plevs,2),real(log10(xpdmat.Ri))); colorbar; colormap jet; shading interp
  set(gca,'yscale','log'); set(gca,'ydir','reverse')
  title('real(log(R))'); ylabel('P [mb]')
  
figure(2); clf
  semilogy(nanmean(xpdmat.Ri(:,ocean),2),nanmean(xpdmat.plevs(:,ocean),2),'b',nanmean(xpdmat.Ri(:,land),2),nanmean(xpdmat.plevs(:,land),2),'r')
  plotaxis2;
  set(gca,'yscale','log'); set(gca,'ydir','reverse')
  xlabel('<Ri>'); ylabel('P [mb]')
  legend('ocean','land','location','best');
  axis([-4 12 500 1050])
  
figure(3); clf
  scatter_coast(xpdmat.rlon,xpdmat.rlat,10,xpdmat.zPBLH_Ri);
  title('PBLH from Ri [meters]')
  colormap jet
  
figure(4); clf
  pcolor(1:length(xpdmat.stemp),nanmean(xpdmat.plevs,2),xpdmat.stable); colorbar; colormap jet; shading interp
  set(gca,'yscale','log'); set(gca,'ydir','reverse')
  caxis([-2 +2]); colormap(usa2)
  title('stability +1=US, 0=CS, -1=AS, -2=N')
  hold on; plot(1:length(xpdmat.stemp),xpdmat.pPBLH_Ri,'k','linewidth',2); hold off
  hold on; plot(1:length(xpdmat.stemp),xpdmat.spres,'g','linewidth',2); hold off  
  ylim([10 1050])

DALR = 10;  %% dry adiatic lapse rate
MALR = 6;   %% moist adiatic lapse rate
figure(5); clf
  semilogy(nanmean(xpdmat.lapse(:,ocean),2),nanmean(xpdmat.plevs(:,ocean),2),'b',nanmean(xpdmat.lapse(:,land),2),nanmean(xpdmat.plevs(:,land),2),'r')
     ylim([10 1000]); set(gca,'ydir','reverse');   xlim([-4 12])
  plotaxis2; xlabel('lapse rate [K/km]'); ylabel('P [mb]');
  line([MALR MALR],[0.005 1100],'color','b'); text(4,50,'Moist LR','color','b');   text(4.00,75,'Stable','color','b');    text(4.00,25, 'S = -1','color','b'); 
  line([DALR DALR],[0.005 1100],'color','r'); text(10.25,50,'Dry LR','color','r'); text(10.25,75,'Unstable','color','r'); text(10.25,25,'S = +1','color','r'); 
  text(6.25,50,'Conditionally','color','g');     text(6.25,75,'Stable','color','g');   text(6.25,25,'S = 0','color','g'); 
  title('Lapse Rate K/km')
  legend('ocean','land','location','best')
  
figure(6); clf
  semilogy(nanmean(xpdmat.Ri(:,ocean),2),nanmean(xpdmat.plevs(:,ocean),2),'b',nanmean(xpdmat.Ri(:,land),2),nanmean(xpdmat.plevs(:,land),2),'r')
  axis([-5 50 500 1000])  
  plotaxis2;
  plotaxis2(0.27,0,'r');       
  title('Bulk Richardson')
  legend('ocean','land','location','best')
  set(gca,'ydir','reverse');   xlim([-4 12])
  
figure(7); clf
  scatter_coast(xpdmat.rlon,xpdmat.rlat,10,rad2bt(1231,xpdmat.robs1))
  title('BT 1231 obs')
  colormap jet

figure(8); clf
  [mm,nn] = size(xpdmat.lapse);
  yyaxis left;  semilogy(nanmean(xpdmat.lapse(:,ocean),2),nanmean(xpdmat.plevs(:,ocean),2),'b');
    ylim([10 1000]); set(gca,'ydir','reverse'); axis([-4 12 500 1050]); set(gca, 'XAxisLocation', 'bottom'); xlabel('Lapse rate K/km')
  yyaxis right; semilogy(nanmean(xpdmat.Ri(:,ocean),2),nanmean(xpdmat.plevs(:,ocean),2),'r');
    ylim([10 1000]); set(gca,'ydir','reverse');    axis([-4 12 500 1050]); set(gca, 'XAxisLocation', 'top');    xlabel('Bulk Ri number'); 
    plotaxis2;
    plotaxis2(0.27,0,'r');     
  legend('Ocean Lapse rate K/km','Ocean Bulk Richardon number','location','best')

figure(8); clf
  [mm,nn] = size(xpdmat.lapse);
  ax1 = gca;
  ax1.Position = [0.15, 0.15, 0.75, 0.75]; 
  
  yyaxis left;  semilogy(nanmean(xpdmat.lapse(:,ocean),2),nanmean(xpdmat.plevs(:,ocean),2),'b'); set(gca,'ydir','reverse'); axis([-4 12 500 1050]);
    xlabel('Ocean Lapse rate K/km','color','blue')
  yyaxis right; semilogy(nanmean(xpdmat.Ri(:,ocean),2),nanmean(xpdmat.plevs(:,ocean),2),'r',0.27*ones(1,mm),nanmean(xpdmat.plevs,2),'r');
    set(gca,'ydir','reverse');    axis([-4 12 500 1050]);

  % 2. Create the secondary top axes
  ax2 = axes('Position', ax1.Position, ...
           'XAxisLocation', 'top', ...
           'YAxisLocation', 'right', ... % Match right side to avoid extra lines
           'Color', 'none', ...          % Make it transparent
           'YTick', []);                 % Hide duplicate Y ticks

  % Add the top label
  xlabel(ax2, 'Ocean Bulk Ri number','color','red'); 

  plotaxis2;
  line([0.27 0.27],log([500 1050]),'color','r')

  % 3. Synchronize both X-axes
  linkaxes([ax1, ax2], 'x');
  legend('(b) Lapse rate K/km','(r) Bulk Richardon number','location','best')
