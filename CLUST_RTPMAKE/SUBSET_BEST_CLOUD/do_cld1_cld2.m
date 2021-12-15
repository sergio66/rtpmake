function [tcld1,tcld2] = do_cld1_cld2(pxx,tobs899,tcal899,iBadDCCOnly);

tcld1 = ones(size(pxx.stemp)) * NaN;
tcld2 = ones(size(pxx.stemp)) * NaN;
for ii = 1 : length(pxx.stemp)
  nlevs = pxx.nlevs(ii);
  xp = pxx.plevs(1:nlevs-1,ii);
  xt = pxx.ptemp(1:nlevs-1,ii);
  if pxx.cfrac(ii) > 0 & pxx.cngwat(ii) > 0
    tcld1(ii) = interp1(log(xp),xt,log(pxx.cprtop(ii)));
  end
  if pxx.cfrac2(ii) > 0 & pxx.cngwat2(ii) > 0
    tcld2(ii) = interp1(log(xp),xt,log(pxx.cprtop2(ii)));
  end
end

theplot = 1 : length(pxx.stemp);

figure(2);
  plot(tobs899(theplot),tcld1(theplot),'bo',tobs899(theplot),tcld2(theplot),'r.',tobs899(theplot),tcal899(theplot),'k.',tobs899(theplot),tobs899(theplot),'k')
  xlabel('BT900 obs');   title('(b) tcld1 (r) tcld2 (k) tcal BT900')

figure(3);
  scatter(tobs899(theplot),tcal899(theplot),10,tcld1(theplot),'filled'); colorbar; colormap jet
  xlabel('BT900 obs');
  ylabel('BT900 cal');
  ax = axis; mna = min(ax(1),ax(3)); mxa = max(ax(2),ax(4));
  axis([mna mxa mna mxa]); ax = axis; line([ax(1) ax(2)],[ax(3) ax(4)],'linewidth',2)
  title('colorbar = tcld1')
  text(200,290,'tcld > tobs not good!')
  text(260,240,'tcld < tobs fine!')  
  
figure(4);
  scatter(tobs899(theplot),tcal899(theplot),10,tcld2(theplot),'filled'); colorbar; colormap jet
  xlabel('BT900 obs');
  ylabel('BT900 cal');
  ax = axis; mna = min(ax(1),ax(3)); mxa = max(ax(2),ax(4));
  axis([mna mxa mna mxa]); ax = axis; line([ax(1) ax(2)],[ax(3) ax(4)],'linewidth',2)
  title('colorbar = tcld2')
  text(200,290,'tcld > tobs not good!')
  text(260,240,'tcld < tobs fine!')  

figure(5)
  wah1 = pxx.cprtop;  wah1(wah1 < 0) = nan; 
  wah2 = pxx.cprtop2; wah2(wah2 < 0) = nan;
  plot(wah1(theplot)-wah2(theplot)); title('cprtop - cprtop2')

figure(6); clf
  scatter(tobs899(theplot),tcld1(theplot),10,tcal899(theplot),'filled'); colorbar; colormap jet
  xlabel('BT900 obs');
  ylabel('TCLD1');
  ax = axis; mna = min(ax(1),ax(3)); mxa = max(ax(2),ax(4));
  axis([mna mxa mna mxa]); ax = axis; line([ax(1) ax(2)],[ax(3) ax(4)],'linewidth',2)
  title('colorbar = tcalc')
  text(200,290,'tcld > tobs not good!')
  text(260,240,'tcld < tobs fine!')  

if iBadDCCOnly == 45
  theplot = find(pxx.xtrack >= 45 & pxx.xtrack <= 46);

  figure(7);
  scatter(tobs899(theplot),tcal899(theplot),10,tcld1(theplot),'filled'); colorbar; colormap jet
  xlabel('BT900 obs');
  ylabel('BT900 cal');
  ax = axis; mna = min(ax(1),ax(3)); mxa = max(ax(2),ax(4));
  axis([mna mxa mna mxa]); ax = axis; line([ax(1) ax(2)],[ax(3) ax(4)],'linewidth',2)
  title('colorbar = tcld1')
  text(200,290,'tcld > tobs not good!')
  text(260,240,'tcld < tobs fine!')  
  
  figure(8);
  scatter(tobs899(theplot),tcal899(theplot),10,tcld2(theplot),'filled'); colorbar; colormap jet
  xlabel('BT900 obs');
  ylabel('BT900 cal');
  ax = axis; mna = min(ax(1),ax(3)); mxa = max(ax(2),ax(4));
  axis([mna mxa mna mxa]); ax = axis; line([ax(1) ax(2)],[ax(3) ax(4)],'linewidth',2)
  title('colorbar = tcld2')
  text(200,290,'tcld > tobs not good!')
  text(260,240,'tcld < tobs fine!')  

  figure(9); clf
  scatter(tobs899(theplot),tcld1(theplot),10,tcal899(theplot),'filled'); colorbar; colormap jet
  xlabel('BT900 obs');
  ylabel('TCLD1');
  ax = axis; mna = min(ax(1),ax(3)); mxa = max(ax(2),ax(4));
  axis([mna mxa mna mxa]); ax = axis; line([ax(1) ax(2)],[ax(3) ax(4)],'linewidth',2)
  title('colorbar = tcalc')
  text(200,290,'tcld > tobs not good!')
  text(260,240,'tcld < tobs fine!')  
  
end

%disp('docld1cld2 ret'); pause