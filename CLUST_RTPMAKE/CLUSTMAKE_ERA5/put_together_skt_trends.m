skt1 = [tm2.skt_trend_Q90_cld   tm1.skt_trend_Q90_cld   t00.skt_trend_Q90_cld   tp1.skt_trend_Q90_cld   tp2.skt_trend_Q90_cld];
  figure(1); clf; aslmap(1,rlat65,rlon73,smoothn((reshape(skt1,72,64)') ,1), [-90 +90],[-180 +180]);, title('Q90 BT1231 cld   stemp K/yr');  caxis([-1 +1]*0.15); colormap(llsmap5);
skt2 = [tm2.skt_trend_Q90_clr   tm1.skt_trend_Q90_clr   t00.skt_trend_Q90_clr   tp1.skt_trend_Q90_clr   tp2.skt_trend_Q90_clr];
  figure(2); clf; aslmap(2,rlat65,rlon73,smoothn((reshape(skt2,72,64)') ,1), [-90 +90],[-180 +180]);, title('Q90 BT1231 clr   stemp K/yr');  caxis([-1 +1]*0.15); colormap(llsmap5);
skt3 = [tm2.skt_trend_Q90       tm1.skt_trend_Q90       t00.skt_trend_Q90       tp1.skt_trend_Q90       tp2.skt_trend_Q90   ];
  figure(3); clf; aslmap(3,rlat65,rlon73,smoothn((reshape(skt3,72,64)') ,1), [-90 +90],[-180 +180]);, title('Q90 SKT          stemp K/yr');  caxis([-1 +1]*0.15); colormap(llsmap5);
skt4 = [tm2.skt_trend_mean      tm1.skt_trend_mean      t00.skt_trend_mean      tp1.skt_trend_mean      tp2.skt_trend_mean  ];
  figure(4); clf; aslmap(4,rlat65,rlon73,smoothn((reshape(skt4,72,64)') ,1), [-90 +90],[-180 +180]);, title('meanSKT          stemp K/yr');  caxis([-1 +1]*0.15); colormap(llsmap5);
skt5 = [tm2.skt_trend_cntr      tm1.skt_trend_cntr      t00.skt_trend_cntr      tp1.skt_trend_cntr      tp2.skt_trend_cntr  ];
  figure(5); clf; aslmap(5,rlat65,rlon73,smoothn((reshape(skt5,72,64)') ,1), [-90 +90],[-180 +180]);, title('cntrSKT          stemp K/yr');  caxis([-1 +1]*0.15); colormap(llsmap5);

%%%%%%%%%%%%%%%%%%%%%%%%%
  figure(1); clf; aslmap(1,rlat65,rlon73,smoothn((reshape(skt1,     72,64)') ,1), [-90 +90],[-180 +180]);, title('Q90 BT1231 cld   stemp K/yr');  caxis([-1 +1]*0.15); colormap(llsmap5);
  figure(2); clf; aslmap(2,rlat65,rlon73,smoothn((reshape(skt1-skt2,72,64)') ,1), [-90 +90],[-180 +180]);, title('Abs diff Q90 BT1231 clr   stemp K/yr');  caxis([-1 +1]*0.15/10); colormap(llsmap5);
  figure(3); clf; aslmap(3,rlat65,rlon73,smoothn((reshape(skt1-skt3,72,64)') ,1), [-90 +90],[-180 +180]);, title('Abs diff Q90 SKT          stemp K/yr');  caxis([-1 +1]*0.15/10); colormap(llsmap5);
  figure(4); clf; aslmap(4,rlat65,rlon73,smoothn((reshape(skt1-skt4,72,64)') ,1), [-90 +90],[-180 +180]);, title('Abs diff meanSKT          stemp K/yr');  caxis([-1 +1]*0.15/10); colormap(llsmap5);
  figure(5); clf; aslmap(5,rlat65,rlon73,smoothn((reshape(skt1-skt5,72,64)') ,1), [-90 +90],[-180 +180]);, title('Abs diff cntrSKT          stemp K/yr');  caxis([-1 +1]*0.15/10); colormap(llsmap5);

boo1 = squeeze(nanmean(reshape(skt1,72,64),1));
boo2 = squeeze(nanmean(reshape(skt2,72,64),1));
boo3 = squeeze(nanmean(reshape(skt3,72,64),1));
boo4 = squeeze(nanmean(reshape(skt4,72,64),1));
boo5 = squeeze(nanmean(reshape(skt5,72,64),1));

%% look at Fig 8 of trends paper submitted to JGR
figure(6); plot(meanvaluebin(rlat65),[boo1; boo2; boo3; boo4; boo5],'linewidth',2); axis([-90 +90 -0.05 +0.10])
  plotaxis2;  legend('Q90 cld BT1231','Q90 clr BT1231','Q90 skt','mean skt','center skt','location','best','fontsize',10);

%% look at Fig 8 of trends paper submitted to JGR
figure(7); plot(meanvaluebin(rlat65),[boo1; boo2; boo3; boo4; boo5]-boo1,'linewidth',2); axis([-90 +90 -0.05 +0.10])
hold on; plot(meanvaluebin(rlat65),boo1,'k.-','linewidth',4); hold off
  plotaxis2;  legend('Q90 cld BT1231','Q90 clr BT1231','Q90 skt','mean skt','center skt','location','best','fontsize',8);
  title('Absolute difference')

figure(8); plot(meanvaluebin(rlat65),100*([boo1; boo2; boo3; boo4; boo5]-boo1)./boo1,'linewidth',2); axis([-90 +90 -50 +50])
  plotaxis2;  legend('Q90 cld BT1231','Q90 clr BT1231','Q90 skt','mean skt','center skt','location','best','fontsize',10);
  title('Percent difference')
