rlatS = [];
stempS = [];
r1231S = [];
solzenS = [];

rlatH = [];
stempH = [];
r1231H = [];
solzenH = [];

for dd = 1 : 16
  rlatSX = []; rlonSX = []; stempSX = []; r1231SX = []; solzenSX = [];
  for hh = 1 : 23
    fprintf(1,'dd = %2i hh = %2i \n',dd,hh)
    fnameH = ['/asl/data/rtprod_airs/2009/09/' num2str(dd,'%02d') '/rnd_nadir6track_cloudy_airs_l1b_era_sarta_baum_ice.2009.09.'];
    fnameH = [fnameH num2str(dd,'%02d') '.' num2str(hh,'%02d') '.rtp'];

    fnameS = ['/asl/data/rtprod_airs/2009/09/' num2str(dd,'%02d') '/sergio_nadir_cloudy_airs_l1b_era_sarta_baum_ice.2009.09.'];
    fnameS = [fnameS num2str(dd,'%02d') '.' num2str(hh,'%02d') '.rtp'];

    if exist(fnameH) & exist(fnameS)
      [h,ha,p,pa] = rtpread(fnameH);
      rlatH = [rlatH p.rlat];
      stempH = [stempH p.stemp];
      r1231H = [r1231H p.robs1(1291,:)];            
      solzenH = [solzenH p.solzen];
      boo = find(p.xtrack < 43 | p.xtrack > 48);
      if length(boo) > 0
        error('ooo bad nadirs!')
      end
      
      [h,ha,p,pa] = rtpread(fnameS);
      rlatSX = [rlatSX p.rlat];
      rlonSX = [rlonSX p.rlon];
      stempSX = [stempSX p.stemp];
      r1231SX = [r1231SX p.robs1(1291,:)];
      solzenSX = [solzenSX p.solzen];      
    end
  end
  [latOUT,lonOUT,aOUT,subsample,latgridMID,nlatOUT,chosenOUT] = subsample_sphere_transcom_plot(rlatSX,rlonSX,stempSX,-90:1:+90);
  rlatS = [rlatS rlatSX(chosenOUT)];
  stempS = [stempS stempSX(chosenOUT)];
  r1231S = [r1231S r1231SX(chosenOUT)];  
  solzenS = [solzenS solzenSX(chosenOUT)];
  
  dr = -90:+90; plot(dr,hist(rlatH,dr)/sum(hist(rlatH,dr)),'b',dr,hist(rlatS,dr)/sum(hist(rlatS,dr)),'r');
  title([num2str(dd) ' ' num2str(hh)]); pause(0.1);  
end

fprintf(1,'overall means : H : stemp vs r1231 = %8.6f %8.6f \n',nanmean(stempH),nanmean(rad2bt(1231,r1231H)))
fprintf(1,'                s : stemp vs r1231 = %8.6f %8.6f \n',nanmean(stempS),nanmean(rad2bt(1231,r1231S)))

wooDH = find(solzenH < 90); wooNH = find(solzenH > 90);
wooDS = find(solzenS < 90); wooNS = find(solzenS > 90);
fprintf(1,'  day means : H : stemp vs r1231 = %8.6f %8.6f \n',nanmean(stempH(wooDH)),nanmean(rad2bt(1231,r1231H(wooDH))))
fprintf(1,'              s : stemp vs r1231 = %8.6f %8.6f \n',nanmean(stempS(wooDS)),nanmean(rad2bt(1231,r1231S(wooDS))))
fprintf(1,'  ngt means : H : stemp vs r1231 = %8.6f %8.6f \n',nanmean(stempH(wooNH)),nanmean(rad2bt(1231,r1231H(wooNH))))
fprintf(1,'              s : stemp vs r1231 = %8.6f %8.6f \n',nanmean(stempS(wooNS)),nanmean(rad2bt(1231,r1231S(wooNS))))
