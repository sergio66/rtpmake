sergio = '/asl/data/rtprod_airs/2002/09/01/sergio_nadir_cloudy_airs_l1b_era_sarta_baum_ice.2002.09.01.';
thedirS = dir([sergio '*.rtp']);

breno = '/asl/data/rtprod_airs/2002/09/01/airs_l1b.era.umw.ctrtrk.sarta_baum_ice.2002.09.01.';
thedirB = dir([breno '*.R1.9n-M1.9l.rtp']);

stempB = [];
obsB   = [];
clrB   = [];
cldB   = [];
rlonB  = [];
rlatB  = [];
landB  = [];

stempS = [];
obsS   = [];
clrS   = [];
cldS   = [];
rlonS  = [];
rlatS  = [];
landS  = [];

iaFound = [];
for ii = 00 : 23
  bfile = ['/asl/data/rtprod_airs/2002/09/01/airs_l1b.era.umw.ctrtrk.sarta_baum_ice.2002.09.01.'];
  bfile = [bfile num2str(ii,'%02d') '0000_' num2str(ii+1,'%02d') '0000.R1.9n-M1.9l.rtp'];
  iaFound(ii+1) = -1;;
  if exist(bfile)
    iaFound(ii+1) = +1;;
  
    [h,ha,p,pa] = rtpread(bfile);
    stempB = [stempB p.stemp];
    rlonB = [rlonB p.rlon];
    rlatB = [rlatB p.rlat];
    landB = [landB p.landfrac];
    obsB  = [obsB  rad2bt(1231,p.robs1(1291,:))];
    clrB  = [clrB  rad2bt(1231,p.sarta_rclearcalc(1291,:))];
    cldB  = [cldB  rad2bt(1231,p.rcalc(1291,:))];

    sfile = ['/asl/data/rtprod_airs/2002/09/01/sergio_nadir_cloudy_airs_l1b_era_sarta_baum_ice.2002.09.01.'];
    sfile = [sfile num2str(ii,'%02d') '.rtp'];
    [h,ha,p,pa] = rtpread(sfile);
    stempS = [stempS p.stemp];
    rlonS = [rlonS p.rlon];
    rlatS = [rlatS p.rlat];
    landS = [landS p.landfrac];    
    obsS  = [obsS  rad2bt(1231,p.robs1(1291,:))];
    clrS  = [clrS  rad2bt(1231,p.sarta_rclearcalc(1291,:))];
    cldS  = [cldS  rad2bt(1231,p.rcalc(1291,:))];

    ooB = find(landB == 0);
    ooS = find(landS == 0);
    xooB = find(landB == 1);
    xooS = find(landS == 1);

    figure(1); clf;
    dbt = -20 : 0.025 : +20;
    plot(dbt,hist(obsB-obsS,dbt),'k',dbt,hist(clrB-clrS,dbt),'b',dbt,hist(cldB-cldS,dbt),'r')
        hl=legend('obs diff','clr B-S','cld B-S'); set(hl,'fontsize',10)
    semilogy(dbt,hist(clrB-clrS,dbt),'b',dbt,hist(cldB-cldS,dbt),'r')
    hl=legend('clr B-S','cld B-S'); set(hl,'fontsize',10)
    
    figure(2); clf;
    dbt = -10 : 0.025 : +10;
    semilogy(dbt,hist(stempB-stempS,dbt),'k',dbt,hist(stempB(ooB)-stempS(ooS),dbt),'b',...
             dbt,hist(stempB(xooB)-stempS(xooS),dbt),'r')    
    hl=legend('all stemp B-S','ocean stemp B-S','land stemp B-S'); set(hl,'fontsize',10)
    
    figure(3); clf
    simplemap(rlatB,rlonB,stempB-stempS)
    title('stemp B-S'); 
    pause(0.1);

    figure(4); clf
    simplemap(rlatB,rlonB,stempB)
    title('stemp B');     
    pause(0.1);

  end
end
iaFound
[sum(rlonS-rlonB) sum(rlatS-rlatB)]
