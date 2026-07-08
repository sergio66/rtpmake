disp('doing trends of wspeed and stemp for 4608 tiles : progress    + is 1000     . is 100')

if ~exist('p0')
  use_this_rtp = 'RTP/summary_19years_all_lat_all_lon_2002_2021_monthlyERA5.rp.rtp';
  [h0,ha,p00,pa] = rtpread(['/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/' use_this_rtp]);
  p0 = p00;
end

for iCnt = 1 : iCntMax
  daysSince2002(iCnt) = change2days(yysave(iCnt),mmsave(iCnt),15,2002);
end

for ii = 1 : 4608
  if mod(ii,1000) == 0
    fprintf(1,'+')
  elseif mod(ii,100) == 0
    fprintf(1,'.')
  end
  [B, err, stats] = Math_tsfit_lin_robust(daysSince2002,wspeedsave(:,ii),4);
  trend.wspeed(ii) = B(2);
  trend_err.wspeed(ii) = err.se(2);

  [B, err, stats] = Math_tsfit_lin_robust(daysSince2002,stempsave(:,ii),4);
  trend.stemp(ii) = B(2);
  trend_err.stemp(ii) = err.se(2);

  [B, err, stats] = Math_tsfit_lin_robust(daysSince2002,tccsave(:,ii),4);
  trend.tcc(ii) = B(2);
  trend_err.tcc(ii) = err.se(2);
end
fprintf(1,'\n')

addpath /home/sergio/MATLABCODE/COLORMAP/COLORBREWER/cbrewer/cbrewer
blues = flipud(cbrewer('seq', 'Blues', 256));

[p0.salti,p0.landfrac] = usgs_deg10_dem(p0.rlat, p0.rlon);
ocean = find(p0.landfrac == 0); whos ocean
plot(-0.1:0.01:+0.1,histc(trend.wspeed(ocean),-0.1:0.01:+0.1)); grid

figure(1); clf; scatter_coast(p0.rlon,p0.rlat,50,nanmean(wspeedsave,1)); caxis([0 1]*10);  colormap(jet); title('<speed> over 20 years');
figure(2); clf; scatter_coast(p0.rlon,p0.rlat,50,nanmean(stempsave,1));  caxis([200 300]); colormap(jet); title('<stemp> over 20 years')
figure(3); clf; scatter_coast(p0.rlon,p0.rlat,50,nanmean(tccsave,1));    caxis([0 1]);     colormap(jet); title('<tcc> over 20 years');; colormap(blues)

junk = [nanmean(trend.wspeed(ocean)) nanstd(trend.wspeed(ocean)) nanmean(abs(trend.wspeed(ocean))) max(abs(trend.wspeed(ocean)))];
fprintf(1,'ocean wspeed trends : %8.4f +/- %8.4f m/s/yr; mean(abs(trend.wspeed)) = %8.4f max(abs(trend.wspeed)) = %8.4f \n',junk); 
figure(4); clf; scatter_coast(p0.rlon,p0.rlat,50,trend.wspeed); caxis([-1 +1]/10); colormap(usa2); title('<speed> trend over 20 years')

junk = [nanmean(trend.stemp(ocean)) nanstd(trend.stemp(ocean)) nanmean(abs(trend.stemp(ocean))) max(abs(trend.stemp(ocean)))];
fprintf(1,'ocean stemp trends : %8.4f +/- %8.4f K/yr; mean(abs(trend.stemp)) = %8.4f max(abs(trend.stemp)) = %8.4f \n',junk); 
figure(5); clf; scatter_coast(p0.rlon,p0.rlat,50,trend.stemp); caxis([-1 +1]/10); colormap(usa2); title('<stemp> trend over 20 years')

junk = [nanmean(trend.tcc(ocean)) nanstd(trend.tcc(ocean)) nanmean(abs(trend.tcc(ocean))) max(abs(trend.tcc(ocean)))];
fprintf(1,'ocean tcc trends : %8.4f +/- %8.4f []/yr; mean(abs(trend.tcc)) = %8.4f max(abs(trend.tcc)) = %8.4f \n',junk); 
figure(6); clf; scatter_coast(p0.rlon,p0.rlat,50,trend.tcc); caxis([-1 +1]/100); colormap(usa2); title('<tcc> trend over 20 years')

% save wspeed_2002_09_2022_08.mat yysave mmsave wspeedsave trend trend_err
saver = ['save wspeed_2002_09_' num2str(yyEnd) '_' num2str(mmEnd,'%02d') '.mat  yysave mmsave wspeedsave stempsave tccsave trend trend_err yyEnd mmEnd'];
eval(saver)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pa = {{'profiles','rtime','seconds since 1993'}};
[pX,pa] = rtp_add_emis(p0,pa);

pnew = p0;
%pnew.wspeed = pnew.wspeed + 0.025; 
pnew.wspeed = pnew.wspeed + trend.wspeed; 
pnew.stemp  = pnew.stemp  + trend.stemp; 
[pnew,pa] = rtp_add_emis(pnew,pa);

figure(7); clf; 
plot(pnew.efreq,pnew.emis - pX.emis);

ocean = find(pnew.landfrac == 0);
plot(pnew.efreq(:,ocean),pnew.emis(:,ocean) - pX.emis(:,ocean));
plot(pnew.efreq(:,ocean),nanmean(pnew.emis(:,ocean) - pX.emis(:,ocean),2));
plot(pnew.efreq(:,ocean),nanmean(abs(pnew.emis(:,ocean) - pX.emis(:,ocean)),2));

%%%%%%%%%%%%%%%%%%%%%%%%%

clear emissivitychange

emissivitychange.landfrac = pnew.landfrac;
emissivitychange.ocean    = ocean;

emissivitychange.p0.stemp = pX.stemp;
emissivitychange.p0.wspeed = pX.wspeed;
emissivitychange.p0.efreq = pX.efreq;
emissivitychange.p0.emis  = pX.emis;
emissivitychange.p0.nemis = pX.nemis;

emissivitychange.pnew.stemp = pnew.stemp;
emissivitychange.pnew.wspeed = pnew.wspeed;
emissivitychange.pnew.efreq = pnew.efreq;
emissivitychange.pnew.emis  = pnew.emis;
emissivitychange.pnew.nemis = pnew.nemis;

saver = ['save wspeed_2002_09_' num2str(yyEnd) '_' num2str(mmEnd,'%02d') '.mat  yysave mmsave wspeedsave stempsave tccsave trend trend_err emissivitychange yyEnd mmEnd'];
eval(saver)
