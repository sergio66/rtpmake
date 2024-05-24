%% cp -a /home/sergio/MATLABCODE/CAMEL_emissivity/driver_wspeed_trends.m .
addpath /home/sergio/MATLABCODE

use_this_rtp = 'RTP/summary_19years_all_lat_all_lon_2002_2021_monthlyERA5.rp.rtp';
[h0,ha,p0,pa] = rtpread(['/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/' use_this_rtp]);

iCnt = 0;
for yy = 2002 : 2022
  mmS = 1;
  mmE = 12;
  if yy == 2002
    mmS = 9;
  elseif yy == 2022
    mmE = 8;
  end
  for mm = mmS : mmE
    iCnt = iCnt + 1;
    fprintf(1,'iCnt = %3i  : %4i/%2i \n',iCnt,yy,mm)
    a = read_netcdf_lls(['/asl/models/era5_monthly/' num2str(yy) '/' num2str(yy) '-' num2str(mm,'%02i') '_sfc.nc']);
    wspd = sqrt(a.u10.^2 + a.v10.^2);
    wspd_max  = max(wspd,3);
    wspd_min  = min(wspd,3);
    wspd_mean = mean(wspd,3);
    %pcolor(a.longitude,a.latitude,wspd_mean'); shading flat; colorbar; colormap jet; caxis([0 15])
  
    %% see ~/MATLABCODE/RTPMAKE/CLUST_RTPMAKE/GRIB/grib_interpolate_era5.m  
    [X,Y] = ndgrid(a.latitude,a.longitude);
    iX = flipud(X); iY = flipud(Y);
  
    F = griddedInterpolant(iX,iY,flipud(wspd_mean'));
  
    p0.wspeed = F(p0.rlat,wrapTo360(p0.rlon));
    scatter_coast(p0.rlon,p0.rlat,50,p0.wspeed); caxis([0 15]); colormap jet; pause(0.1)
  
    yysave(iCnt) = yy;
    mmsave(iCnt) = mm;
    wspeedsave(iCnt,:) = p0.wspeed;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iCnt = 1 : 240
  daysSince2002(iCnt) = change2days(yysave(iCnt),mmsave(iCnt),15,2002);
end

scatter_coast(p0.rlon,p0.rlat,50,nanmean(wspeedsave,1)); caxis([0 15]); colormap jet; title('<speed> over 20 years')

% save wspeed_2002_09_2022_08.mat yysave mmsave wspeedsave

for ii = 1 : 4608
  [B, err, stats] = Math_tsfit_lin_robust(daysSince2002,wspeedsave(:,ii),4);
  trend(ii) = B(2);
  trend_err(ii) = err.se(2);
end

addpath /asl/rtp_prod2/emis/
[p0.salti,p0.landfrac] = usgs_deg10_dem(p0.rlat, p0.rlon);
ocean = find(p0.landfrac == 0); whos ocean
plot(-0.1:0.01:+0.1,histc(trend(ocean),-0.1:0.01:+0.1)); grid

junk = [nanmean(trend(ocean)) nanstd(trend(ocean)) nanmean(abs(trend(ocean))) max(abs(trend(ocean)))];
fprintf(1,'ocean wspeed trends : %8.4f +/- %8.4f m/s; mean(abs(trend)) = %8.4f max(abs(trend)) = %8.4f \n',junk); 

scatter_coast(p0.rlon,p0.rlat,50,trend); caxis([-1 +1]/10); colormap(usa2); title('<speed> trend over 20 years')

% save wspeed_2002_09_2022_08.mat yysave mmsave wspeedsave trend trend_err

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


pa = {{'profiles','rtime','seconds since 1993'}};
[pX,pa] = rtp_add_emis(p0,pa);

pnew = p0;
%pnew.wspeed = pnew.wspeed + 0.025; 
pnew.wspeed = pnew.wspeed + trend; 
[pnew,pa] = rtp_add_emis(pnew,pa);

plot(pnew.efreq,pnew.emis - pX.emis);

ocean = find(pnew.landfrac == 0);
plot(pnew.efreq(:,ocean),pnew.emis(:,ocean) - pX.emis(:,ocean));
plot(pnew.efreq(:,ocean),nanmean(pnew.emis(:,ocean) - pX.emis(:,ocean),2));
plot(pnew.efreq(:,ocean),nanmean(abs(pnew.emis(:,ocean) - pX.emis(:,ocean)),2));
