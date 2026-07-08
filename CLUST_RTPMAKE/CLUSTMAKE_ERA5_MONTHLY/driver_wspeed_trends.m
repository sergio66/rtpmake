%% cp -a /home/sergio/MATLABCODE/CAMEL_emissivity/driver_wspeed_trends.m .
addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/TIME
addpath /home/sergio/MATLABCODE/TIME/MOTTELER_CCAST_MOTMSC_TIME
addpath /home/sergio/MATLABCODE/COLORMAP
addpath /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/StrowCodeforTrendsAndAnomalies/
addpath /asl/rtp_prod2/emis/
addpath /home/sergio/MATLABCODE/matlib/rtp_prod2/emis/
addpath /home/sergio/MATLABCODE/matlib/rtp_prod2/util/

use_this_rtp = 'RTP/summary_19years_all_lat_all_lon_2002_2021_monthlyERA5.rp.rtp';
[h0,ha,p00,pa] = rtpread(['/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/' use_this_rtp]);
p0 = p00;

yyEnd = 2022; mmEnd = 08;
yyEnd = 2024; mmEnd = 06;

iCnt = 0;
for yy = 2002 : yyEnd
  mmS = 1;
  mmE = 12;

  if yy == 2002
    mmS = 9;
  %elseif yy == 2022
  %  mmE = 8;
  %elseif yy == 2024
  %  mmE = 6;
  elseif yy == yyEnd
    mmE = mmEnd;
  end

  for mm = mmS : mmE
    iCnt = iCnt + 1;
    fprintf(1,'iCnt = %3i  : %4i/%2i \n',iCnt,yy,mm)
    a = read_netcdf_lls(['/asl/models/era5_monthly/' num2str(yy) '/' num2str(yy) '-' num2str(mm,'%02i') '_sfc.nc']);
    %% see ~/MATLABCODE/RTPMAKE/CLUST_RTPMAKE/GRIB/grib_interpolate_era5.m  
    [X,Y] = ndgrid(a.latitude,a.longitude);
    iX = flipud(X); iY = flipud(Y);

    stemp = a.skt;
    stemp = mean(stemp,3);
    F = griddedInterpolant(iX,iY,flipud(stemp'));
    p0.stemp = F(p0.rlat,wrapTo360(p0.rlon));
    scatter_coast(p0.rlon,p0.rlat,50,p0.stemp); caxis([200 320]); colormap jet; pause(0.1)

    tcc = a.tcc;
    tcc = mean(tcc,3);
    F = griddedInterpolant(iX,iY,flipud(tcc'));
    p0.tcc = F(p0.rlat,wrapTo360(p0.rlon));
    scatter_coast(p0.rlon,p0.rlat,50,p0.tcc); caxis([200 320]); colormap jet; pause(0.1)

    wspd = sqrt(a.u10.^2 + a.v10.^2);
    wspd_max  = max(wspd,3);
    wspd_min  = min(wspd,3);
    wspd_mean = mean(wspd,3);
    %pcolor(a.longitude,a.latitude,wspd_mean'); shading flat; colorbar; colormap jet; caxis([0 15])  
    F = griddedInterpolant(iX,iY,flipud(wspd_mean'));
    p0.wspeed = F(p0.rlat,wrapTo360(p0.rlon));
    scatter_coast(p0.rlon,p0.rlat,50,p0.wspeed); caxis([0 15]); colormap jet; pause(0.1)
  
    yysave(iCnt) = yy;
    mmsave(iCnt) = mm;
    wspeedsave(iCnt,:) = p0.wspeed;
    stempsave(iCnt,:)  = p0.stemp;
    tccsave(iCnt,:)    = p0.tcc;

  end
end

iCntMax = iCnt;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

scatter_coast(p0.rlon,p0.rlat,50,nanmean(wspeedsave,1)); caxis([0 15]); colormap jet; title('<speed> over 20 years')

% save wspeed_2002_09_2022_08.mat yysave mmsave wspeedsave
saver = ['save wspeed_2002_09_' num2str(yyEnd) '_' num2str(mmEnd,'%02d') '.mat  yysave mmsave wspeedsave stempsave tccsave iCntMax yyEnd mmEnd'];
eval(saver)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

compute_wspeed_trends

disp('now look at /home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/wind_speed_changes.m')
disp('now look at /home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/wind_speed_changes.m')
disp('now look at /home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/wind_speed_changes.m')
