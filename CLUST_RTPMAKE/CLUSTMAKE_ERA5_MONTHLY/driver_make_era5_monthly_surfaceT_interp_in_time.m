addpath /home/sergio/MATLABCODE/TIME
addpath /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/StrowCodeforTrendsAndAnomalies/
addpath /home/sergio/MATLABCODE/COLORMAP/LLS
load llsmap5

disp('this is after     clustbatch_make_era5_monthly_surfaceT      has run')

iDorA = +1;   %% desc
iDorA = -1;   %% asc

disp('the overpass time code is from clust_loop_make_16day_tile_center.m')
%% from comment, see /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_For_HowardObs_TimeSeries/driver_loop_get_asc_desc_solzen_time.m
%% now we need to get overpass times and solzen angles, just set scanang to 22 deg
%% VERY VERY IMPORTANT : see driver_fix_thedata_asc_desc_solzen_time_412_64x72.m
oldwrong = load('/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_For_HowardObs_TimeSeries/asc_desc_solzen_time_412_64x72.mat');  %% this had THREE skips of data so bad offsets
if iDorA == 1
  load('/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_For_HowardObs_TimeSeries/asc_desc_solzen_time_412_64x72_fix_desc.mat');
elseif iDorA == -1
  load('/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_For_HowardObs_TimeSeries/asc_desc_solzen_time_412_64x72_fix_asc.mat');
end

if iDorA > 0
  rtime = squeeze(nanmean(thedata.rtime_desc,3));
  hour  = squeeze(nanmean(thedata.hour_desc,3));  
elseif iDorA < 0
  rtime = squeeze(nanmean(thedata.rtime_asc,3));
  hour  = squeeze(nanmean(thedata.hour_asc,3));  
end
figure(4); clf; pcolor(hour'); colorbar; title('overpass hh'); colormap jet
disp('ret to ontinue'); pause

times = 0 : 3 : 21;
yy = 2002;
mm = 08;
for iStep = 1 : 12*20
  mm = mm + 1;
  if mm > 12
    yy = yy + 1;
    mm = 1;
  end
  ysave(iStep) = yy;
  msave(iStep) = mm;
  doy(iStep) = change2days(yy,mm,01,18.00);;
  fprintf(1,'%3i %4i/%2i \n',iStep,ysave(iStep),msave(iStep))
  fsavename = ['SKT_TIMESERIES_2002_09_to_2022_08/surfaceT_' num2str(yy) '_' num2str(mm,'%02d') '.mat'];
  a = load(fsavename);

  for jj = 1 : 64
    for ii = 1 : 72
      data = squeeze(a.cntr_stempIJ(ii,jj,:));
      cntr_stempIJ(ii,jj,iStep) = interp1(times,data,hour(ii,jj),[],'extrap');

      data = squeeze(a.mean_stempIJ(ii,jj,:));
      mean_stempIJ(ii,jj,iStep) = interp1(times,data,hour(ii,jj),[],'extrap');
    end
  end
  rlatIJ = a.rlatIJ;
  rlonIJ = a.rlonIJ;  
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for jj = 1 : 64
  for ii = 1 : 72
    boo = squeeze(cntr_stempIJ(ii,jj,:));
    good = 1 : 240;
    [B stats] = Math_tsfit_lin_robust(doy(good),double(boo(good)),4);
    trendERA5.cntr(ii,jj) = B(2);
    trendERA5.cntr_err(ii,jj) = stats.se(2);

    boo = squeeze(mean_stempIJ(ii,jj,:));
    good = 1 : 240;
    [B stats] = Math_tsfit_lin_robust(doy(good),double(boo(good)),4);
    trendERA5.mean(ii,jj) = B(2);
    trendERA5.mean_err(ii,jj) = stats.se(2);

  end
end

figure(1); clf; scatter_coast(rlonIJ,rlatIJ,100,trendERA5.cntr); title('ERA5 CNTR'); caxis([-1 +1]*0.15); colormap(llsmap5)
figure(2); clf; scatter_coast(rlonIJ,rlatIJ,100,trendERA5.mean); title('ERA5 MEAN'); caxis([-1 +1]*0.15); colormap(llsmap5)
figure(3); clf; scatter_coast(rlonIJ,rlatIJ,100,trendERA5.cntr - trendERA5.mean); title('ERA5 MEAN-CNTR'); caxis([-1 +1]*0.15/10); colormap(llsmap5)
figure(4); clf
  delta = trendERA5.cntr - trendERA5.mean;
  dbt = -15:0.25:+15; dbt = dbt/100; 
  plot(dbt,histc(delta(:),dbt)); xlim([-1 +1]*0.02); plotaxis2;
figure(5); clf
  delta = trendERA5.cntr - trendERA5.mean;
  bad = find(abs(trendERA5.mean) <= 0.15/100);
  delta = 100 * delta./trendERA5.mean;
  delta(bad) = NaN;
  dbt = -100:2:+100;
  plot(dbt,histc(delta(:),dbt)); xlim([-1 +1]*50); plotaxis2;
  
if iDorA > 0
  saver = ['save skt_timeseries_20years_desc.mat cntr_stempIJ mean_stempIJ rlatIJ rlonIJ trendERA5'];
else
  saver = ['save skt_timeseries_20years_asc.mat  cntr_stempIJ mean_stempIJ rlatIJ rlonIJ  trendERA5'];  
end
%% eval(saver)
