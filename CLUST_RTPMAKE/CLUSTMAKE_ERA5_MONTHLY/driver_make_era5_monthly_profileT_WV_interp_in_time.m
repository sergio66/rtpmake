addpath /home/sergio/MATLABCODE/TIME
addpath /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/StrowCodeforTrendsAndAnomalies/
addpath /home/sergio/MATLABCODE/COLORMAP/LLS
load llsmap5

disp('this is after     clustbatch_make_era5_monthly_profileT_WV      has run')

iDorA = -1;   %% asc
iDorA = +1;   %% desc

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
  fsavename = ['SKT_TIMESERIES_2002_09_to_2022_08/profileT_WV_' num2str(yy) '_' num2str(mm,'%02d') '.mat'];
  a = load(fsavename);

  for jj = 1 : 64
    for ii = 1 : 72 
      for ll = 1 : 37
        data = squeeze(a.cntr_TIJ(ii,jj,ll,:));
        cntr_TIJ(ii,jj,ll,iStep) = interp1(times,data,hour(ii,jj),[],'extrap');
        data = squeeze(a.cntr_QIJ(ii,jj,ll,:));
        cntr_QIJ(ii,jj,ll,iStep) = interp1(times,data,hour(ii,jj),[],'extrap');

        data = squeeze(a.mean_TIJ(ii,jj,ll,:));
        mean_TIJ(ii,jj,ll,iStep) = interp1(times,data,hour(ii,jj),[],'extrap');
        data = squeeze(a.mean_QIJ(ii,jj,ll,:));
        mean_QIJ(ii,jj,ll,iStep) = interp1(times,data,hour(ii,jj),[],'extrap');
      end
    end
  end
  rlatIJ = a.rlatIJ;
  rlonIJ = a.rlonIJ;  
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('read in data, now fitting ...')

for jj = 1 : 64
  fprintf(1,'lat bin jj = %2i of 64 \n',jj);
  for ii = 1 : 72
    for ll = 1 : 37
      boo = squeeze(cntr_TIJ(ii,jj,ll,:));
      good = 1 : 240;
      [B stats] = Math_tsfit_lin_robust(doy(good),double(boo(good)),4);
      trendT_ERA5.cntr(ii,jj,ll) = B(2);
      trendT_ERA5.cntr_err(ii,jj,ll) = stats.se(2);
  
      boo = squeeze(mean_TIJ(ii,jj,ll,:));
      good = 1 : 240;
      [B stats] = Math_tsfit_lin_robust(doy(good),double(boo(good)),4);
      trendT_ERA5.mean(ii,jj,ll) = B(2);
      trendT_ERA5.mean_err(ii,jj,ll) = stats.se(2);

      boo = squeeze(cntr_QIJ(ii,jj,ll,:));
      boo = boo./mean(boo);
      good = 1 : 240;
      [B stats] = Math_tsfit_lin_robust(doy(good),double(boo(good)),4);
      trendQ_ERA5.cntr(ii,jj,ll) = B(2);
      trendQ_ERA5.cntr_err(ii,jj,ll) = stats.se(2);
  
      boo = squeeze(mean_QIJ(ii,jj,ll,:));
      boo = boo./mean(boo);
      good = 1 : 240;
      [B stats] = Math_tsfit_lin_robust(doy(good),double(boo(good)),4);
      trendQ_ERA5.mean(ii,jj,ll) = B(2);
      trendQ_ERA5.mean_err(ii,jj,ll) = stats.se(2);
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plevs = load('/home/sergio/MATLABCODE/RTPMAKE/CLUST_RTPMAKE/CLUSTMAKE_ERA5/era5plevs.mat');
plevs = plevs.era5plevs;

if iDorA > 0
  saver = ['save profileT_timeseries_20years_desc.mat cntr_TIJ mean_TIJ cntr_QIJ mean_QIJ rlatIJ rlonIJ trend*_ERA5 plevs'];
  saver = ['save profileT_timeseries_20years_desc.mat rlatIJ rlonIJ trend*_ERA5 plevs'];
else
  saver = ['save profileT_timeseries_20years_asc.mat  cntr_TIJ mean_TIJ cntr_QIJ mean_QIJ rlatIJ rlonIJ trend*_ERA5 plevs'];  
  saver = ['save profileT_timeseries_20years_asc.mat  rlatIJ rlonIJ trend*_ERA5 plevs'];  
end
eval(saver)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1); clf; pcolor(rlatIJ(1,:),plevs,squeeze(nanmean(trendT_ERA5.cntr,1))');                    
  shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','log'); title('ERA5 CNTR dT/dt'); caxis([-1 +1]*0.15);  colormap(llsmap5); colorbar; ylim([10 1000])
figure(2); clf; pcolor(rlatIJ(1,:),plevs,squeeze(nanmean(trendT_ERA5.mean,1))');                    
  shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','log'); title('ERA5 MEAN dT/dt'); caxis([-1 +1]*0.15);  colormap(llsmap5); colorbar; ylim([10 1000])
figure(3); clf; pcolor(rlatIJ(1,:),plevs,squeeze(nanmean(trendT_ERA5.mean - trendT_ERA5.cntr,1))'); 
  shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','log'); title('ERA5 MEAN-CNTR dT/dt'); caxis([-1 +1]*0.015); colormap(llsmap5); colorbar; ylim([10 1000])
figure(4); clf; pcolor(rlatIJ(1,:),plevs,100*squeeze(nanmean(trendT_ERA5.mean - trendT_ERA5.cntr,1))'./squeeze(nanmean(trendT_ERA5.cntr,1))'); 
  shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','linear'); title('ERA5 percent MEAN-CNTR dT/dt'); caxis([-1 +1]*10); colormap(llsmap5); colorbar; ylim([100 1000])
  
figure(1); clf; pcolor(rlatIJ(1,:),plevs,squeeze(nanmean(trendQ_ERA5.cntr,1))');                    
  shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','linear'); title('ERA5 CNTR dfracWV/dt'); caxis([-1 +1]*0.01);  colormap(llsmap5); colorbar; ylim([100 1000])
figure(2); clf; pcolor(rlatIJ(1,:),plevs,squeeze(nanmean(trendQ_ERA5.mean,1))');                    
  shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','linear'); title('ERA5 MEAN dfracWV/dt'); caxis([-1 +1]*0.01);  colormap(llsmap5); colorbar; ylim([100 1000])
figure(3); clf; pcolor(rlatIJ(1,:),plevs,squeeze(nanmean(trendQ_ERA5.mean - trendQ_ERA5.cntr,1))'); 
  shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','linear'); title('ERA5 MEAN-CNTR dfracWV/dt'); caxis([-1 +1]*0.001); colormap(llsmap5); colorbar; ylim([100 1000])
figure(4); clf; pcolor(rlatIJ(1,:),plevs,100*squeeze(nanmean(trendQ_ERA5.mean - trendQ_ERA5.cntr,1))'./squeeze(nanmean(trendQ_ERA5.cntr,1))'); 
  shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','linear'); title('ERA5 percent MEAN-CNTR dfracWV/dt'); caxis([-1 +1]*10); colormap(llsmap5); colorbar; ylim([100 1000])
  
figure(5); 
tpercent = 100*squeeze(nanmean(trendT_ERA5.mean - trendT_ERA5.cntr,1))'./squeeze(nanmean(trendT_ERA5.cntr,1))';
  tbad = find(abs(squeeze(nanmean(trendT_ERA5.cntr,1))) < 0.15/100);
  tpercent(tbad) = NaN;
qpercent = 100*squeeze(nanmean(trendQ_ERA5.mean - trendQ_ERA5.cntr,1))'./squeeze(nanmean(trendQ_ERA5.cntr,1))';
  qbad = find(abs(squeeze(nanmean(trendQ_ERA5.cntr,1))) < 0.01/100);
  qpercent(qbad) = NaN;
delta = -100:0.1:100;
plot(delta,histc(qpercent(:),delta),'b',delta,histc(tpercent(:),delta),'r'); 
plotaxis2; xlim([-1 +1]*10); legend('WV','T','location','best'); title('Zonal Avg Percent diff MEAN-CNTR'); 

figure(6); 
tpercent = 100*(trendT_ERA5.mean - trendT_ERA5.cntr)./trendT_ERA5.cntr;
  tbad = find(abs(trendT_ERA5.cntr) < 0.15/100);
  tpercent(tbad) = NaN;
qpercent = 100*(trendQ_ERA5.mean - trendQ_ERA5.cntr)./trendQ_ERA5.cntr;
  qbad = find(abs(trendQ_ERA5.cntr) < 0.01/100);
  qpercent(tbad) = NaN;
delta = -100:0.1:100;
plot(delta,histc(qpercent(:),delta),'b',delta,histc(tpercent(:),delta),'r'); 
plotaxis2; xlim([-1 +1]*20); legend('WV','T','location','best'); title('Percent diff MEAN-CNTR'); 

figure(7);
tt = trendT_ERA5.mean - trendT_ERA5.cntr; tt = squeeze(nanmean(squeeze(nanmean(tt,1))));
qq = trendQ_ERA5.mean - trendQ_ERA5.cntr; qq = squeeze(nanmean(squeeze(nanmean(qq,1))));
  plot(qq,plevs,'b',tt,plevs,'r'); set(gca,'ydir','reverse')
tt0 = trendT_ERA5.cntr; tt0 = squeeze(nanmean(squeeze(nanmean(tt0,1))));
qq0 = trendQ_ERA5.cntr; qq0 = squeeze(nanmean(squeeze(nanmean(qq0,1))));
  plot(100*qq./qq0,plevs,'b',100*tt./tt0,plevs,'r'); set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([10 1000])
plotaxis2; xlim([-1 +1]*2); legend('WV','T','location','best'); title('Percent diff MEAN-CNTR'); 
