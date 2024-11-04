addpath /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/StrowCodeforTrendsAndAnomalies
addpath /asl/matlib/maps
addpath /asl/matlib/h4tools
addpath /home/sergio/MATLABCODE/COLORMAP
addpath /home/sergio/MATLABCODE/COLORMAP/LLS

bt1231cld = nan(240,21,12,4608);
bt1231clr = nan(240,21,12,4608);
landfrac  = nan(240,21,12,4608);
stemp     = nan(240,21,12,4608);

%% clust_loop_make_monthly_tile_273points.m made these files, then clust_find_hottest_10percent_from_ERA5clearcalcs.m finds the hottest 10 percent and avg
%% 20 years = 240 months : divide into 24 files, each 10 steps long
for ii = 1 : 24
  fin = ['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Monthly/summary_' num2str(ii,'%02d') '.mat'];
  loader = ['a = load(''' fin ''');'];
  fprintf(1,'%2i %s \n',ii,fin)
  eval(loader)
  jobs = a.JOBB;
  bt1231cld(jobs,:,:,:) = a.bt1231cld;
  bt1231clr(jobs,:,:,:) = a.bt1231clr;
  landfrac(jobs,:,:,:)  = a.landfrac;
  stemp(jobs,:,:,:)     = a.stemp;
end

timeseries_btcld = squeeze(bt1231cld(:,10,6,2000));
timeseries_btclr = squeeze(bt1231clr(:,10,6,2000));
timeseries_stemp = squeeze(stemp(:,10,6,2000));
plot(1:240,timeseries_stemp,'b',1:240,timeseries_btclr,'g',1:240,timeseries_btcld,'r')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[h,ha,p,pa] = rtpread('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/summary_17years_all_lat_all_lon_2002_2019_palts_startSept2002_CLEAR.rtp');

%% figure 1 in trends paper shows BT1231 histogram for 2012/08 to 2012/09 = 10 years on = 120 month

xbt1231cld = squeeze(bt1231cld(120,:,:,:)); whos xbt1231cld
xbt1231cld = cat(1,xbt1231cld,squeeze(bt1231cld(121,:,:,:))); whos xbt1231cld
xbt1231cld = cat(1,xbt1231cld,squeeze(bt1231cld(119,:,:,:))); whos xbt1231cld
xbt1231cld = cat(1,xbt1231cld,squeeze(bt1231cld(122,:,:,:))); whos xbt1231cld
xbt1231cld = cat(1,xbt1231cld,squeeze(bt1231cld(118,:,:,:))); whos xbt1231cld
xbt1231cld = cat(1,xbt1231cld,squeeze(bt1231cld(123,:,:,:))); whos xbt1231cld
xbt1231cld = cat(1,xbt1231cld,squeeze(bt1231cld(117,:,:,:))); whos xbt1231cld
[mm,nn,oo] = size(xbt1231cld);
dbt = 200 : 1 : 320; clear histbt1231cld;
for ii = 1 : 64
  ind = (1:72) + (ii-1)*72;
  junk = xbt1231cld(:,:,ind);
  junk = junk(:);
  histbt1231cld(ii,:) = histc(junk,dbt)/(72*mm*nn);
end 
pcolor(dbt,1:64,(histbt1231cld)); colormap jet; shading interp; colorbar;
pcolor(dbt,1:64,log10(histbt1231cld)); colormap jet; shading interp; colorbar;
contourf(dbt,1:64,log10(histbt1231cld),10); colormap jet; shading interp; colorbar;
contourf(dbt,1:64,histbt1231cld,10); colormap jet; shading interp; colorbar;
  jett = jet(128); jett(1,:) = 1; colormap(jett)

xbt1231clr = squeeze(bt1231clr(120,:,:,:));
dbt = 200 : 1 : 320; clear histbt1231clr;
for ii = 1 : 64
  ind = (1:72) + (ii-1)*72;
  junk = xbt1231clr(:,:,ind);
  junk = junk(:);
  histbt1231clr(ii,:) = histc(junk,dbt)/(72*252);
end 
pcolor(dbt,1:64,(histbt1231clr)); colormap jet; shading interp; colorbar;

xstemp = squeeze(stemp(120,:,:,:));
dbt = 200 : 1 : 320; clear histstemp;
for ii = 1 : 64
  ind = (1:72) + (ii-1)*72;
  junk = xstemp(:,:,ind);
  junk = junk(:);
  histstemp(ii,:) = histc(junk,dbt)/(72*252);
end 
pcolor(dbt,1:64,(histstemp)); colormap jet; shading interp; colorbar;

pause

% xbt1231cld = squeeze(bt1231cld(120,:,:,:));
% xbt1231cld = reshape(xbt1231cld,12*21,4608);
% xbt1231cld = xbt1231cld';
% xbt1231cld = reshape(xbt1231cld',72,64,252);;
% xbt1231cld = permute(xbt1231cld,[1 3 2]);
% xbt1231cld = reshape(xbt1231cld,72*252,64);
% for ii = 1 : 64
%   ind = (1:72) + (ii-1)*72;
%   junk = xbt1231cld(:,:,ind);
%   histbt1231cld(ii,:) = histc(xbt1231cld(:,ii),dbt)/(72*252);
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ii = 1 : 240
  junk = squeeze(bt1231clr(ii,:,:,1000));
  raaClr(ii,:) = sort(junk(:));

  junk = squeeze(bt1231cld(ii,:,:,1000));
  raaCld(ii,:) = sort(junk(:));

  junk = squeeze(stemp(ii,:,:,1000));
  raaSKT(ii,:) = sort(junk(:));
end
figure(1); clf; imagesc(raaClr); colorbar; colormap jet; xlabel('Grid pts'); ylabel('Time'); title('Clr')
figure(2); clf; imagesc(raaCld); colorbar; colormap jet; xlabel('Grid pts'); ylabel('Time'); title('Cld')
figure(3); clf; imagesc(raaSKT); colorbar; colormap jet; xlabel('Grid pts'); ylabel('Time'); title('SKT')

meanstemp = squeeze(nanmean(squeeze(nanmean(squeeze(nanmean(stemp,1)),1)),1));    figure(4); pcolor(reshape(meanstemp,72,64)'); colorbar; shading flat; colormap jet; title('Stemp mean')
meanCld = squeeze(nanmean(squeeze(nanmean(squeeze(nanmean(bt1231cld,1)),1)),1));  figure(5); pcolor(reshape(meanCld,72,64)'); colorbar; shading flat; colormap jet; title('BT1231 Cld mean')
meanClr = squeeze(nanmean(squeeze(nanmean(squeeze(nanmean(bt1231clr,1)),1)),1));  figure(6); pcolor(reshape(meanClr,72,64)'); colorbar; shading flat; colormap jet; title('BT1231 Clr mean')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

indQ90 = 21*12;
indQ90 = ceil(0.9*max(indQ90)) : indQ90;
for ii = 1 : 240
  junk = squeeze(bt1231clr(ii,:,:,1000));
  junk = sort(junk(:));
  bt1231clr_mean(ii) = mean(junk);
  bt1231clr_Q90(ii)  = mean(junk(indQ90));
end
figure(7); plot(1:240,bt1231clr_mean,1:240,bt1231clr_Q90)

%%%%%%%%%%%%%%%%%%%%%%%%%
find_trends_summary_10percent_from_ERA5clearcalc

