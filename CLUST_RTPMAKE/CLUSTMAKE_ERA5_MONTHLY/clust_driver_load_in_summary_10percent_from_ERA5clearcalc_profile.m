addpath /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/StrowCodeforTrendsAndAnomalies
addpath /asl/matlib/maps
addpath /asl/matlib/h4tools
addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/PLOTTER
addpath /home/sergio/MATLABCODE/COLORMAP
addpath /home/sergio/MATLABCODE/COLORMAP/LLS

%% same as driver_load_in_summary_10percent_from_ERA5clearcal.m except we can handle T or WV profiles

%% clust_loop_make_monthly_tile_273points.m made these files, then clust_find_hottest_10percent_from_ERA5clearcalcs.m finds the hottest 10 percent and avg
%% 20 years = 240 months : divide into 24 files, each 10 steps long

JOB = str2num(getenv('SLURM_ARRAY_TASK_ID'));  %% about 5 regions, two profiles (T,WV) so JOB = 1 : 10;
if length(JOB) == 0
  JOB = 01;
end

if JOB <= 5
  iTorWV = +1;  %% T profile trends
else
  iTorWV = -1;  %% WV profile trends
end

%%%%%%%%%%%%%%%%%%%%%%%%%

%% this cannot do ALL profiles so just do a set at a time ef 0001-1536, 1537-3072, 3073-4608
load latB64.mat
rlat65 = latB2; rlon73 = -180 : 5 : +180;
rlon = -180 : 5 : +180;  rlat = latB2;
rlon = 0.5*(rlon(1:end-1)+rlon(2:end));
rlat = 0.5*(rlat(1:end-1)+rlat(2:end));

[Y,X] = meshgrid(rlat,rlon);
Y0 = Y;
X0 = X;
Y = Y(:); X = X(:);

%{
  ii = 1;
  fin = ['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Monthly/summary_with_Tprofile_' num2str(ii,'%02d') '.mat'];
  loader = ['a = load(''' fin ''');'];
  fprintf(1,'%2i %s \n',ii,fin)
%% 1 /asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Monthly/summary_with_Tprofile_01.mat
  eval(loader)

moo = squeeze(nanmean(a.stemp,1));
moo = squeeze(nanmean(moo,1));
moo = squeeze(nanmean(moo,1));
scatter_coast(X,Y,moo)

colormap jet
simplemap_thickcoast(Y,X,moo',5)

%}
%%%%%%%%%%%%%%%%%%%%%%%%%
Y000 = Y0(1,:);

iSub = -2; thepointsM2 = find(Y <= -60);
iSub = -1; thepointsM1 = find(Y > -60 & Y <= -30);
iSub = +0; thepoints00 = find(Y > -30 & Y <= +30);
iSub = +1; thepointsP1 = find(Y > +30 & Y <= +60);
iSub = +2; thepointsP2 = find(Y > +60);
whos thepoints*; 
fprintf(1,'Sum of SP + SML + T + NML + NP points = %4i \n',length(thepointsM2) + length(thepointsM1) + length(thepoints00) + length(thepointsP1) + length(thepointsP2))

iSub = 0; %% tropics

iSub = mod(JOB,5);
iSub(iSub == 0) = 5;
iSub = iSub - 3;

if iSub == -2
  thepoints = thepointsM2;
  disp('processing Southern Polar')
elseif iSub == -1
  thepoints = thepointsM1;
  disp('processing Southern MidLat')
elseif iSub == 00
  thepoints = thepoints00;
  disp('processing Tropics')
elseif iSub == +1
  thepoints = thepointsP1;
  disp('processing Northern MidLat')
elseif iSub == +2
  thepoints = thepointsP2;
  disp('processing Northern Polar')
end

mooY = unique(Y(thepoints));
[Ylat,I1,yindex] = intersect(mooY,Y000);

fprintf(1,'JOB iSub iTorWV = %2i %2i %2i \n',JOB,iSub,iTorWV)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bt1231cld = nan(240,21,12,length(thepoints));
bt1231clr = nan(240,21,12,length(thepoints));
landfrac  = nan(240,21,12,length(thepoints));
stemp     = nan(240,21,12,length(thepoints));
if iTorWV == 1
  ptemp = nan(240,21,12,101,length(thepoints));
elseif iTorWV == -1
  gas_1 = nan(240,21,12,101,length(thepoints));
end

for ii = 1 : 24
  if iTorWV == 1
    fin = ['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Monthly/summary_with_Tprofile_' num2str(ii,'%02d') '.mat'];
  elseif iTorWV == -1
    fin = ['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Monthly/summary_with_WVprofile_' num2str(ii,'%02d') '.mat'];
  end
  loader = ['a = load(''' fin ''');'];
  fprintf(1,'%2i %s \n',ii,fin)
  eval(loader)
  jobs = a.JOBB;

  bt1231cld(jobs,:,:,:) = a.bt1231cld(:,:,:,thepoints);
  bt1231clr(jobs,:,:,:) = a.bt1231clr(:,:,:,thepoints);;
  landfrac(jobs,:,:,:)  = a.landfrac(:,:,:,thepoints);;
  stemp(jobs,:,:,:)     = a.stemp(:,:,:,thepoints);;
  if iTorWV == 1
    ptemp(jobs,:,:,:,:)     = a.ptemp(:,:,:,:,thepoints);;
  elseif iTorWV == -1
    gas_1(jobs,:,:,:,:)     = a.gas_1(:,:,:,:,thepoints);;
  end
end

midpt = floor(length(thepoints)/2);
timeseries_btcld = squeeze(bt1231cld(:,10,6,midpt));
timeseries_btclr = squeeze(bt1231clr(:,10,6,midpt));
timeseries_stemp = squeeze(stemp(:,10,6,midpt));
plot(1:240,timeseries_stemp,'b',1:240,timeseries_btclr,'g',1:240,timeseries_btcld,'r')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[h,ha,p,pa] = rtpread('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/summary_17years_all_lat_all_lon_2002_2019_palts_startSept2002_CLEAR.rtp');

%% figure 1 in trends paper shows BT1231 histogram for 2012/08 to 2012/09 = 10 years on = 120 month

xbt1231cld = squeeze(bt1231cld(120,:,:,:));                   % whos xbt1231cld
xbt1231cld = cat(1,xbt1231cld,squeeze(bt1231cld(121,:,:,:))); % whos xbt1231cld
xbt1231cld = cat(1,xbt1231cld,squeeze(bt1231cld(119,:,:,:))); % whos xbt1231cld
xbt1231cld = cat(1,xbt1231cld,squeeze(bt1231cld(122,:,:,:))); % whos xbt1231cld
xbt1231cld = cat(1,xbt1231cld,squeeze(bt1231cld(118,:,:,:))); % whos xbt1231cld
xbt1231cld = cat(1,xbt1231cld,squeeze(bt1231cld(123,:,:,:))); % whos xbt1231cld
xbt1231cld = cat(1,xbt1231cld,squeeze(bt1231cld(117,:,:,:))); % whos xbt1231cld
[mm,nn,oo] = size(xbt1231cld);
dbt = 200 : 1 : 320; clear histbt1231cld;
for ii = 1 : length(yindex)
  yii = yindex(ii);
  ind = (1:72) + (ii-1)*72;
  junk = xbt1231cld(:,:,ind);
  junk = junk(:);
  histbt1231cld(ii,:) = histc(junk,dbt)/(72*mm*nn);
end 
pcolor(dbt,1:length(yindex),(histbt1231cld)); colormap jet; shading interp; colorbar;
pcolor(dbt,1:length(yindex),log10(histbt1231cld)); colormap jet; shading interp; colorbar;
contourf(dbt,1:length(yindex),log10(histbt1231cld),10); colormap jet; shading interp; colorbar;
contourf(dbt,1:length(yindex),histbt1231cld,10); colormap jet; shading interp; colorbar;
  jett = jet(128); jett(1,:) = 1; colormap(jett)

xbt1231clr = squeeze(bt1231clr(120,:,:,:));
dbt = 200 : 1 : 320; clear histbt1231clr;
for ii = 1 : length(yindex)
  ind = (1:72) + (ii-1)*72;
  junk = xbt1231clr(:,:,ind);
  junk = junk(:);
  histbt1231clr(ii,:) = histc(junk,dbt)/(72*252);
end 
pcolor(dbt,1:length(yindex),histbt1231clr); colormap jet; shading interp; colorbar;
pcolor(dbt,yindex,histbt1231clr); colormap jet; shading interp; colorbar;

xstemp = squeeze(stemp(120,:,:,:));
dbt = 200 : 1 : 320; clear histstemp;
for ii = 1 : length(yindex)
  ind = (1:72) + (ii-1)*72;
  junk = xstemp(:,:,ind);
  junk = junk(:);
  histstemp(ii,:) = histc(junk,dbt)/(72*252);
end 
pcolor(dbt,1:length(yindex),histstemp); colormap jet; shading interp; colorbar;
pcolor(dbt,yindex,histstemp); colormap jet; shading interp; colorbar;

pause

% xbt1231cld = squeeze(bt1231cld(120,:,:,:));
% xbt1231cld = reshape(xbt1231cld,12*21,length(thepoints));
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
  junk = squeeze(bt1231clr(ii,:,:,midpt));
  raaClr(ii,:) = sort(junk(:));

  junk = squeeze(bt1231cld(ii,:,:,midpt));
  raaCld(ii,:) = sort(junk(:));

  junk = squeeze(stemp(ii,:,:,midpt));
  raaSKT(ii,:) = sort(junk(:));
end
figure(1); clf; imagesc(raaClr); colorbar; colormap jet; xlabel('Grid pts'); ylabel('Time'); title('Clr')
figure(2); clf; imagesc(raaCld); colorbar; colormap jet; xlabel('Grid pts'); ylabel('Time'); title('Cld')
figure(3); clf; imagesc(raaSKT); colorbar; colormap jet; xlabel('Grid pts'); ylabel('Time'); title('SKT')

meanstemp = squeeze(nanmean(squeeze(nanmean(squeeze(nanmean(stemp,1)),1)),1));    figure(4); pcolor(reshape(meanstemp,72,length(yindex))'); colorbar; shading flat; colormap jet; title('Stemp mean')
meanCld = squeeze(nanmean(squeeze(nanmean(squeeze(nanmean(bt1231cld,1)),1)),1));  figure(5); pcolor(reshape(meanCld,72,length(yindex))'); colorbar; shading flat; colormap jet; title('BT1231 Cld mean')
meanClr = squeeze(nanmean(squeeze(nanmean(squeeze(nanmean(bt1231clr,1)),1)),1));  figure(6); pcolor(reshape(meanClr,72,length(yindex))'); colorbar; shading flat; colormap jet; title('BT1231 Clr mean')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

indQ90 = 21*12;
indQ90 = ceil(0.9*max(indQ90)) : indQ90;
for ii = 1 : 240
  junk = squeeze(bt1231clr(ii,:,:,midpt));
  junk = sort(junk(:));
  bt1231clr_mean(ii) = mean(junk);
  bt1231clr_Q90(ii)  = mean(junk(indQ90));

  junk = squeeze(bt1231cld(ii,:,:,midpt));
  junk = sort(junk(:));
  bt1231cld_mean(ii) = mean(junk);
  bt1231cld_Q90(ii)  = mean(junk(indQ90));
end
figure(7); plot(1:240,bt1231clr_mean,'c',1:240,bt1231clr_Q90,'b',1:240,bt1231cld_mean,'m',1:240,bt1231cld_Q90,'rx-','linewidth',2); 
  xlabel('Time in months starting Sept 2002')
  legend('Clr mean','Clr Q90','Cld mean','Cld Q90','location','best','fontsize',10);
  title('Need to partition using Q90 cld == rx')

%%%%%%%%%%%%%%%%%%%%%%%%%
find_trends_summary_10percent_from_ERA5clearcalc_T_WV

