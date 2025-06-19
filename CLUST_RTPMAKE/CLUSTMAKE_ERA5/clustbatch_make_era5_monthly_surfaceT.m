addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/PLOTTER
addpath /home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/

disp('also see /home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/clustbatch_make_L3_airsV7_climcaps_surfaceT.m')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

do_XX_YY_from_X_Y

JOB = str2num(getenv('SLURM_ARRAY_TASK_ID'));
if length(JOB) == 0
  JOB = 120;
end

yy = 2002;
mm = 08;
for ii = 1 : JOB
  mm = mm + 1;
  if mm > 12
    yy = yy + 1;
    mm = 1;
  end
  ysave(ii) = yy;
  msave(ii) = mm;
  fprintf(1,'%3i %4i/%2i \n',ii,ysave(ii),msave(ii))
end

fname = ['/home/sergio//asl/asl/ERA5_monthly/' num2str(yy) '/' num2str(yy) '-' num2str(mm,'%02d') '_sfc.nc'];
a = read_netcdf_lls(fname);
[Y,I] = sort(a.latitude);
a.latitude = a.latitude(I);
a.skt = a.skt(:,I,:);

a.longitude = wrapTo180(a.longitude);
[Y,I] = sort(a.longitude);
a.longitude = a.longitude(I);
a.skt = a.skt(I,:,:);

for jj = 1 : 64
  fprintf(1,'latbin %2i of 64 \n',jj);
  for ii = 1 : 72
    tile = (jj-1)*72 + ii;

    rlonIJ(ii,jj)  = rlon(ii);    
    rlatIJ(ii,jj)  = rlat(jj);    

    xlim1 = rlon73(ii);
    xlim2 = rlon73(ii+1);
    ylim1 = rlat65(jj);
    ylim2 = rlat65(jj+1);
    booX = find(a.longitude >= xlim1 & a.longitude < xlim2);
    booY = find(a.latitude >= ylim1 & a.latitude < ylim2);
    mean_rlonIJ(ii,jj)  = nanmean(a.longitude(booX));
    mean_rlatIJ(ii,jj)  = nanmean(a.latitude(booY));    
    wah = a.skt(booX,booY,:);
    for tt = 1 : 8
      wahx = squeeze(wah(:,:,tt));
      wahx = wahx(:);
      mean_stempIJ(ii,jj,tt) = nanmean(wahx);
    end
    
    xctr1 = rlon(ii);
    yctr1 = rlat(jj);
    booX = find(a.longitude >= xctr1,1);
    booY = find(a.latitude >= yctr1,1);
    cntr_rlonIJ(ii,jj)  = a.longitude(booX);
    cntr_rlatIJ(ii,jj)  = a.latitude(booY);    
    wah = a.skt(booX,booY,:);
    for tt = 1 : 8
      wahx = squeeze(wah(:,:,tt));
      wahx = wahx(:);    
      cntr_stempIJ(ii,jj,tt) = wahx;
    end
  end
end

addpath /home/sergio/MATLABCODE/COLORMAP
figure(1); scatter_coast(rlonIJ,rlatIJ,50,squeeze(nanmean(mean_stempIJ,3))); colormap jet
figure(2); scatter_coast(rlonIJ,rlatIJ,50,squeeze(nanmean(cntr_stempIJ,3))); colormap jet
figure(3); scatter_coast(rlonIJ,rlatIJ,50,squeeze(nanmean(mean_stempIJ - cntr_stempIJ,3))); colormap(usa2); caxis([-1 +1])
figure(4); scatter_coast(rlonIJ,rlatIJ,50,rlatIJ - mean_rlatIJ); colormap(usa2); caxis([-1 +1])
figure(4); scatter_coast(rlonIJ,rlatIJ,50,rlatIJ - cntr_rlatIJ); colormap(usa2); caxis([-1 +1])
figure(4); scatter_coast(rlonIJ,rlatIJ,50,rlonIJ - mean_rlonIJ); colormap(usa2); caxis([-1 +1])
figure(4); scatter_coast(rlonIJ,rlatIJ,50,rlonIJ - cntr_rlonIJ); colormap(usa2); caxis([-1 +1])

comment = 'see clustbatch_make_merra2_monthly_surfaceT.m';
saver = ['save SKT_TIMESERIES_2002_09_to_2022_08/surfaceT_' num2str(yy) '_' num2str(mm,'%02d') '.mat cntr_stempIJ cntr_rlonIJ cntr_rlatIJ mean_stempIJ mean_rlonIJ mean_rlatIJ rlonIJ rlatIJ comment'];
eval(saver)
