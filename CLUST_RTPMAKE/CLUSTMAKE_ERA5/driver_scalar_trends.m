%% cp -a /home/sergio/MATLABCODE/CAMEL_emissivity/driver_wspeed_trends.m .
addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/COLORMAP
addpath /home/sergio/MATLABCODE/COLORMAP/LLS
addpath /asl/matlib/maps
addpath /asl/matlib/h4tools

load latB64.mat
  rlat65 = latB2; rlon73 = -180 : 5 : +180;
  rlon = -180 : 5 : +180;  rlat = latB2;
  rlon = 0.5*(rlon(1:end-1)+rlon(2:end));
  rlat = 0.5*(rlat(1:end-1)+rlat(2:end));

load llsmap5

iWgt = +1;  %% use the [ii ii+1] and [ii+1 ii+2] bins 
iWgt = -1;  %% use only the [ii ii+1] bin, wieght 1, which believe it or not is apparently what fill_era5_monthly does

iWgt

use_this_rtp = 'RTP/summary_19years_all_lat_all_lon_2002_2021_monthlyERA5.rp.rtp';
[h0,ha,p0,pa] = rtpread(['/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/' use_this_rtp]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% now we need to get overpass times and solzen angles, just set scanang to 22 deg, 
%% but remeber this code is for monthly (so 12 sets/year while this file is 16 day so 23 sets/year ... and I only got stuff till 2019 ... so do some averaging!!!!
load('/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_For_HowardObs_TimeSeries/asc_desc_solzen_time_412_64x72.mat');
  utchhx = squeeze(nanmean(thedata.hour_desc,3))'; rlon = squeeze(nanmean(thedata.rlon_desc,3))'; lshDx = utchour2localtime(utchhx,rlon); lshDx = mod(lshDx,24); 
  utchhx = squeeze(nanmean(thedata.hour_asc,3))';  rlon = squeeze(nanmean(thedata.rlon_asc,3))';  lshAx = utchour2localtime(utchhx,rlon); lshAx = mod(lshAx,24); 
lshDx = lshDx';
lshAx = lshAx';
lshDx = lshDx(:);
lshAx = lshAx(:);

load('/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_For_HowardObs_TimeSeries/asc_desc_solzen_time_457_64x72.mat');
  utchh = squeeze(nanmean(thedata.hour_desc,3))'; rlon = squeeze(nanmean(thedata.rlon_desc,3))'; lshD = utchour2localtime(utchh,rlon); lshD = mod(lshD,24); 
  utchh = squeeze(nanmean(thedata.hour_asc,3))';  rlon = squeeze(nanmean(thedata.rlon_asc,3))';  lshA = utchour2localtime(utchh,rlon); lshA = mod(lshA,24); 
lshD = lshD';
lshA = lshA';
lshD = lshD(:);
lshA = lshA(:);
figure(1); clf; scatter_coast(p0.rlon,p0.rlat,50,lshD-lshDx); colorbar; title('Desc local hour New-Old'); colormap(usa2); caxis([-1 +1]/100); shading interp
figure(2); clf; scatter_coast(p0.rlon,p0.rlat,50,lshA-lshAx); colorbar; title('Asc local hour New-Old'); colormap(usa2); caxis([-1 +1]/100); shading interp

figure(1); clf; scatter_coast(p0.rlon,p0.rlat,50,lshD); colorbar; title('Desc local hour'); colormap jet; shading interp
figure(2); clf; scatter_coast(p0.rlon,p0.rlat,50,lshA); colorbar; title('Asc local hour');  colormap jet; shading interp
figure(3); clf; hhx = -1:0.25:+25; plot(hhx,histc(lshD(:),hhx),hhx,histc(lshA(:),hhx),'linewidth',2); hl = legend('Desc','Asc','location','best'); grid; xlim([-1 25])

eraHrs = 0 : 3: 24;
for ii = 1 : length(eraHrs) - 1
  boo = find(lshD >= eraHrs(ii) & lshD < eraHrs(ii+1));
  indD(boo) = 1;
  indD(boo) = ii;
  if iWgt == -1
    wgtD(boo) = 1;
  else
    wgtD(boo) = 1 - (lshD(boo)-eraHrs(ii))/3; %% counterintuitive but suppose lshD(ii) = 1 ... so we are looking at [0 3] bin and weight should be 2/3 for utc=0 and 1/3 for utc=3
  end

  boo = find(lshA >= eraHrs(ii) & lshA < eraHrs(ii+1));
  indA(boo) = 1;
  indA(boo) = ii;
  wgtA(boo) = ii;
  if iWgt == -1
    wgtA(boo) = +1;
  else
    wgtA(boo) = 1 - (lshA(boo)-eraHrs(ii))/3; %% counterintuitive but suppose lshA(ii) = 8 ... so we are looking at [6 9] bin and weight should be 1/3 for utc=6 and 2/3 for utc=9
  end
end

for ii = 1 : 8
  boo = find(indA == ii);
  cntA(ii) = length(boo);
  boo = find(indD == ii);
  cntD(ii) = length(boo);
end

% a.siconc = zeros(1440,721);
% a.sp     = zeros(1440,721);
% a.tcc    = zeros(1440,721);
% a.u10    = zeros(1440,721);
% a.v10    = zeros(1440,721);
% a.skt    = zeros(1440,721);
% a.t2m    = zeros(1440,721);
% a.d2m    = zeros(1440,721); 

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
    xa = read_netcdf_lls(['/asl/models/era5_monthly/' num2str(yy) '/' num2str(yy) '-' num2str(mm,'%02i') '_sfc.nc']);
    %% see ~/MATLABCODE/RTPMAKE/CLUST_RTPMAKE/GRIB/grib_interpolate_era5.m  
    [X,Y] = ndgrid(xa.latitude,xa.longitude);
    iX = flipud(X); iY = flipud(Y);
    yysave(iCnt) = yy;
    mmsave(iCnt) = mm;

    for ii = 1 : 8
      a1.siconc = squeeze(xa.siconc(:,:,ii));
      a1.sp     = squeeze(xa.sp(:,:,ii));
      a1.tcc    = squeeze(xa.tcc(:,:,ii));
      a1.u10    = squeeze(xa.u10(:,:,ii));
      a1.v10    = squeeze(xa.v10(:,:,ii));
      a1.skt    = squeeze(xa.skt(:,:,ii));
      a1.d2m    = squeeze(xa.d2m(:,:,ii));
      a1.t2m    = squeeze(xa.t2m(:,:,ii));

      if iWgt < 0
        a2.siconc = 0*squeeze(xa.siconc(:,:,ii));
        a2.sp     = 0*squeeze(xa.sp(:,:,ii));
        a2.tcc    = 0*squeeze(xa.tcc(:,:,ii));
        a2.u10    = 0*squeeze(xa.u10(:,:,ii));
        a2.v10    = 0*squeeze(xa.v10(:,:,ii));
        a2.skt    = 0*squeeze(xa.skt(:,:,ii));
        a2.d2m    = 0*squeeze(xa.d2m(:,:,ii));
        a2.t2m    = 0*squeeze(xa.t2m(:,:,ii));
      else
        if ii <= 7
          a2.siconc = squeeze(xa.siconc(:,:,ii+1));
          a2.sp     = squeeze(xa.sp(:,:,ii+1));
          a2.tcc    = squeeze(xa.tcc(:,:,ii+1));
          a2.u10    = squeeze(xa.u10(:,:,ii+1));
          a2.v10    = squeeze(xa.v10(:,:,ii+1));
          a2.skt    = squeeze(xa.skt(:,:,ii+1));
          a2.d2m    = squeeze(xa.d2m(:,:,ii+1));
          a2.t2m    = squeeze(xa.t2m(:,:,ii+1));
        elseif ii == 8
          a2.siconc = squeeze(xa.siconc(:,:,1));
          a2.sp     = squeeze(xa.sp(:,:,1));
          a2.tcc    = squeeze(xa.tcc(:,:,1));
          a2.u10    = squeeze(xa.u10(:,:,1));
          a2.v10    = squeeze(xa.v10(:,:,1));
          a2.skt    = squeeze(xa.skt(:,:,1));
          a2.d2m    = squeeze(xa.d2m(:,:,1));
          a2.t2m    = squeeze(xa.t2m(:,:,1));
        end
      end
  
      nanaA = find(indA == ii);
      nanaD = find(indD == ii);
      
      wspd1 = sqrt(a1.u10.^2 + a1.v10.^2);
      wspd2 = sqrt(a2.u10.^2 + a2.v10.^2);
      F1 = griddedInterpolant(iX,iY,flipud(wspd1'));  
      F2 = griddedInterpolant(iX,iY,flipud(wspd2'));  
      junk1 = F1(p0.rlat,wrapTo360(p0.rlon));
      junk2 = F2(p0.rlat,wrapTo360(p0.rlon));
      wspeedsaveA(iCnt,nanaA) = wgtA(nanaA).*junk1(nanaA) + (1-wgtA(nanaA)).*junk2(nanaA);
      wspeedsaveD(iCnt,nanaD) = wgtD(nanaD).*junk1(nanaD) + (1-wgtD(nanaD)).*junk2(nanaD);

      F1 = griddedInterpolant(iX,iY,flipud(a1.siconc'));
      F2 = griddedInterpolant(iX,iY,flipud(a2.siconc'));
      junk1 = F1(p0.rlat,wrapTo360(p0.rlon));
      junk2 = F2(p0.rlat,wrapTo360(p0.rlon));
      sisaveA(iCnt,nanaA) = wgtA(nanaA).*junk1(nanaA) + (1-wgtA(nanaA)).*junk2(nanaA);
      sisaveD(iCnt,nanaD) = wgtD(nanaD).*junk1(nanaD) + (1-wgtD(nanaD)).*junk2(nanaD);
      
      F1 = griddedInterpolant(iX,iY,flipud(a1.sp'));
      F2 = griddedInterpolant(iX,iY,flipud(a2.sp'));
      junk1 = F1(p0.rlat,wrapTo360(p0.rlon));
      junk2 = F2(p0.rlat,wrapTo360(p0.rlon));
      spsaveA(iCnt,nanaA) = wgtA(nanaA).*junk1(nanaA) + (1-wgtA(nanaA)).*junk2(nanaA);
      spsaveD(iCnt,nanaD) = wgtD(nanaD).*junk1(nanaD) + (1-wgtD(nanaD)).*junk2(nanaD);
      
      F1 = griddedInterpolant(iX,iY,flipud(a1.tcc'));
      F2 = griddedInterpolant(iX,iY,flipud(a2.tcc'));
      junk1 = F1(p0.rlat,wrapTo360(p0.rlon));
      junk2 = F2(p0.rlat,wrapTo360(p0.rlon));
      tccsaveA(iCnt,nanaA) = wgtA(nanaA).*junk1(nanaA) + (1-wgtA(nanaA)).*junk2(nanaA);
      tccsaveD(iCnt,nanaD) = wgtD(nanaD).*junk1(nanaD) + (1-wgtD(nanaD)).*junk2(nanaD);
      
      F1 = griddedInterpolant(iX,iY,flipud(a1.skt'));
      F2 = griddedInterpolant(iX,iY,flipud(a2.skt'));
      junk1 = F1(p0.rlat,wrapTo360(p0.rlon));
      junk2 = F2(p0.rlat,wrapTo360(p0.rlon));
      sktsaveA(iCnt,nanaA) = wgtA(nanaA).*junk1(nanaA) + (1-wgtA(nanaA)).*junk2(nanaA);
      sktsaveD(iCnt,nanaD) = wgtD(nanaD).*junk1(nanaD) + (1-wgtD(nanaD)).*junk2(nanaD);
      
      F1 = griddedInterpolant(iX,iY,flipud(a1.d2m'));
      F2 = griddedInterpolant(iX,iY,flipud(a2.d2m'));
      junk1 = F1(p0.rlat,wrapTo360(p0.rlon));
      junk2 = F2(p0.rlat,wrapTo360(p0.rlon));
      d2msaveA(iCnt,nanaA) = wgtA(nanaA).*junk1(nanaA) + (1-wgtA(nanaA)).*junk2(nanaA);
      d2msaveD(iCnt,nanaD) = wgtD(nanaD).*junk1(nanaD) + (1-wgtD(nanaD)).*junk2(nanaD);
      
      F1 = griddedInterpolant(iX,iY,flipud(a1.t2m'));
      F2 = griddedInterpolant(iX,iY,flipud(a2.t2m'));
      junk1 = F1(p0.rlat,wrapTo360(p0.rlon));
      junk2 = F2(p0.rlat,wrapTo360(p0.rlon));
      t2msaveA(iCnt,nanaA) = wgtA(nanaA).*junk1(nanaA) + (1-wgtA(nanaA)).*junk2(nanaA);
      t2msaveD(iCnt,nanaD) = wgtD(nanaD).*junk1(nanaD) + (1-wgtD(nanaD)).*junk2(nanaD);
      
    end         
    %pause(0.1)

  end
end

%% eair2 = 6.1079 *exp(17.269.*(TdewSurf0-273)./(237.3+(TdewSurf0-273)));  %% Eqn 7 of paper
e2msaveD = 6.1079 *exp(17.269.*(t2msaveD-273)./(237.3+(t2msaveD-273)));
e2msaveA = 6.1079 *exp(17.269.*(t2msaveA-273)./(237.3+(t2msaveA-273)));

for iCnt = 1 : 240
  daysSince2002(iCnt) = change2days(yysave(iCnt),mmsave(iCnt),15,2002);
end

figure(7); clf; scatter_coast(p0.rlon,p0.rlat,50,nanmean(wspeedsaveD-wspeedsaveA,1)); caxis([0 15]); colormap jet; title('D-A <speed> over 20 years')
pause(0.1);

%{
if iWgt < 0
  save scalar_2002_09_2022_08.mat *save* 
else
  save scalar_2002_09_2022_08_wgtT1T2.mat *save* 
end
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('doing 4608 treds for D and A ... "+" = 1000, "." = 100')
warning off
for ii = 1 : 4608
  if mod(ii,1000) == 0
    fprintf(1,'+');
  elseif mod(ii,100) == 0
    fprintf(1,'.');
  end

  [B, err, stats] = Math_tsfit_lin_robust(daysSince2002,wspeedsaveA(:,ii),4);
  wspeedtrendA(ii) = B(2);
  wspeedtrendA_err(ii) = err.se(2);

  [B, err, stats] = Math_tsfit_lin_robust(daysSince2002,d2msaveA(:,ii),4);
  d2mtrendA(ii) = B(2);
  d2mtrendA_err(ii) = err.se(2);

  [B, err, stats] = Math_tsfit_lin_robust(daysSince2002,t2msaveA(:,ii),4);
  t2mtrendA(ii) = B(2);
  t2mtrendA_err(ii) = err.se(2);

  woom = sisaveA(:,ii);
  boo = find(isfinite(woom));
  if length(boo) > 20
    [B, err, stats] = Math_tsfit_lin_robust(daysSince2002(boo),sisaveA(boo,ii),4);
    sitrendA(ii) = B(2);
    sitrendA_err(ii) = err.se(2);
  else
    sitrendA(ii) = NaN;
    sitrendA_err(ii) = NaN;
  end

  [B, err, stats] = Math_tsfit_lin_robust(daysSince2002,sktsaveA(:,ii),4);
  skttrendA(ii) = B(2);
  skttrendA_err(ii) = err.se(2);

  [B, err, stats] = Math_tsfit_lin_robust(daysSince2002,spsaveA(:,ii),4);
  sptrendA(ii) = B(2);
  sptrendA_err(ii) = err.se(2);

  [B, err, stats] = Math_tsfit_lin_robust(daysSince2002,tccsaveA(:,ii),4);
  tcctrendA(ii) = B(2);
  tcctrendA_err(ii) = err.se(2);

  [B, err, stats] = Math_tsfit_lin_robust(daysSince2002,e2msaveA(:,ii),4);
  e2mtrendA(ii) = B(2);
  e2mtrendA_err(ii) = err.se(2);

  %%%%%%%%%%%%%%%%%%%%%%%%%

  [B, err, stats] = Math_tsfit_lin_robust(daysSince2002,wspeedsaveD(:,ii),4);
  wspeedtrendD(ii) = B(2);
  wspeedtrendD_err(ii) = err.se(2);

  [B, err, stats] = Math_tsfit_lin_robust(daysSince2002,d2msaveD(:,ii),4);
  d2mtrendD(ii) = B(2);
  d2mtrendD_err(ii) = err.se(2);

  [B, err, stats] = Math_tsfit_lin_robust(daysSince2002,t2msaveD(:,ii),4);
  t2mtrendD(ii) = B(2);
  t2mtrendD_err(ii) = err.se(2);

  woom = sisaveD(:,ii);
  boo = find(isfinite(woom));
  if length(boo) > 20
    [B, err, stats] = Math_tsfit_lin_robust(daysSince2002(boo),sisaveD(boo,ii),4);
    sitrendD(ii) = B(2);
    sitrendD_err(ii) = err.se(2);
  else
    sitrendD(ii) = NaN;
    sitrendD_err(ii) = NaN;
  end

  [B, err, stats] = Math_tsfit_lin_robust(daysSince2002,sktsaveD(:,ii),4);
  skttrendD(ii) = B(2);
  skttrendD_err(ii) = err.se(2);

  [B, err, stats] = Math_tsfit_lin_robust(daysSince2002,spsaveD(:,ii),4);
  sptrendD(ii) = B(2);
  sptrendD_err(ii) = err.se(2);

  [B, err, stats] = Math_tsfit_lin_robust(daysSince2002,tccsaveD(:,ii),4);
  tcctrendD(ii) = B(2);
  tcctrendD_err(ii) = err.se(2);

  [B, err, stats] = Math_tsfit_lin_robust(daysSince2002,e2msaveD(:,ii),4);
  e2mtrendD(ii) = B(2);
  e2mtrendD_err(ii) = err.se(2);

end
warning on
fprintf(1,'\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath /asl/rtp_prod2/emis/
[p0.salti,p0.landfrac] = usgs_deg10_dem(p0.rlat, p0.rlon);
%ocean = find(p0.landfrac == 0); whos ocean
%plot(-0.1:0.01:+0.1,histc(trend(ocean),-0.1:0.01:+0.1)); grid
%junk = [nanmean(trend(ocean)) nanstd(trend(ocean)) nanmean(abs(trend(ocean))) max(abs(trend(ocean)))];
%fprintf(1,'ocean windspeed trends : %8.4f +/- %8.4f m/s; mean(abs(trend)) = %8.4f max(abs(trend)) = %8.4f \n',junk); 

%{
if iWgt < 0
  save scalar_2002_09_2022_08.mat *save* *trend*
else
  save scalar_2002_09_2022_08_wgtT1T2.mat *save* *trend*
end
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(8); clf; scatter_coast(p0.rlon,p0.rlat,50,e2mtrendD); caxis([-1 +1]/5);     colormap(llsmap5); title('D e2m trend over 20 years')
figure(7); clf; scatter_coast(p0.rlon,p0.rlat,50,wspeedtrendD); caxis([-1 +1]/10); colormap(llsmap5); title('D wspeed trend over 20 years')
figure(1); clf; scatter_coast(p0.rlon,p0.rlat,50,sitrendD); caxis([-1 +1]*0.003);  colormap(llsmap5); title('D si trend over 20 years')
figure(2); clf; scatter_coast(p0.rlon,p0.rlat,50,sptrendD); caxis([-1 +1]*0.25);   colormap(llsmap5); title('D sp trend over 20 years')
figure(3); clf; scatter_coast(p0.rlon,p0.rlat,50,tcctrendD); caxis([-1 +1]*0.005); colormap(llsmap5); title('D tcc trend over 20 years')
figure(4); clf; scatter_coast(p0.rlon,p0.rlat,50,d2mtrendD); caxis([-1 +1]*0.15);  colormap(llsmap5); title('D d2m trend over 20 years')
figure(5); clf; scatter_coast(p0.rlon,p0.rlat,50,t2mtrendD); caxis([-1 +1]*0.15);  colormap(llsmap5); title('D t2m trend over 20 years')
figure(6); clf; scatter_coast(p0.rlon,p0.rlat,50,skttrendD); caxis([-1 +1]*0.15);  colormap(llsmap5); title('D skt trend over 20 years')

figure(8); clf; scatter_coast(p0.rlon,p0.rlat,50,e2mtrendA); caxis([-1 +1]/5);     colormap(llsmap5); title('A e2m trend over 20 years')
figure(7); clf; scatter_coast(p0.rlon,p0.rlat,50,wspeedtrendA); caxis([-1 +1]/10); colormap(llsmap5); title('A wspeed trend over 20 years')
figure(1); clf; scatter_coast(p0.rlon,p0.rlat,50,sitrendA); caxis([-1 +1]*0.003);  colormap(llsmap5); title('A si trend over 20 years')
figure(2); clf; scatter_coast(p0.rlon,p0.rlat,50,sptrendA); caxis([-1 +1]*0.25);   colormap(llsmap5); title('A sp trend over 20 years')
figure(3); clf; scatter_coast(p0.rlon,p0.rlat,50,tcctrendA); caxis([-1 +1]*0.005); colormap(llsmap5); title('A tcc trend over 20 years')
figure(4); clf; scatter_coast(p0.rlon,p0.rlat,50,d2mtrendA); caxis([-1 +1]*0.15);  colormap(llsmap5); title('A d2m trend over 20 years')
figure(5); clf; scatter_coast(p0.rlon,p0.rlat,50,t2mtrendA); caxis([-1 +1]*0.15);  colormap(llsmap5); title('A t2m trend over 20 years')
figure(6); clf; scatter_coast(p0.rlon,p0.rlat,50,skttrendA); caxis([-1 +1]*0.15);  colormap(llsmap5); title('A skt trend over 20 years')

figure(1); clf; scatter_coast(p0.rlon,p0.rlat,50,tcctrendD); caxis([-1 +1]*0.005); colormap(llsmap5); title('D tcc trend over 20 years')
figure(2); clf; scatter_coast(p0.rlon,p0.rlat,50,tcctrendA); caxis([-1 +1]*0.005); colormap(llsmap5); title('A tcc trend over 20 years')
figure(3); clf; scatter_coast(p0.rlon,p0.rlat,50,skttrendD); caxis([-1 +1]*0.15);  colormap(llsmap5); title('D skt trend over 20 years')
figure(4); clf; scatter_coast(p0.rlon,p0.rlat,50,skttrendA); caxis([-1 +1]*0.15);  colormap(llsmap5); title('A skt trend over 20 years')
figure(5); clf; scatter_coast(p0.rlon,p0.rlat,50,e2mtrendD); caxis([-1 +1]/5);     colormap(llsmap5); title('D e2m trend over 20 years')
figure(6); clf; scatter_coast(p0.rlon,p0.rlat,50,e2mtrendA); caxis([-1 +1]/5);     colormap(llsmap5); title('A e2m trend over 20 years')

aslmap(1,rlat65,rlon73,smoothn((reshape(tcctrendD',72,64)') ,1), [-90 +90],[-180 +180]); caxis([-1 +1]*0.005); colormap(llsmap5); title('D tcc trend over 20 years')
aslmap(2,rlat65,rlon73,smoothn((reshape(tcctrendA',72,64)') ,1), [-90 +90],[-180 +180]); caxis([-1 +1]*0.005); colormap(llsmap5); title('A tcc trend over 20 years')
aslmap(3,rlat65,rlon73,smoothn((reshape(skttrendD',72,64)') ,1), [-90 +90],[-180 +180]); caxis([-1 +1]*0.15);  colormap(llsmap5); title('D skt trend over 20 years')
aslmap(4,rlat65,rlon73,smoothn((reshape(skttrendA',72,64)') ,1), [-90 +90],[-180 +180]); caxis([-1 +1]*0.15);  colormap(llsmap5); title('A skt trend over 20 years')
aslmap(5,rlat65,rlon73,smoothn((reshape(e2mtrendD',72,64)') ,1), [-90 +90],[-180 +180]); caxis([-1 +1]/5);     colormap(llsmap5); title('D e2m trend over 20 years')
aslmap(6,rlat65,rlon73,smoothn((reshape(e2mtrendA',72,64)') ,1), [-90 +90],[-180 +180]); caxis([-1 +1]/5);     colormap(llsmap5); title('A e2m trend over 20 years')

iTile = 2788;
fprintf(1,'results for India tile %4i \n',iTile);
fprintf(1,'d/dt E2M  D  A  = %8.4f %8.4f mb/yr \n',e2mtrendD(iTile),e2mtrendA(iTile))
fprintf(1,'d/dt SKT  D  A  = %8.4f %8.4f K /yr \n',skttrendD(iTile),skttrendA(iTile))
fprintf(1,'d/dt T2M  D  A  = %8.4f %8.4f K /yr \n',t2mtrendD(iTile),t2mtrendA(iTile))
fprintf(1,'d/dt D2M  D  A  = %8.4f %8.4f K /yr \n',d2mtrendD(iTile),d2mtrendA(iTile))
fprintf(1,'d/dt TCC  D  A  = %8.4f %8.4f []/yr \n',tcctrendD(iTile),tcctrendA(iTile))

disp('from ERA5 trend fies in /home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/');
mooA = load('/home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/ERA5_atm_N_cld_data_2002_09_to_2022_08_trends_asc_surf.mat');
mooD = load('/home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/ERA5_atm_N_cld_data_2002_09_to_2022_08_trends_desc_surf.mat');
fprintf(1,'d/dt E2M  D  A  = %8.4f %8.4f mb/yr \n',mooD.trend_e2m(iTile),mooA.trend_e2m(iTile))
fprintf(1,'d/dt SKT  D  A  = %8.4f %8.4f K /yr \n',mooD.trend_stemp(iTile),mooA.trend_stemp(iTile))
fprintf(1,'d/dt T2M  D  A  = %8.4f %8.4f K /yr \n',mooD.trend_t2m(iTile),mooA.trend_t2m(iTile))
%fprintf(1,'d/dt D2M  D  A  = %8.4f %8.4f K /yr \n',mooD.trend_t2m(iTile),mooA.trend_t2m(iTile))
fprintf(1,'d/dt colW D  A  = %8.4f %8.4f mm/yr \n',mooD.trend_mmw(iTile),mooA.trend_mmw(iTile))

error(';js;js')
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OLD JUNK
      %wspd = sqrt(a1.u10.^2 + a1.v10.^2);
      %wspd_max  = max(wspd,3);
      %wspd_min  = min(wspd,3);
      %wspd_mean = mean(wspd,3);
      %F = griddedInterpolant(iX,iY,flipud(wspd_mean'));  
      %junk = F(p0.rlat,wrapTo360(p0.rlon));
      %pcolor(a1.longitude,a1.latitude,wspd_mean'); shading flat; colorbar; colormap jet; caxis([0 15])
      %figure(7); clf; scatter_coast(p0.rlon,p0.rlat,50,junk); caxis([0 15]); colormap jet; title('Wspeed')
      %wspeedsaveA(iCnt,nanaA) = junk(nanaA);
      %wspeedsaveD(iCnt,nanaD) = junk(nanaD);

      %F =  griddedInterpolant(iX,iY,flipud(a1.siconc'));
      %junk = F(p0.rlat,wrapTo360(p0.rlon));
      %%figure(1); clf; scatter_coast(p0.rlon,p0.rlat,50,junk); caxis([0 1]); colormap jet; title('SeaIce')
      %sisaveA(iCnt,nanaA) = junk(nanaA);
      %sisaveD(iCnt,nanaD) = junk(nanaD);

      %F =  griddedInterpolant(iX,iY,flipud(a1.sp'));
      %junk = F(p0.rlat,wrapTo360(p0.rlon))/100;
      %spsaveA(iCnt,nanaA) = junk(nanaA);
      %spsaveD(iCnt,nanaD) = junk(nanaD);
      
      %F =  griddedInterpolant(iX,iY,flipud(a1.tcc'));
      %junk = F(p0.rlat,wrapTo360(p0.rlon));
      %tccsaveA(iCnt,nanaA) = junk(nanaA);
      %tccsaveD(iCnt,nanaD) = junk(nanaD);
      
      %F =  griddedInterpolant(iX,iY,flipud(a1.t2m'));
      %junk = F(p0.rlat,wrapTo360(p0.rlon));
      %t2msaveA(iCnt,nanaA) = junk(nanaA);
      %t2msaveD(iCnt,nanaD) = junk(nanaD);
      
      %F =  griddedInterpolant(iX,iY,flipud(a1.d2m'));
      %junk = F(p0.rlat,wrapTo360(p0.rlon));
      %d2msaveA(iCnt,nanaA) = junk(nanaA);
      %d2msaveD(iCnt,nanaD) = junk(nanaD);
      
      %F =  griddedInterpolant(iX,iY,flipud(a1.skt'));
      %junk = F(p0.rlat,wrapTo360(p0.rlon));
      %sktsaveA(iCnt,nanaA) = junk(nanaA);
      %sktsaveD(iCnt,nanaD) = junk(nanaD);

      %figure(7); clf; scatter_coast(p0.rlon,p0.rlat,50,junk); caxis([0 15]); colormap jet; title('Wspeed')
      %figure(1); clf; scatter_coast(p0.rlon,p0.rlat,50,junk); caxis([0 1]); colormap jet; title('SeaIce')
      %figure(2); clf; scatter_coast(p0.rlon,p0.rlat,50,junk); caxis([700 1050]); colormap jet; title('SPRES')
      %figure(3); clf; scatter_coast(p0.rlon,p0.rlat,50,junk); caxis([0 1]); colormap jet; title('TCC')
      %figure(4); clf; scatter_coast(p0.rlon,p0.rlat,50,junk); caxis([200 320]); colormap jet; title('T2M')
      %figure(5); clf; scatter_coast(p0.rlon,p0.rlat,50,junk); caxis([200 320]); colormap jet; title('D2M')
      %figure(6); clf; scatter_coast(p0.rlon,p0.rlat,50,junk); caxis([200 320]); colormap jet; title('SKT')
      
      %%%%%%%%%%%%%%%%%%%%%%%%%
