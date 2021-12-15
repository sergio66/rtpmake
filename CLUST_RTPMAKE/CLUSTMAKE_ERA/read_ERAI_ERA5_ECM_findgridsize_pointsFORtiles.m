addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/TIME
addpath /home/sergio/MATLABCODE/PLOTTER
addpath /home/sergio/MATLABCODE/COLORMAP
addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS
addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS/IDL_WV_ROUTINES/atmos_phys/MATLAB/
addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS/Strow_humidity/convert_humidity/

addpath /asl/matlib/aslutil
addpath /asl/packages/time
addpath /asl/matlib/science

%% see ~/MATLABCODE/RTPMAKE/CLUST_RTPMAKE/GRIB/fill_era_interp.m

clear all

iERAorECM = input('Enter (5) ERA5 (+1) ERA (-1) ECM : ');
if length(iERAorECM) == 0
  iERAorECM = +1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if iERAorECM == 5
  disp('reading era5')
  
  f1 = '/asl/models/era5_avg/2012/2012-07_sfc.nc';
  f2 = '/asl/models/era5_avg/2012/2012-07_lev.nc';
  
  s1 = read_netcdf_lls(f1);
  s2 = read_netcdf_lls(f2);
  
  yres = 0.25;
  [mm,nn,oo] = size(s1.skt); fprintf(1,'ERA5 size(skt) = %4i %4i %4i \n',mm,nn,oo);
  [(360-0.25)/(mm-1) 180/(nn-1)]    %% remember N/S is unique hence (180)/(nn-1), while E-W wraps around hence (360-0.25)/(mm-1)
  era5_rlon = (1:1440); era5_rlon = era5_rlon-1; era5_rlon = 0 + era5_rlon*0.25;
  era5_rlat = (1:721); era5_rlat = era5_rlat-1; era5_rlat = -90 + era5_rlat*0.25;
  [era5_Y,era5_X] = meshgrid(era5_rlat,era5_rlon);
  
  [salti, landfrac] = usgs_deg10_dem(era5_Y,era5_X);
  
  figure(1); clf; 
  figure(1); scatter_coast(era5_X(:),era5_Y(:),10,landfrac(:));
  for ii = 1 : 8
    figure(1)
      skt = squeeze(s1.skt(:,:,ii)); skt = fliplr(skt); scatter_coast(era5_X(:),era5_Y(:),10,skt(:)); title(['stemp ERA5 ' num2str(ii) ' = ' num2str((ii-1)*3) ' GMT']); pause(0.1);
    figure(2)
      tcc = squeeze(s1.tcc(:,:,ii)); tcc = fliplr(tcc); scatter_coast(era5_X(:),era5_Y(:),10,tcc(:)); title(['tcc ERA5 ' num2str(ii) ' = ' num2str((ii-1)*3) ' GMT']); pause(0.1);
  end
  
  tccall = s1.tcc(:); histc(tccall,100); dn = 0:0.002:1; dtccall = histc(tccall,dn); 
  figure(3); plot(dn,dtccall/length(tccall),'.-'); 
  figure(4); plot(dn,cumsum(dtccall)/length(tccall),'.-')
  
  ocean = find(landfrac == 0);
  otccall = [];
  for ii = 1 : 8;
    junk = squeeze(s1.tcc(:,:,ii)); junk = junk(ocean); junk = junk(:); otccall = [otccall junk];
  end
  histc(otccall,100); dotccall = histc(otccall,dn); 
  figure(3); plot(dn,dtccall/length(tccall),'k.-',dn,dotccall/length(otccall),'b.-'); 
  figure(4); plot(dn,cumsum(dtccall)/length(tccall),'k.-',dn,cumsum(dotccall)/length(otccall),'b.-')
  
  land = find(landfrac > 0);
  ltccall = [];
  for ii = 1 : 8;
    junk = squeeze(s1.tcc(:,:,ii)); junk = junk(land); junk = junk(:); ltccall = [ltccall junk];
  end
  histc(ltccall,100); dltccall = histc(ltccall,dn); 
  figure(3); plot(dn,dtccall/length(tccall),'k.-',dn,dotccall/length(otccall),'b.-',dn,dltccall/length(ltccall),'r.-'); 
    hl = legend('all','ocean','land','location','best'); title('ERA5 hist tcc')
  figure(4); plot(dn,cumsum(dtccall)/length(tccall),'k.-',dn,cumsum(dotccall)/length(otccall),'b.-',dn,cumsum(dltccall)/length(ltccall),'r.-')
    xlim([0 0.1]); grid; hl = legend('all','ocean','land','location','best'); title('ERA5 cumusm tcc')
  
  disp('ret'); pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif iERAorECM == -1
  disp('reading ecmwf')
  
  f1 = '/asl/models/ecmwf/2021/07/UAD07040000070400001-1.nc';
  f2 = '/asl/models/ecmwf/2021/07/UAD07040000070400001-2.nc';
  f3 = '/asl/models/ecmwf/2021/07/UAD07040000070400001-d.nc';
  
  s1 = read_netcdf_lls(f1);
  s2 = read_netcdf_lls(f2);
  s3 = read_netcdf_lls(f3);
  
  yres = 0.25;
  [mm,nn,oo] = size(s1.skt); fprintf(1,'ECM size(skt) = %4i %4i %4i \n',mm,nn,oo);
  [(360-0.25)/(mm-1) 180/(nn-1)]    %% remember N/S is unique hence (180)/(nn-1), while E-W wraps around hence (360-0.25)/(mm-1)
  ecm_rlon = (1:1440); ecm_rlon = ecm_rlon-1; ecm_rlon = 0 + ecm_rlon*0.25;
  ecm_rlat = (1:721); ecm_rlat = ecm_rlat-1; ecm_rlat = -90 + ecm_rlat*0.25;
  [ecm_Y,ecm_X] = meshgrid(ecm_rlat,ecm_rlon);
  
  [salti, landfrac] = usgs_deg10_dem(ecm_Y,ecm_X);
  
  figure(1); clf; 
  figure(1); scatter_coast(ecm_X(:),ecm_Y(:),10,landfrac(:));
  for ii = 1 : 1
    figure(1)
      skt = s1.skt(:,:); skt = fliplr(skt); scatter_coast(ecm_X(:),ecm_Y(:),10,skt(:)); title(['stemp ECM ' num2str(ii) ' = ' num2str((ii-1)*6) ' GMT']); pause(0.1);
    figure(2)
      tcc = s1.tcc(:,:); tcc = fliplr(tcc); scatter_coast(ecm_X(:),ecm_Y(:),10,tcc(:)); title(['tcc ECM ' num2str(ii) ' = ' num2str((ii-1)*6) ' GMT']); pause(0.1);
  end
  
  tccall = s1.tcc(:); histc(tccall,100); dn = 0:0.002:1; dtccall = histc(tccall,dn); 
  figure(3); plot(dn,dtccall/length(tccall),'.-'); 
  figure(4); plot(dn,cumsum(dtccall)/length(tccall),'.-')
  
  ocean = find(landfrac == 0);
  otccall = s1.tcc(ocean); otccall = otccall(:); histc(otccall,100); dotccall = histc(otccall,dn); 
  figure(3); plot(dn,dtccall/length(tccall),'k.-',dn,dotccall/length(otccall),'b.-'); 
  figure(4); plot(dn,cumsum(dtccall)/length(tccall),'k.-',dn,cumsum(dotccall)/length(otccall),'b.-')
  
  land = find(landfrac > 0);
  ltccall = s1.tcc(land); ltccall = ltccall(:); histc(ltccall,100); dltccall = histc(ltccall,dn); 
  figure(3); plot(dn,dtccall/length(tccall),'k.-',dn,dotccall/length(otccall),'b.-',dn,dltccall/length(ltccall),'r.-'); 
    hl = legend('all','ocean','land','location','best'); title('ECM hist tcc')
  figure(4); plot(dn,cumsum(dtccall)/length(tccall),'k.-',dn,cumsum(dotccall)/length(otccall),'b.-',dn,cumsum(dltccall)/length(ltccall),'r.-')
    xlim([0 0.1]); grid; hl = legend('all','ocean','land','location','best'); title('ECM cumusm tcc')
  
  disp('ret'); pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif iERAorECM == +1
  disp('reading era')
  
  f1 = '/asl/models/era/2007/11/20071125_sfc.nc';
  f2 = '/asl/models/era/2007/11/20071125_lev.nc';
  
  s1 = read_netcdf_lls(f1);
  s2 = read_netcdf_lls(f2);
  
  yres = 0.75;
  [mm,nn,oo] = size(s1.skt); fprintf(1,'ERA-I size(skt) = %4i %4i %4i \n',mm,nn,oo);
  [(360-0.75)/(mm-1) 180/(nn-1)]  %% remember N/S is unique hence (180)/(nn-1), while E-W wraps around hence (360-0.75)/(mm-1)
  era_rlon = (1:480); era_rlon = era_rlon-1; era_rlon = 0 + era_rlon*0.75;
  era_rlat = (1:241); era_rlat = era_rlat-1; era_rlat = -90 + era_rlat*0.75;
  [era_Y,era_X] = meshgrid(era_rlat,era_rlon);
  
  [salti, landfrac] = usgs_deg10_dem(era_Y,era_X);
  
  figure(1); clf; 
  figure(1); scatter_coast(era_X(:),era_Y(:),10,landfrac(:));
  for ii = 1 : 4
    figure(1)
      skt = squeeze(s1.skt(:,:,ii)); skt = fliplr(skt); scatter_coast(era_X(:),era_Y(:),10,skt(:)); title(['stemp ERA ' num2str(ii) ' = ' num2str((ii-1)*6) ' GMT']); pause(0.1);
    figure(2)
      tcc = squeeze(s1.tcc(:,:,ii)); tcc = fliplr(tcc); scatter_coast(era_X(:),era_Y(:),10,tcc(:)); title(['tcc ERA ' num2str(ii) ' = ' num2str((ii-1)*6) ' GMT']); pause(0.1);
  end
  
  tccall = s1.tcc(:); histc(tccall,100); dn = 0:0.002:1; dtccall = histc(tccall,dn); 
  figure(3); plot(dn,dtccall/length(tccall),'.-'); 
  figure(4); plot(dn,cumsum(dtccall)/length(tccall),'.-')
  
  ocean = find(landfrac == 0);
  otccall = s1.tcc(ocean); otccall = otccall(:); histc(otccall,100); dotccall = histc(otccall,dn); 
  figure(3); plot(dn,dtccall/length(tccall),'k.-',dn,dotccall/length(otccall),'b.-'); 
  figure(4); plot(dn,cumsum(dtccall)/length(tccall),'k.-',dn,cumsum(dotccall)/length(otccall),'b.-')
  
  land = find(landfrac > 0);
  ltccall = s1.tcc(land); ltccall = ltccall(:); histc(ltccall,100); dltccall = histc(ltccall,dn); 
  figure(3); plot(dn,dtccall/length(tccall),'k.-',dn,dotccall/length(otccall),'b.-',dn,dltccall/length(ltccall),'r.-'); 
    hl = legend('all','ocean','land','location','best'); title('ERA hist tcc')
  figure(4); plot(dn,cumsum(dtccall)/length(tccall),'k.-',dn,cumsum(dotccall)/length(otccall),'b.-',dn,cumsum(dltccall)/length(ltccall),'r.-')
    xlim([0 0.1]); grid; hl = legend('all','ocean','land','location','best'); title('ERA cumusm tcc')
  
  disp('ret'); pause
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('loading in latB64, lon72 for tiles');

load /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_For_HowardObs_TimeSeries/latB64.mat
rlon = -180 : +180;          rlat = -90 : +90;
drlon = 5;                   drlat = 3;

rlon = -180 : drlon : +180;  rlat = -90 : drlat : +90;
rlon = rlon;                 rlat = latB2;

nx = floor((5/yres)+1);  %% longitude is 5 degrees and xres = yres
fprintf(1,'distribution of latbin widths from latB64 .. and how many points will fit there at res %8.4f degrees (num_X_NWP pts = %4i) \n',yres,nx)
  disp('diff(rlatB64)  NumTimes   num_Y_NWPpts num_Y_NWPpts*num_X_NWPpts')
 dx = 0:0.25:9; dn = hist(diff(rlat),dx); junk = [dx; dn; floor(dx/yres + 1); floor(dx/yres + 1)*nx;]'; 
 boo = find(junk(:,2) > 0); junk = junk(boo,:); fprintf(1,'%8.4f         %5i         %5i       %5i \n',junk');
disp(' ')

x = 0.5*(rlon(1:end-1)+rlon(2:end));  y = 0.5*(rlat(1:end-1)+rlat(2:end));

[Y,X] = meshgrid(y,x);
Y = Y(:);
X = X(:);
plot(X(1:73), Y(1:73),'o')   %% yay so I loop     do outer 1 : 64; do inner 1 : 72;    ?????   .....

if iERAorECM == +1
  ii = 4; skt = squeeze(s1.skt(:,:,ii)); skt = fliplr(skt); scatter_coast(era_X(:),era_Y(:),10,skt(:)); title([num2str(ii) ' = ' num2str((ii-1)*6) ' GMT']); pause(0.1);
elseif iERAorECM == 5
  ii = 4; skt = squeeze(s1.skt(:,:,ii)); skt = fliplr(skt); scatter_coast(era5_X(:),era5_Y(:),10,skt(:)); title([num2str(ii) ' = ' num2str((ii-1)*6) ' GMT']); pause(0.1);
elseif iERAorECM == -1
  ii = 4; skt = s1.skt(:,:); skt = fliplr(skt); scatter_coast(ecm_X(:),ecm_Y(:),10,skt(:)); title([num2str(ii) ' = ' num2str((ii-1)*6) ' GMT']); pause(0.1);
end

line([-180 +180],[rlat(23) rlat(23)],'color','k')
line([-180 +180],[rlat(24) rlat(24)],'color','k')
for ii = 1 : 73
  line([rlon(ii) rlon(ii)],[rlat(23) rlat(24)],'color','k')
end

%%%%%%%%%%%%%%%%%%%%%%%%%

if iERAorECM == 1
  nwp_grid_X = era_X;    nwp_grid_Y = era_Y; 
elseif iERAorECM == -1
  nwp_grid_X = ecm_X;    nwp_grid_Y = ecm_Y; 
elseif iERAorECM == 5
  nwp_grid_X = era5_X;    nwp_grid_Y = era5_Y; 
end

figure(2); clf
plot(nwp_grid_X(:)-180,nwp_grid_Y(:),'.');
line([-180 +180],[rlat(23) rlat(23)],'color','k')
line([-180 +180],[rlat(24) rlat(24)],'color','k')
line([-180 +180],[rlat(25) rlat(25)],'color','k')
line([-180 +180],[rlat(26) rlat(26)],'color','k')
for ii = 1 : 73
  line([rlon(ii) rlon(ii)],[rlat(23) rlat(26)],'color','k')
end
axis([-190 +190 rlat(22) rlat(27)]); title('blue dots = era grid; black lines = tiles')
axis([-7.5 +7.5 rlat(22) rlat(27)]); title('blue dots = era grid; black lines = tiles')
axis([-181 -160 rlat(22) rlat(27)]); title('blue dots = era grid; black lines = tiles')

%%%%%%%%%%%%%%%%%%%%%%%%%

figure(2); clf
plot(nwp_grid_X(:)-180,nwp_grid_Y(:),'.');
line([-180 +180],[rlat(01) rlat(01)],'color','k')
line([-180 +180],[rlat(02) rlat(02)],'color','k')
line([-180 +180],[rlat(03) rlat(03)],'color','k')
line([-180 +180],[rlat(04) rlat(04)],'color','k')
for ii = 1 : 73
  line([rlon(ii) rlon(ii)],[rlat(01) rlat(04)],'color','k')
end
axis([-190 +190 rlat(01)-1 rlat(05)]); title('blue dots = era grid; black lines = tiles')
axis([-7.5 +7.5 rlat(01)-1 rlat(05)]); title('blue dots = era grid; black lines = tiles')
axis([-181 -160 rlat(01)-1 rlat(05)]); title('blue dots = era grid; black lines = tiles')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for jj = 1 : length(rlat)-1
  ii = 1;
  oo = find(nwp_grid_X-180 >= rlon(ii)-0.01 & nwp_grid_X-180 <= rlon(ii+1)+0.01 & nwp_grid_Y >= rlat(jj)-0.01 & nwp_grid_Y <= rlat(jj+1)+0.01);
  lenfound(jj) = length(oo);
  mooY = unique(nwp_grid_Y(oo));
  plot(nwp_grid_X(oo)-180,nwp_grid_Y(oo),'o'); title([num2str(jj) ' : length found = ' num2str(length(oo))]); 
  %rectangle('Position',[x0 y0 xlen ylen])
  rectangle('Position',[rlon(ii) rlat(jj) rlon(ii+1)-rlon(ii) rlat(jj+1)-rlat(jj)],'linewidth',2,'edgecolor','r')
  %fprintf(1,' latbin %2i of %2i : unique Y nwp points found %3i <ret to continue> \n',jj,length(rlat)-1,length(mooY)); pause
  fprintf(1,' latbin %2i of %2i : unique Y nwp points found %3i <ret to continue> \n',jj,length(rlat)-1,length(mooY));
  pause(0.1)
end

%% for ERA  0.75 deg ==> 3x5 --> 4 points in lat,  7 points in lon
%% for ECM  0.25 deg ==> 3x5 --> 13 points in lat, 21 points in lon
%% for ERA5 0.25 deg ==> 3x5 --> 13 points in lat, 21 points in lon
for jj = 1 : length(rlat)-1
  for ii = 1 : length(rlon)-1
   %% TRUE NUMBER : latitude : 12 at poles, 3 just off poles, 4 everywehere else
   %% TRUE NUMBER : longitude : 7 everywehre excpet at last tile where it is 6 (because ERA/ECM wraps to 0)
    oo = find(nwp_grid_X-180 >= rlon(ii)-0.01 & nwp_grid_X-180 <= rlon(ii+1)+0.01 & nwp_grid_Y >= rlat(jj)-0.01 & nwp_grid_Y <= rlat(jj+1)+0.01);  
    theNWPgridptslist.true_lenfound(ii,jj) = length(oo);
    theNWPgridptslist.true_indices{ii,jj} = oo;

    oox = unique(nwp_grid_X(oo)-180);
    theNWPgridptslist.adjusted_xlenfound(ii,jj) = length(oox);
    theNWPgridptslist.adjusted_xvalues{ii,jj} = oox;
    if ii == length(rlon)-1
      %% recall things can wrap to 0
      theNWPgridptslist.adjusted_xlenfound(ii,jj) = length(oox) + 1;
      theNWPgridptslist.adjusted_xvalues{ii,jj} = [oox; 0.00];
    end

    if iERAorECM == 1
      ooy = unique(nwp_grid_Y(oo));
      theNWPgridptslist.adjusted_ylenfound(ii,jj) = length(ooy);
      theNWPgridptslist.adjusted_yvalues{ii,jj} = ooy;
      if length(ooy) == 3
        %% some of the polar regions have these problems, so use lower boundary value
        ooy = [rlat(jj); ooy];
        ooy = sort(ooy);
        theNWPgridptslist.adjusted_ylenfound(ii,jj) = length(ooy);
        theNWPgridptslist.adjusted_yvalues{ii,jj} = ooy;
      end
      if length(ooy) == 12
        %% the poles have these problems, so use one every four
        ooy = ooy([1 4 7 10]+1);
        theNWPgridptslist.adjusted_ylenfound(ii,jj) = length(ooy);
        theNWPgridptslist.adjusted_yvalues{ii,jj} = ooy;
      end

    elseif iERAorECM == -1 | iERAorECM == +5
      ooy = unique(nwp_grid_Y(oo));
      theNWPgridptslist.adjusted_ylenfound(ii,jj) = length(ooy);
      theNWPgridptslist.adjusted_yvalues{ii,jj} = ooy;
      if length(ooy) >= 9 & length(ooy) <= 11
        %% some of the polar regions have these problems, so use lower boundary value
        ooyN(1) = 0.5*(ooy(1)+ooy(2));
        ooyN(2) = 0.5*(ooy(2)+ooy(3));
        ooyN(3) = 0.5*(ooy(3)+ooy(4));
        if length(ooy) == 9
          ooy = sort([ooy; ooyN']);
        elseif length(ooy) == 10
          ooy = sort([ooy; ooyN(1:2)']);
        elseif length(ooy) == 11
          ooy = sort([ooy; ooyN(1)]);
        end
        theNWPgridptslist.adjusted_ylenfound(ii,jj) = length(ooy);
        theNWPgridptslist.adjusted_yvalues{ii,jj} = ooy;
      end
      if length(ooy) == 34
        %% the poles have these problems, need 12 so one every 3
        woo = [1:3:34];
        woo = sort(woo);
        ooy = ooy(woo);
        theNWPgridptslist.adjusted_ylenfound(ii,jj) = length(ooy);
        theNWPgridptslist.adjusted_yvalues{ii,jj} = ooy;
      end
    end


  end
end

theNWPgridptslist.adjusted_yvalues{15,1:3}

figure(3); clf;
woo = nan(73,65); xnan = [x nan]; ynan = [y; nan]; woo(1:72,1:64) = theNWPgridptslist.true_lenfound; pcolor(xnan,ynan,woo');
woo = zeros(73,65); xnan = [x x(end)+1]; ynan = [y; y(end)+1]; woo(1:72,1:64) = theNWPgridptslist.true_lenfound; pcolor(woo');

figure(3); clf; pcolor(x,y,theNWPgridptslist.true_lenfound'); colorbar; colormap jet; title('true number of ERA grid points per tile');
figure(3); clf; imagesc(x,y,theNWPgridptslist.true_lenfound'); colorbar; colormap jet; title('true number of ERA grid points per tile');
figure(4); clf; dn = 0:1:100; plot(dn,histc(theNWPgridptslist.true_lenfound(:),dn)); title('true number of ERA grid points per tile');

figure(5); clf; imagesc(x,y,theNWPgridptslist.adjusted_xlenfound'); colorbar; colormap jet; title('adjusted X number of ERA grid points per tile');
figure(6); clf; imagesc(x,y,theNWPgridptslist.adjusted_ylenfound'); colorbar; colormap jet; title('adjusted Y number of ERA grid points per tile');
figure(6); clf; plot(0:15,hist(theNWPgridptslist.adjusted_ylenfound(:),0:15),'o-'); title('Y adjusted number of ERA grid points per tile'); %% so usually 4, sometimes 3, and 12 at the poles
figure(6); clf; plot(theNWPgridptslist.adjusted_ylenfound(15,:),'o-'); xlim([1 64]); title('Y adjusted number of ERA grid points per tile'); %% so usually 4, sometimes 3, and 12 at the poles

comment = 'see /home/sergio/MATLABCODE/RTPMAKE/CLUST_RTPMAKE/CLUSTMAKE_ERA/test_read_eraI.m';
error('popo')

if iERAorECM == 1
  theERAgridptslist = theNWPgridptslist;
  save /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_For_HowardObs_TimeSeries/theERAgridptslist_ERA_tile_points.mat theERAgridptslist comment
elseif iERAorECM == -1
  theECMgridptslist = theNWPgridptslist;
  save /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_For_HowardObs_TimeSeries/theECMgridptslist_ECM_tile_points.mat theECMgridptslist comment
elseif iERAorECM == 5
  theERA5gridptslist = theNWPgridptslist;
  save /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_For_HowardObs_TimeSeries/theERA5gridptslist_ERA5_tile_points.mat theERA5gridptslist comment
end
