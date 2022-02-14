addpath /asl/matlib/rtptools/
addpath /asl/matlib/aslutil
addpath /asl/matlib/h4tools
addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/TIME
addpath /home/sergio/MATLABCODE/PLOTTER
addpath ../GRIB
addpath /home/sergio/MATLABCODE/matlib/clouds/sarta
addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS

system_slurm_stats

[h,ha,p,pa] = rtpread('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/summary_17years_all_lat_all_lon_2002_2019_palts_startSept2002_CLEAR.rtp');

%aE = read_netcdf_lls('/asl/models/era5/2010/11/20101110_sfc.nc');
aE = read_netcdf_lls('/asl/models/era5_avg/2010/2010-11_sfc.nc')
aM = read_netcdf_lls('/asl/models/merra2_monthly/2010/merra2_201011_sfc.nc')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1); plot(p.rlon,p.rlat,'.'); title('rtp');

[X,Y] = ndgrid(aE.latitude,aE.longitude); iXE = flipud(X); iYE = flipud(Y);
figure(2); plot(iYE(:),iXE(:),'.'); title('ERA5')

[X,Y] = ndgrid(aM.latitude,aM.longitude); iXM = flipud(X); iYM = flipud(Y);
figure(3); plot(iYM(:),iXM(:),'.'); title('MERA2')

disp('ret to continue'); pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1); scatter_coast(p.rlon,p.rlat,50,p.stemp); title('rtp');

[X,Y] = ndgrid(aE.latitude,aE.longitude); iXE = flipud(X); iYE = flipud(Y);
junk = squeeze(aE.skt(:,:,1));
junk = junk';
junk = flipud(junk);
whos X Y junk
figure(2); scatter_coast(iYE(:),iXE(:),10,junk(:)); title('ERA5')

[X,Y] = ndgrid(aM.latitude,aM.longitude); iXM = flipud(X); iYM = flipud(Y);
junk = squeeze(aM.skt(:,:,1));
junk = junk';
junk = flipud(junk);
whos X Y junk
figure(3); scatter_coast(iYM(:),iXM(:),10,junk(:)); title('MERRA2')

