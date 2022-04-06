function [aT,aW,aO3] = mls_reader_L3(fT,fW,fO3,mm);

%{
outout units T, vmr, vmr

fT  = '/asl/xfs3/mls/ML3MBT_004/MLS-Aura_L3MB-Temperature_v04-23-c02_2006.nc';
fW  = '/asl/xfs3/mls/ML3MBH2O_004/MLS-Aura_L3MB-H2O_v04-23-c02_2006.nc';
fO3 = '/asl/xfs3/mls/ML3MBO3_004/MLS-Aura_L3MB-O3_v04-23-c02_2006.nc';
[aT,aW,a03] = mls_reader_L3(fT,fW,fO3);
%% T,O3  works from level 08 to 55 (261 to 0,001 mb)
%% WV works from level 07 to 55 (316 to 0,001 mb)

aO3.O3_PressureGrid.lev ./ aT.Temperature_PressureGrid.lev = 1 = aW.H2O_PressureGrid.lev ./ aT.Temperature_PressureGrid.lev
%}

if nargin == 3
  mm = 8;
end

%aT = read_netcdf_lls_nospace('/asl/s1/sergio/JUNK/MLS-Aura_L3MB-Temperature_v04-23-c02_2004.nc');
aT = read_netcdf_lls_nospace(fT);
plev = aT.Temperature_PressureGrid.lev;
lon  = aT.Temperature_PressureGrid.lon;
lat  = aT.Temperature_PressureGrid.lat;
miaow = squeeze(aT.Temperature_PressureGrid.value(:,:,16,mm))';
  figure(1); clf; pcolor(lon,lat,miaow'); colormap jet; colorbar; shading flat; title('T at level 16 (56 mb)')
miaow = squeeze(aT.Temperature_PressureGrid.value(:,:,08,mm))';
  figure(2); clf; pcolor(lon,lat,miaow'); colormap jet; colorbar; shading flat; title('T at level 08 (261 mb)')

%aW = read_netcdf_lls_nospace('/asl/s1/sergio/JUNK/MLS-Aura_L3MB-H2O_v04-23-c02_2004.nc');
aW = read_netcdf_lls_nospace(fW);
plev = aW.H2O_PressureGrid.lev;
lon  = aW.H2O_PressureGrid.lon;
lat  = aW.H2O_PressureGrid.lat;
miaow = squeeze(aW.H2O_PressureGrid.value(:,:,16,mm))';
  figure(3); clf; pcolor(lon,lat,miaow'); colormap jet; colorbar;  shading flat; title('WV VMR at level 16 (56 mb)')
  figure(3); clf; pcolor(lon,lat,1e6*miaow'); colormap jet; colorbar;  shading flat; title('WV PPMV at level 16 (56 mb)')
miaow = squeeze(aW.H2O_PressureGrid.value(:,:,08,mm))';
  figure(4); clf; pcolor(lon,lat,miaow'); colormap jet; colorbar;  shading flat; title('WV VMR at level 08 (261 mb)')
  figure(4); clf; pcolor(lon,lat,1e6*miaow'); colormap jet; colorbar;  shading flat; title('WV PPMV at level 08 (261 mb)')

%aO3 = read_netcdf_lls_nospace('/asl/s1/sergio/JUNK/MLS-Aura_L3MB-O3_v04-23-c02_2004.nc');
aO3 = read_netcdf_lls_nospace(fO3);
plev = aO3.O3_PressureGrid.lev;
lon  = aO3.O3_PressureGrid.lon;
lat  = aO3.O3_PressureGrid.lat;
miaow = squeeze(aO3.O3_PressureGrid.value(:,:,16,mm))';
  figure(5); clf; pcolor(lon,lat,miaow'); colormap jet; colorbar;  shading flat; title('O3 VMR at level 16 (56 mb)')
  figure(5); clf; pcolor(lon,lat,1e6*miaow'); colormap jet; colorbar;  shading flat; title('O3 PPMV at level 16 (56 mb)')
miaow = squeeze(aO3.O3_PressureGrid.value(:,:,08,mm))';
  figure(6); clf; pcolor(lon,lat,miaow'); colormap jet; colorbar;  shading flat; title('O3 VMR at level 08 (261 mb)')
  figure(6); clf; pcolor(lon,lat,1e6*miaow'); colormap jet; colorbar;  shading flat; title('O3 PPMV at level 08 (261 mb)')
