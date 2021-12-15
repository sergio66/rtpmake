addpath /home/sergio/MATLABCODE/PLOTTER

figure(1);
clear lat
lat = -90 : 5 : 90;
[A,V] = surfacearea_sphericalzones(lat);
fprintf(1,'latbins spaced 5 deg apart : %8.6f \n',sum(A))
latx = 0.5 * (lat(1:end-1) + lat(2:end));
plot(latx,A,'ro-')
plot(latx,A/max(A),'ro-')
title('5 deg latbins')

figure(2)
clear lat
load /home/sergio/MATLABCODE/oem_pkg_run/Sergio_AllSky/lat.mat
[A,V] = surfacearea_sphericalzones(lat);
fprintf(1,'clear sky latbins  : %8.6f \n',sum(A))
latx = 0.5 * (lat(1:end-1) + lat(2:end));
plot(latx,A,'ro-')
plot(latx,A/max(A),'ro-')
title('36 clear sky latbins')

figure(3)
clear lat
N = 20;
lat(1) = 0;
for ii = 2 : N
  sinx = 1/N + sin(lat(ii-1)*pi/180);
  x = asin(sinx)*180/pi;
  lat(ii) = x;
end
lat(end+1) = 90;
lat = [-fliplr(lat(2:end)) lat];
[A,V] = surfacearea_sphericalzones(lat);
fprintf(1,'clear sky latbins  : %8.6f \n',sum(A))
latx = 0.5 * (lat(1:end-1) + lat(2:end));
plot(latx,A,'ro-')
plot(latx,A/max(A),'ro-')
title('Equal Area latbins')