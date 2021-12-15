function F = grib_interpolate_era5_monthly2m_mat(x,hindex);

%% see MATLABCODE/RTPMAKE/CLUST_RTPMAKE/GRIB/grib_interpolate_era.m
F.s_longitude = x.longitude;
F.s_latitude  = x.latitude;
F.s_time      = x.time;
F.s_mtime     = datenum(1900,0,0,double(F.s_time),0,0);
[X,Y] = ndgrid(F.s_latitude,F.s_longitude);
iX = flipud(X); iY = flipud(Y);

%F.t2m.ig  = griddedInterpolant(iX,iY,flipud(single(ncread(fn_s,'t2m',[1 1 hindex],[Inf Inf 1]))'),'linear');
%F.d2m.ig  = griddedInterpolant(iX,iY,flipud(single(ncread(fn_s,'d2m',[1 1 hindex],[Inf Inf 1]))'),'linear');

F.t2m.ig    = griddedInterpolant(iX,iY,flipud(single(squeeze(x.t2m(:,:,hindex)))'),'linear');
F.d2m.ig    = griddedInterpolant(iX,iY,flipud(single(squeeze(x.d2m(:,:,hindex)))'),'linear');
F.skt_2m.ig = griddedInterpolant(iX,iY,flipud(single(squeeze(x.skt(:,:,hindex)))'),'linear');  %% remember this is my monthly average, so I give it a bizarre name


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

