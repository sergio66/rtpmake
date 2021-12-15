function [prof,head] = fill_era5_monthly2m_mat(prof,head,yyM,mmM);

%% based on /umbc/xfs2/strow/asl/s1/sergio/home/MATLABCODE/RTPMAKE/CLUST_RTPMAKE/GRIB/fill_era5_monthly.m

f2m = ['/asl/s1/sergio/ERA5_monthlyavg_2m/monthly2m_ERA5average_' num2str(yyM,'%04d') '_' num2str(mmM,'%02d') '_2meter.mat'];
x = load(f2m);

addpath /asl/matlib/aslutil
addpath /asl/packages/time

ename = '';  % This should be placed outside a rtp file loop
mtime = tai2dnum(prof.rtime);

% Get a cell array of era grib files for each time
% Round to get 4 forecast hours per day

rmtime = round(mtime*8)/8;
timestr = datestr(rmtime,'yyyymmddhh');
ystr = timestr(:,1:4);
mstr = timestr(:,5:6);
dstr = timestr(:,7:8);
hstr = timestr(:,9:10);
yearindex = str2num(ystr);
dayindex = str2num(dstr);
hourindex = str2num(hstr);
  figure(1); hist(hourindex,25);

F(1) = grib_interpolate_era5_monthly2m_mat(x,1);
F(2) = grib_interpolate_era5_monthly2m_mat(x,2);
F(3) = grib_interpolate_era5_monthly2m_mat(x,3);
F(4) = grib_interpolate_era5_monthly2m_mat(x,4);
F(5) = grib_interpolate_era5_monthly2m_mat(x,5);
F(6) = grib_interpolate_era5_monthly2m_mat(x,6);
F(7) = grib_interpolate_era5_monthly2m_mat(x,7);
F(8) = grib_interpolate_era5_monthly2m_mat(x,8);

ic = ones(size(mtime));
i = 1;
m = find( ic == i );  % indices of first era file
u_hour = unique(hourindex);
nn = length(u_hour);
% Only loop over hours needed
for jj = 1:nn
  % index for this hour (1:8);  u_hour = [0 3 6 9 12 15 18 21]
  fhi = (u_hour(jj)/3) + 1;
  l = find( hourindex == u_hour(jj));
  k = intersect(l,m);

  if k > 0
    % Assume rtp lat/lon are +-180??  Need to be 0-360 for grib interpolation
    rlat = prof.rlat(k);
    rlon = prof.rlon(k);
    rlon(rlon<0) = rlon(rlon<0) + 360;

    try
      %% 2m air and dewpoint temperatures
      prof.d2m(k)    = F(fhi).d2m.ig(rlat,rlon);
      prof.t2m(k)    = F(fhi).t2m.ig(rlat,rlon);
      prof.skt_2m(k) = F(fhi).skt_2m.ig(rlat,rlon); %% bizarre name!!!
    end
  end
end
