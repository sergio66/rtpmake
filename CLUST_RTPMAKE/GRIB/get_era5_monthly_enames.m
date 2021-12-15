function enames = get_era5_monthly_enames(mtime);

mtime1 = mtime(1);

% /asl/models/era5_avg/2002/2002-08_lev.nc
% /asl/models/era5_avg/2002/2002-08_sfc.nc
enames.fn  = ['/asl/models/era5_avg/' datestr(mtime1,'yyyy') '/' datestr(mtime1,'yyyy') '-' datestr(mtime1,'mm') ];
enames.sfc = ['/asl/models/era5_avg/' datestr(mtime1,'yyyy') '/' datestr(mtime1,'yyyy') '-' datestr(mtime1,'mm') '_sfc.nc'];
enames.lev = ['/asl/models/era5_avg/' datestr(mtime1,'yyyy') '/' datestr(mtime1,'yyyy') '-' datestr(mtime1,'mm') '_lev.nc'];
enames.twm = ['/asl/models/era5_avg/' datestr(mtime1,'yyyy') '/' datestr(mtime1,'yyyy') '-' datestr(mtime1,'mm') '_2meter.nc'];

%enames
