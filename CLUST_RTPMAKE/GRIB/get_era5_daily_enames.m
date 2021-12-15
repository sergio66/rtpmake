function enames = get_era5_daily_enames(mtime);

mtime1 = mtime(1);

% ls /asl/models/era5/2014/01
% 20140101_lev_test.nc  20140105_lev_test.nc  20140109_lev_test.nc  20140113_lev_test.nc  20140117_lev_test.nc  20140121_lev_test.nc  20140125_lev_test.nc  20140129_lev_test.nc
% 20140101_sfc.nc       20140105_sfc.nc       20140109_sfc.nc       20140113_sfc.nc       20140117_sfc.nc       20140121_sfc.nc       20140125_sfc.nc       20140129_sfc.nc
% 20140102_lev_test.nc  20140106_lev_test.nc  20140110_lev_test.nc  20140114_lev_test.nc  20140118_lev_test.nc  20140122_lev_test.nc  20140126_lev_test.nc  20140130_lev_test.nc

enames.fn  = ['/asl/models/era5/' datestr(mtime1,'yyyy') '/' datestr(mtime1,'mm') '/' datestr(mtime1,'yyyy') datestr(mtime1,'mm') datestr(mtime1,'dd') ];
enames.sfc = ['/asl/models/era5/' datestr(mtime1,'yyyy') '/' datestr(mtime1,'mm') '/' datestr(mtime1,'yyyy') datestr(mtime1,'mm') datestr(mtime1,'dd') '_sfc.nc'];
enames.lev = ['/asl/models/era5/' datestr(mtime1,'yyyy') '/' datestr(mtime1,'mm') '/' datestr(mtime1,'yyyy') datestr(mtime1,'mm') datestr(mtime1,'dd') '_lev_test.nc'];
enames.twm = ['/asl/models/era5/' datestr(mtime1,'yyyy') '/' datestr(mtime1,'mm') '/' datestr(mtime1,'yyyy') datestr(mtime1,'mm') datestr(mtime1,'dd') '_2meter.nc'];

%enames
