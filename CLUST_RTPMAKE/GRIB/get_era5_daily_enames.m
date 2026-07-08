function enames = get_era5_daily_enames(mtime,fhdr);

mtime1 = mtime(1);

if nargin == 1
  fhdr = '/asl/data/era5/';            %% from Oct 2019 -
  fhdr = '/asl/models/era5/';          %% from Oct 2019 -
  fhdr = '/umbc/rs/strow/asl/era5/';   %% Sad Nov2025-Jan2026
  fhdr = '/umbc/rs/strow/asl/ERA5/';   %% A new Beginning
end

% ls /asl/models/era5/2014/01
% 20140101_lev_test.nc  20140105_lev_test.nc  20140109_lev_test.nc  20140113_lev_test.nc  20140117_lev_test.nc  20140121_lev_test.nc  20140125_lev_test.nc  20140129_lev_test.nc
% 20140101_sfc.nc       20140105_sfc.nc       20140109_sfc.nc       20140113_sfc.nc       20140117_sfc.nc       20140121_sfc.nc       20140125_sfc.nc       20140129_sfc.nc
% 20140102_lev_test.nc  20140106_lev_test.nc  20140110_lev_test.nc  20140114_lev_test.nc  20140118_lev_test.nc  20140122_lev_test.nc  20140126_lev_test.nc  20140130_lev_test.nc

%% steve
enames.fn  = [fhdr datestr(mtime1,'yyyy') '/' datestr(mtime1,'mm') '/' datestr(mtime1,'yyyy') datestr(mtime1,'mm') datestr(mtime1,'dd') ];
enames.sfc = [fhdr datestr(mtime1,'yyyy') '/' datestr(mtime1,'mm') '/' datestr(mtime1,'yyyy') datestr(mtime1,'mm') datestr(mtime1,'dd') '_sfc.nc'];
enames.lev = [fhdr datestr(mtime1,'yyyy') '/' datestr(mtime1,'mm') '/' datestr(mtime1,'yyyy') datestr(mtime1,'mm') datestr(mtime1,'dd') '_lev_test.nc'];
enames.twm = [fhdr datestr(mtime1,'yyyy') '/' datestr(mtime1,'mm') '/' datestr(mtime1,'yyyy') datestr(mtime1,'mm') datestr(mtime1,'dd') '_2meter.nc'];

%% sergio
%-rw-rw-r-- 1 sergio pi_sergio 20482411256 Jul  2 12:13 /umbc/rs/strow/asl/ERA5/2024/11/era5_ml137_20241114.nc
%-rw-rw-r-- 1 sergio pi_sergio   118654061 Jul  2 00:39 /umbc/rs/strow/asl/ERA5/2024/11/era5_sfc_20241114.nc
enames.fn  = [fhdr datestr(mtime1,'yyyy') '/' datestr(mtime1,'mm') '/era5_' datestr(mtime1,'yyyy') datestr(mtime1,'mm') datestr(mtime1,'dd') ];
enames.sfc = [fhdr datestr(mtime1,'yyyy') '/' datestr(mtime1,'mm') '/era5_sfc_' datestr(mtime1,'yyyy') datestr(mtime1,'mm') datestr(mtime1,'dd') '.nc'];
enames.lev = [fhdr datestr(mtime1,'yyyy') '/' datestr(mtime1,'mm') '/era5_ml137_' datestr(mtime1,'yyyy') datestr(mtime1,'mm') datestr(mtime1,'dd') '.nc'];

%enames
