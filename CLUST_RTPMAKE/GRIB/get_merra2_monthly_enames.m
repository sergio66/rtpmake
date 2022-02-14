function enames = get_merra2_monthly_enames(mtime);

mtime1 = mtime(1);

% /asl/models/merra2_monthly/2002/2002-08_lev.nc
% /asl/models/merra2_monthly/2002/2002-08_sfc.nc
enames.fn  = ['/asl/models/merra2_monthly/' datestr(mtime1,'yyyy') '/merra2_' datestr(mtime1,'yyyy') datestr(mtime1,'mm') ];
enames.sfc = ['/asl/models/merra2_monthly/' datestr(mtime1,'yyyy') '/merra2_' datestr(mtime1,'yyyy') datestr(mtime1,'mm') '_sfc.nc'];
enames.lev = ['/asl/models/merra2_monthly/' datestr(mtime1,'yyyy') '/merra2_' datestr(mtime1,'yyyy') datestr(mtime1,'mm') '_lev.nc'];

%enames.twm = ['/asl/models/merra2_monthly/' datestr(mtime1,'yyyy') '/merra2_' datestr(mtime1,'yyyy') datestr(mtime1,'mm') '_2meter.nc'];
enames.twm = ['/asl/models/merra2_monthly/' datestr(mtime1,'yyyy') '/merra2_' datestr(mtime1,'yyyy') datestr(mtime1,'mm') '_sfc.nc'];

%enames
