function enames = get_era5_monthly_enames(mtime,iType);

iDebugINCOMING = -1;
iDebugINCOMING = +1;

if nargin == 1
  iType = 1;
end

if iType == 0
  mtime1 = mtime(1);
  
  % /asl/models/era5_avg/2002/2002-08_lev.nc
  % /asl/models/era5_avg/2002/2002-08_sfc.nc
  enames.fn  = ['/asl/models/era5_avg/' datestr(mtime1,'yyyy') '/' datestr(mtime1,'yyyy') '-' datestr(mtime1,'mm') ];
  enames.sfc = ['/asl/models/era5_avg/' datestr(mtime1,'yyyy') '/' datestr(mtime1,'yyyy') '-' datestr(mtime1,'mm') '_sfc.nc'];
  enames.lev = ['/asl/models/era5_avg/' datestr(mtime1,'yyyy') '/' datestr(mtime1,'yyyy') '-' datestr(mtime1,'mm') '_lev.nc'];
  enames.twm = ['/asl/models/era5_avg/' datestr(mtime1,'yyyy') '/' datestr(mtime1,'yyyy') '-' datestr(mtime1,'mm') '_2meter.nc'];
 
  if iDebugINCOMING > 0
    enames.fn  = ['/asl/models/era5_avg/INCOMING/'  datestr(mtime1,'yyyy') '-' datestr(mtime1,'mm') ];
    enames.sfc = ['/asl/models/era5_avg/INCOMING/'  datestr(mtime1,'yyyy') '-' datestr(mtime1,'mm') '_sfc.nc'];
    enames.lev = ['/asl/models/era5_avg/INCOMING/'  datestr(mtime1,'yyyy') '-' datestr(mtime1,'mm') '_lev.nc'];
    enames.twm = ['/asl/models/era5_avg/INCOMING/'  datestr(mtime1,'yyyy') '-' datestr(mtime1,'mm') '_2meter.nc'];
  end

  %enames
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
  mtime1 = mtime;
  
  clear enames
  % round to get 8 forecast hours per day
  rmtime = round(mtime*8)/8;
  
  timestr = datestr(rmtime,'yyyymmddhh');
  ystr = timestr(:,1:4);
  mstr = timestr(:,5:6);
  dstr = timestr(:,7:8);
  hstr = timestr(:,9:10);
  
  % index for fhrstr
  h = str2num(hstr);
  
  %enames = [mstr dstr mstr dstr hstr];
  %enames = cellstr(xenames);

  for ii = 1 : length(mtime)
    enames.fn{ii}  = ['/asl/models/era5_avg/' datestr(mtime1(ii),'yyyy') '/' datestr(mtime1(ii),'yyyy') '-' datestr(mtime1(ii),'mm') ];  
    enames.sfc{ii} = ['/asl/models/era5_avg/' datestr(mtime1(ii),'yyyy') '/' datestr(mtime1(ii),'yyyy') '-' datestr(mtime1(ii),'mm') '_sfc.nc'];
    enames.lev{ii} = ['/asl/models/era5_avg/' datestr(mtime1(ii),'yyyy') '/' datestr(mtime1(ii),'yyyy') '-' datestr(mtime1(ii),'mm') '_lev.nc'];
    enames.twm{ii} = ['/asl/models/era5_avg/' datestr(mtime1(ii),'yyyy') '/' datestr(mtime1(ii),'yyyy') '-' datestr(mtime1(ii),'mm') '_2meter.nc'];

    if iDebugINCOMING > 0
      enames.fn{ii}  = ['/asl/models/era5_avg/INCOMING/'  datestr(mtime1(ii),'yyyy') '-' datestr(mtime1(ii),'mm') ];
      enames.sfc{ii} = ['/asl/models/era5_avg/INCOMING/'  datestr(mtime1(ii),'yyyy') '-' datestr(mtime1(ii),'mm') '_sfc.nc'];
      enames.lev{ii} = ['/asl/models/era5_avg/INCOMING/'  datestr(mtime1(ii),'yyyy') '-' datestr(mtime1(ii),'mm') '_lev.nc'];
      enames.twm{ii} = ['/asl/models/era5_avg/INCOMING/'  datestr(mtime1(ii),'yyyy') '-' datestr(mtime1(ii),'mm') '_2meter.nc'];
    end

  end
end
