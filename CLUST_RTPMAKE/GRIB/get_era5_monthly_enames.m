function enames = get_era5_monthly_enames(mtime,iType,fhdr);

iDebugINCOMING = +1;
iDebugINCOMING = -1;

if nargin == 1
  fhdr = '/asl/models/era5_avg/';         %% from July 2021
  fhdr = '/asl/models/era5_monthly/';     %% from Jan 2024
  iType = 1;
elseif nargin == 2
  fhdr = '/asl/models/era5_avg/';         %% from July 2021
  fhdr = '/asl/models/era5_monthly/';     %% from Jan 2024
end

if iType == 0
  mtime1 = mtime(1);
  
  % /asl/models/era5_avg/2002/2002-08_lev.nc
  % /asl/models/era5_avg/2002/2002-08_sfc.nc
  enames.fn  = [fhdr datestr(mtime1,'yyyy') '/' datestr(mtime1,'yyyy') '-' datestr(mtime1,'mm') ];
  enames.sfc = [fhdr datestr(mtime1,'yyyy') '/' datestr(mtime1,'yyyy') '-' datestr(mtime1,'mm') '_sfc.nc'];
  enames.lev = [fhdr datestr(mtime1,'yyyy') '/' datestr(mtime1,'yyyy') '-' datestr(mtime1,'mm') '_lev.nc'];
  enames.twm = [fhdr datestr(mtime1,'yyyy') '/' datestr(mtime1,'yyyy') '-' datestr(mtime1,'mm') '_2meter.nc'];
 
  if iDebugINCOMING > 0
    enames.fn  = [fhdr '/INCOMING/'  datestr(mtime1,'yyyy') '-' datestr(mtime1,'mm') ];
    enames.sfc = [fhdr '/INCOMING/'  datestr(mtime1,'yyyy') '-' datestr(mtime1,'mm') '_sfc.nc'];
    enames.lev = [fhdr '/INCOMING/'  datestr(mtime1,'yyyy') '-' datestr(mtime1,'mm') '_lev.nc'];
    enames.twm = [fhdr '/INCOMING/'  datestr(mtime1,'yyyy') '-' datestr(mtime1,'mm') '_2meter.nc'];
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
    enames.fn{ii}  = [fhdr datestr(mtime1(ii),'yyyy') '/' datestr(mtime1(ii),'yyyy') '-' datestr(mtime1(ii),'mm') ];  
    enames.sfc{ii} = [fhdr datestr(mtime1(ii),'yyyy') '/' datestr(mtime1(ii),'yyyy') '-' datestr(mtime1(ii),'mm') '_sfc.nc'];
    enames.lev{ii} = [fhdr datestr(mtime1(ii),'yyyy') '/' datestr(mtime1(ii),'yyyy') '-' datestr(mtime1(ii),'mm') '_lev.nc'];
    enames.twm{ii} = [fhdr datestr(mtime1(ii),'yyyy') '/' datestr(mtime1(ii),'yyyy') '-' datestr(mtime1(ii),'mm') '_2meter.nc'];
    enames.olr{ii} = [fhdr datestr(mtime1(ii),'yyyy') '/' datestr(mtime1(ii),'yyyy') '-' datestr(mtime1(ii),'mm') '_rad.nc'];

    if iDebugINCOMING > 0
      enames.fn{ii}  = [fhdr '/INCOMING/'  datestr(mtime1(ii),'yyyy') '-' datestr(mtime1(ii),'mm') ];
      enames.sfc{ii} = [fhdr '/INCOMING/'  datestr(mtime1(ii),'yyyy') '-' datestr(mtime1(ii),'mm') '_sfc.nc'];
      enames.lev{ii} = [fhdr '/INCOMING/'  datestr(mtime1(ii),'yyyy') '-' datestr(mtime1(ii),'mm') '_lev.nc'];
      enames.twm{ii} = [fhdr '/INCOMING/'  datestr(mtime1(ii),'yyyy') '-' datestr(mtime1(ii),'mm') '_2meter.nc'];
      enames.olr{ii} = [fhdr '/INCOMING/'  datestr(mtime1(ii),'yyyy') '-' datestr(mtime1(ii),'mm') '_rad.nc'];
    end

  end
end
