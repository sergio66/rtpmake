JOB = str2num(getenv('SLURM_ARRAY_TASK_ID'));
%JOB = 25
%JOB = 15

%% one per month, 19 years of AIRS data so 19x12 = 228 sets of data
%JOB = 1 == 2002/01
%JOB = 8 == 2002/08
%------------------
%JOB = 9 == 2002/09 etc

addpath /asl/matlib/rtptools/
addpath /asl/matlib/aslutil
addpath /asl/matlib/h4tools
addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/TIME
addpath /home/sergio/MATLABCODE/PLOTTER
addpath ../GRIB
addpath /home/sergio/MATLABCODE/matlib/clouds/sarta
addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS

system_slurm_stats


%for JOB=1:15
  yy = floor((JOB-1)/12);
  mm = JOB - (yy)*12;
  if mm == 0
    mm = 12;
  end
  yy = yy + 2002;
  fprintf(1,'JOB = %4i --> %4i/%2i \n',JOB,yy,mm)
%end

if JOB < 8
  JOB
  error('not worth processing JOB < 8')
end

dirout = '/asl/s1/sergio/ERA5_monthlyavg_2m/';
fout   = [dirout '/monthly2m_ERA5average_' num2str(yy,'%04d') '_' num2str(mm,'%02d') '_2meter.mat'];
if exist(fout)
  fprintf(1,'%3i %s already exists, quitting \n',JOB,fout);
  error('file already exists')
end

daysINmonth = [31 28 31 30 31 30 31 31 30 31 30 31];
if mod(yy,4) == 0 
  daysINmonth(2) = 29;
end

data_skt = [];
data_t2m = [];
data_d2m = [];
data_time = [];
iCnt = 0;
for dd = 1 : daysINmonth(mm)
  if mod(dd,10) == 0
    fprintf(1,'+');
  else
    fprintf(1,'.');
  end

  fnameIN = ['/asl/models/era5/' num2str(yy,'%04d') '/' num2str(mm,'%02d') '/' num2str(yy,'%04d') num2str(mm,'%02d') num2str(dd,'%02d') '_2meter.nc'];
  sfnameIN = ['/asl/models/era5/' num2str(yy,'%04d') '/' num2str(mm,'%02d') '/' num2str(yy,'%04d') num2str(mm,'%02d') num2str(dd,'%02d') '_sfc.nc'];
  if exist(fnameIN) & exist(sfnameIN)
    iCnt = iCnt + 1;
    s = read_netcdf_lls(sfnameIN);
    x = read_netcdf_lls(fnameIN);
    [mmm,nnn,ooo] = size(x.t2m);
    if ooo ~= 8
      fprintf(1,'JOB %4i/%2i %s size(d2) = %5i %5i %5i error needed 8 steps/day, not %5i \n',JOB,yy,mm,fnameIN,mmm,nnn,ooo,ooo);
      error('wrong array size')
    end
    if iCnt > 1
      data_skt = data_skt + s.skt;
      data_d2m = data_d2m + x.d2m;
      data_t2m = data_t2m + x.t2m;
      data_time = data_time + x.time;
    else
      data_skt = s.skt;
      data_d2m = x.d2m;
      data_t2m = x.t2m;
      data_time = x.time;
    end
  end
end
fprintf(1,'for %4i/%2i expected %2i days and read in data for %2i days \n',yy,mm,daysINmonth(mm),iCnt);

if daysINmonth(mm) ~= iCnt 
  error('iCnt = days files found and daysINmonth(mm) are different')
end

data_skt = data_skt/iCnt;
data_d2m = data_d2m/iCnt;
data_t2m = data_t2m/iCnt;
data_rh2m = airtemp_dewpointtemp_2_RH(data_t2m,data_d2m,2);
data_time = data_time/iCnt;

lat = x.latitude;
lon = x.longitude;

datajunk = squeeze(data_t2m(:,:,1));
[mX,mY] = meshgrid(lat,lon);
figure(1); simplemap(mX(:),mY(:),datajunk); title('t2m');

datajunk = squeeze(data_d2m(:,:,1));
[mX,mY] = meshgrid(lat,lon);
figure(2); simplemap(mX(:),mY(:),datajunk); title('d2m');

datajunk = squeeze(data_rh2m(:,:,1));
[mX,mY] = meshgrid(lat,lon);
figure(3); simplemap(mX(:),mY(:),datajunk); title('RH2m');

for dd = 1 : 3
  figure(dd); colormap jet
end

skt = data_skt;
t2m = data_t2m;
d2m = data_d2m;
rh2m = data_rh2m;
time = data_time;
latitude = lat;
longitude = lon;
comment = 'see /home/sergio/MATLABCODE/RTPMAKE/CLUST_RTPMAKE/CLUSTMAKE_ERA5/clust_make_monthlyavg_2m.m';
saver = ['save ' fout ' skt t2m d2m rh2m longitude latitude time comment'];
eval(saver);
