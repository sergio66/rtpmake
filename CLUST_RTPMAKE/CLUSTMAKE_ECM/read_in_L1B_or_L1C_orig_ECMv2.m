if iv5or6 == 5
  %% L1B
  % filename = ['/strowdataN/data/airs/Aqua_AIRS_Level1/AIRIBRAD.005/' ystr '/'];
  filename = ['/asl/data/airs/AIRIBRAD/' ystr '/'];
  filename = [filename num2str(days_so_far,'%03d') '/'];
  dir0 = filename;
  filename = [filename 'AIRS.' ystr '.' mstr '.' dstr '.' gstr];
  filename = [filename '.L1B.AIRS_Rad.v5*.hdf'];
elseif iv5or6 == 6
  %% L1C
  filename = ['/asl/data/airs/L1C/' ystr '/'];
  filename = ['/asl/data/airs/L1C_v672/' ystr '/'];
  if yymmddgg(1) <= 2021
    if days_so_far <= 266
      filename = ['/asl/airs/l1c_v672/' ystr '/'];
    else
      filename = ['/asl/airs/l1c_v674/' ystr '/'];
    end
  elseif yymmddgg(1) > 2021
    filename = ['/asl/airs/l1c_v674/' ystr '/'];
  end
  filename = [filename num2str(days_so_far,'%03d') '/'];

  dir0 = filename;
  filename = [filename 'AIRS.' ystr '.' mstr '.' dstr '.' gstr];
  filename = [filename '.L1C.AIRS_Rad.v6*.hdf'];
end

thedir = dir(filename);
iSimulateData = -1;
if length(thedir) == 1
  fname = [dir0 thedir.name];
else
  fprintf(1,'%s \n',filename);
  disp('file does not exist');

  %% excess wet bulb, 2019_06_23
  xjunk = '/home/sergio/MATLABCODE/WetBulbTemperatures/usa_2019_06_23.mat';  

  %% MLS by Werner say a big system and I am STRETCHING it
  xjunk = '/home/sergio/MATLABCODE/WetBulbTemperatures/usa_2019_08_27.mat';

  %% https://earthsky.org/earth/study-predicts-deadly-heat-in-persian-gulf/ hot day in the Persian Gulf
  xjunk = '/home/sergio/MATLABCODE/WetBulbTemperatures/middle_east_2015_07_31.mat';
  xjunk = '/home/sergio/MATLABCODE/WetBulbTemperatures/middle_east_2020_08_23.mat';
  xjunk = '/home/sergio/MATLABCODE/WetBulbTemperatures/middle_east_2020_07_29.mat';
  xjunk = '/home/sergio/MATLABCODE/WetBulbTemperatures/middle_east_2020_07_29pm2.mat';

  %% excess wet bulb, 2020_08_23
  xjunk = '/home/sergio/MATLABCODE/WetBulbTemperatures/usa_2020_08_21.mat';
  xjunk = '/home/sergio/MATLABCODE/WetBulbTemperatures/usa_2020_08_23.mat';    %% hmm, not as excessive as I thought!
  xjunk = '/home/sergio/MATLABCODE/WetBulbTemperatures/usa_2020_08_20.mat';    %% this one looked better
  xjunk = '/home/sergio/MATLABCODE/WetBulbTemperatures/usa_2020_08_23pm2.mat'; %% hmm, not as excessive as I thought!

  xjunk = '/home/sergio/MATLABCODE/WetBulbTemperatures/india_2019_05_11.mat';
  xjunk = '/home/sergio/MATLABCODE/WetBulbTemperatures/india_2019_05_11pm2.mat';


  %%%%%%%%%%%%%%%%%%%%%%%%%

  strjunk = ['load ' xjunk '  as needed ? ']; fprintf(1,'xjunk = %s \n',xjunk);
  %iSimulateData = input(strjunk);
  iSimulateData = 1

  if iSimulateData < 0
    return
  end
end

if iv5or6 == 5
  [a,b,c] = sdload_quiet(fname);
elseif iv5or6 == 6
  %% find v6_readl2cc.m  /asl/*/rtp_prod2/airs/readers
  addpath /home/sergio/MATLABCODE/matlib/rtp_prod2/airs/readers
  addpath /asl/matlib/time

  if iSimulateData > 0
    fprintf(1,'loading %s \n',xjunk);
    pjunk = load(xjunk); 

    gdata.rtime = pjunk.p0.rtime;
    gdata.rlat = pjunk.p0.rlat;
    gdata.rlon = pjunk.p0.rlon;
    gdata.solazi = ones(size(pjunk.p0.rlon)) * 00;
    gdata.solzen = ones(size(pjunk.p0.rlon)) * 150;
      wonk = length(gdata.solzen);
      gdata.solzen(1:wonk/2) = 40; 
    gdata.satzen = ones(size(pjunk.p0.rlon)) * 22;
    gdata.satazi = ones(size(pjunk.p0.rlon)) * 0;
    gdata.robs1 = zeros(2645,length(gdata.rtime));
   
    [gdata.salti, gdata.landfrac] = usgs_deg10_dem(gdata.rlat, gdata.rlon);
  else
    [eq_x_tai, f, gdata, attr, opt] = read_airicrad(fname);  % Steve
    f0_2645 = f;

    hdffile = '/home/sergio/MATLABCODE/airs_l1c_srf_tables_lls_20181205.hdf';   % what he gave in Dec 2018
    vchan2834 = hdfread(hdffile,'freq');
    f = vchan2834;
    load sarta_chans_for_l1c.mat
    theinds2645 = ichan;
    f2645 = f(ichan);
  
    f = f2645;
    theinds = theinds2645;
    
    %      a = read_airs_l1c(fname);   %% Chris Hepplewhite
    %      theinds2645 = cell2mat(a.chanID); theinds2645 = theinds2645';
    %      f2645   = cell2mat(a.freq);   f2645   = f2645';
    % plot(f2645)
    % plot(chanID)

  end  
end
