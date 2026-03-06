fprintf(1,'days_so_far = %3i \n',days_so_far)

airsdatapathV5 = '/strowdataN/data/airs/Aqua_AIRS_Level1/AIRIBRAD.005/';
airsdatapathV5 = '/asl/data/airs/AIRIBRAD/';

airsdatapathV6 = '/asl/data/airs/L1C/';
airsdatapathV6 = '/asl/data/airs/L1C_v672/';
airsdatapathV6 = '/umbc/xfs3/strow/asl/airs/';

if iv5or6 == 5
  %% L1B
  filename = [airsdatapathV5 ystr '/'];
  filename = [filename num2str(days_so_far,'%03d') '/'];
  dir0 = filename;
  filename = [filename 'AIRS.' ystr '.' mstr '.' dstr '.' gstr];
  filename = [filename '.L1B.AIRS_Rad.v5*.hdf'];
elseif iv5or6 == 6
  %% L1C

  filename = [airsdatapathV6 '/L1C/' ystr '/'];
  filename = [airsdatapathV6 '/L1C_v672/' ystr '/'];
  if yymmddgg(1) <= 2021
    if days_so_far <= 266
      filename = [airsdatapathV6 '/l1c_v672/' ystr '/'];
    else
      filename = [airsdatapathV6 '/l1c_v674/' ystr '/'];
      filename = [airsdatapathV6 '/l1c_v672/' ystr '/'];
    end
  elseif yymmddgg(1) > 2021
    filename = [airsdatapathV6 '/l1c_v674/' ystr '/'];
  end

  %% from Chris H, June 2023
  junkdatestr = [ystr '/' mstr '/' dstr];
  dnum = datenum(junkdatestr,'yyyy/mm/dd');
  if(dnum <= datenum('2021/11/30','yyyy/mm/dd'));
    filename = [airsdatapathV6 '/l1c_v672/'];
  elseif(dnum <= datenum('2023/04/01','yyyy/mm/dd') & dnum > datenum('2021/11/30','yyyy/mm/dd'))
    filename = [airsdatapathV6 '/l1c_v674/'];
  else
    filename = [airsdatapathV6 '/l1c_v675/'];
  end
  filename = [filename '/' ystr '/'];
  filename = [filename num2str(days_so_far,'%03d') '/'];

  dir0 = filename;
  filename = [filename 'AIRS.' ystr '.' mstr '.' dstr '.' gstr];
  filename = [filename '.L1C.AIRS_Rad.v6*.hdf'];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

thedir = dir(filename);
iSimulateData = -1;
if length(thedir) == 1
  fname = [dir0 thedir.name];
  fprintf(1,'reading in %s \n',filename)

else
  fprintf(1,'%s L1B/LiC file does not exist \n',filename);

  %% excess wet bulb, 2020_08_23
  xjunk = '/home/sergio/MATLABCODE/WetBulbTemperatures/usa_2020_08_21.mat';
  xjunk = '/home/sergio/MATLABCODE/WetBulbTemperatures/usa_2020_08_23.mat';  %% hmm, not as excessive as I thought!

  %% excess wet bulb, 2020_08_20
  xjunk = '/home/sergio/MATLABCODE/WetBulbTemperatures/usa_2020_08_20.mat';  %% this one looked better

  %% excess wet bulb, 2019_06_23
  xjunk = '/home/sergio/MATLABCODE/WetBulbTemperatures/usa_2019_06_23.mat';  

  %% MLS by Werner say a big system and I am STRETCHING it
  xjunk = '/home/sergio/MATLABCODE/WetBulbTemperatures/usa_2019_08_27.mat';

  %% https://earthsky.org/earth/study-predicts-deadly-heat-in-persian-gulf/ hot day in the Persian Gulf
  xjunk = '/home/sergio/MATLABCODE/WetBulbTemperatures/middle_east_2015_07_31.mat';
  xjunk = '/home/sergio/MATLABCODE/WetBulbTemperatures/middle_east_2020_07_29.mat';
  xjunk = '/home/sergio/MATLABCODE/WetBulbTemperatures/middle_east_2020_08_23.mat';

  xjunk = '/home/sergio/MATLABCODE/WetBulbTemperatures/middle_east_2015_07_31.mat';
  
  strjunk = ['load ' xjunk '  as needed ? ']; fprintf(1,'xjunk = %s \n',xjunk);
  
  %iSimulateData = input(strjunk);
  %iSimulateData = 1;

  disp(' ')
  error('why should I simulate data? is the AIRS L1B/AIRS L1C .hdf file missing ???????')
  
  if iSimulateData < 0
    return
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if iv5or6 == 5
  [a,b,c] = sdload_quiet(fname);
elseif iv5or6 == 6
  %% find v6_readl2cc.m  /asl/*/rtp_prod2/airs/readers

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
    fprintf(1,'reading in this AIRS L1C data %s \n',fname)
    [eq_x_tai, f, gdata, attr, opt] = read_airicrad(fname);  % Steve
    f0_2645 = f;
  end

  i_vchan_vers = 0;
  if i_vchan_vers == 0

    %% [sergio@c24-52 CLUSTMAKE_ECM]$ ls -lt /home/sergio/git/matlabcode/matlibSergio/matlab2012/srftest/*.hdf
    %% -rw-r--r--+ 1 sergio pi_sergio 10755004 Dec  5  2018 /home/sergio/git/matlabcode/matlibSergio/matlab2012/srftest/airs_l1c_srf_tables_lls_20181205_copy.hdf
    %% -rw-r--r--+ 1 sergio pi_sergio 10755004 Dec  5  2018 /home/sergio/git/matlabcode/matlibSergio/matlab2012/srftest/airs_l1c_srf_tables_lls_20181205.hdf
    %% -rw-rw-r--+ 1 sergio pi_sergio  5385780 Mar 31  2008 /home/sergio/git/matlabcode/matlibSergio/matlab2012/srftest/srftables_m140f_withfake_mar08.hdf

    hdffile = '/home/sergio/MATLABCODE/airs_l1c_srf_tables_lls_20181205.hdf';                                     %% what Chris gave in Dec 2018
    hdffile = '/home/sergio/git/matlabcode/matlibSergio/matlab2012/srftest/airs_l1c_srf_tables_lls_20181205.hdf'; %% what Chris gave in Dec 2018
    vchan2834 = hdfread(hdffile,'freq');
    f = vchan2834;
    
    %% load /home/sergio/git/matlabcode/matlibSergio/matlib/rtp_prod2_Aug11_2020/airs/static/sarta_chans_for_l1c.mat
    load sarta_chans_for_l1c.mat
    theinds2645 = ichan;
    f2645 = f(ichan);
  elseif i_vchan_vers == 1
    a = read_airs_l1c(fname);   %% Chris Hepplewhite
    theinds2645 = cell2mat(a.chanID); theinds2645 = theinds2645';
    f2645   = cell2mat(a.freq);   f2645   = f2645';
  end
  plot(f2645)
  %plot(chanID)
end
