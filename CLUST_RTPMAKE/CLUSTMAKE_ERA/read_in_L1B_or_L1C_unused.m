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
  filename = ['/asl/airs/l1c_v672/' ystr '/'];
  filename = [filename num2str(days_so_far,'%03d') '/'];
  dir0 = filename;
  filename = [filename 'AIRS.' ystr '.' mstr '.' dstr '.' gstr];
  filename = [filename '.L1C.AIRS_Rad.v6*.hdf'];
end

thedir = dir(filename);
if length(thedir) == 1
  fname = [dir0 thedir.name];
else
  fprintf(1,'%s \n',filename);
  disp('file does not exist');
  return
end

if iv5or6 == 5
  [a,b,c] = sdload_quiet(fname);
elseif iv5or6 == 6
  %% find v6_readl2cc.m  /asl/*/rtp_prod2/airs/readers
  %addpath /asl/packages/rtp_prod2/airs/readers/
  addpath /home/sergio/MATLABCODE/matlib/rtp_prod2/airs/readers
  addpath /asl/matlib/time
  [eq_x_tai, f, gdata, attr, opt] = read_airicrad(fname);  % Steve
  f0_2645 = f;

  hdffile = '/home/sergio/MATLABCODE/airs_l1c_srf_tables_lls_20181205.hdf';   % what he gave in Dec 2018
  vchan2834 = hdfread(hdffile,'freq');
  f = vchan2834;
  load sarta_chans_for_l1c.mat
  theinds2645 = ichan;
  f2645 = f(ichan);
  
  %      a = read_airs_l1c(fname);   %% Chris Hepplewhite
  %      theinds2645 = cell2mat(a.chanID); theinds2645 = theinds2645';
  %      f2645   = cell2mat(a.freq);   f2645   = f2645';
  % plot(f2645)
  % plot(chanID)
end
