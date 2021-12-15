%% creates an rtp file for ONE granule
%% can be modified for more!

klayers = '/asl/packages/klayersV205/BinV201/klayers_airs';
sarta   = '/asl/packages/sartaV108/Bin/sarta_crisg4_nov09_wcon_nte_nh3';

addpath /home/sergio/Backup_asl_matlab_Feb2013/aslutil
addpath /home/sergio/Backup_asl_matlab_Feb2013//iasi/utils

addpath /asl/matlib/science
addpath /asl/matlib/rtptoolsV201
addpath /asl/matlib/h4tools/
addpath /asl/matlib/rtptools/
addpath /asl/matlib/gribtools/
addpath /asl/matlib/airs/readers/

addpath /home/shared/sergio/MATLABCODE
addpath /asl/matlib/cris/readers/
addpath /asl/matlib/science
addpath /asl/matlib/science/
addpath /home/sergio/MATLABCODE/PLOTTER
addpath /home/sergio/MATLABCODE/DUSTFLAG_CRiS/

addpath /home/sergio/MATLABCODE/DUSTFLAG_CRiS/

addpath /asl/rtp_prod/cris/readers
addpath /asl/rtp_prod/cris/utils
addpath /home/sergio/MATLABCODE/
addpath /asl/matlib/aslutil

freqsNchansCRIS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

warning off    %% else it cries about /asl/matlab/ instead of /asl/matlib or matlib

yymmddgg = input('Enter [YY MM DD HHMM] : ');
ystr = num2str(yymmddgg(1));
mstr = num2str(yymmddgg(2),'%02d');
dstr = num2str(yymmddgg(3),'%02d');
gstr = num2str(yymmddgg(4),'%03d');

year  = yymmddgg(1);
month = yymmddgg(2);
day   = yymmddgg(3); 
gran  = yymmddgg(4);
if mod(year,4) == 0
  mos = [31 29 31 30 31 30 31 31 30 31 30 31];  %% leap year
else
  mos = [31 28 31 30 31 30 31 31 30 31 30 31];  %% normal year
  end
days_so_far = 0;
if month > 1
  days_so_far = sum(mos(1:month-1));
  end
days_so_far = days_so_far + day;
filename = ['/asl/data/cris/sdr60/hdf/' ystr '/'];
filename = [filename num2str(days_so_far,'%03d') '/'];
dir0 = filename;
hstr = num2str(yymmddgg(4),'%02d');
filename = [filename 'SCRIS_npp_d' ystr mstr dstr '_t' hstr '*.h5']

thedir = dir(filename);
if length(thedir) == 1
  fname = [dir0 thedir.name];
elseif length(thedir) > 1
  for ii = 1 : length(thedir)
    xname =  thedir(ii).name;
    fprintf(1,'%3i %s \n',ii,xname);
  end
  ii = input('enter your choice : ');
  fname = [dir0 thedir(ii).name];  
else
  fprintf(1,'%s \n',filename);
  error('file does not exist');
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[p,pa] = readsdr_rtp(fname);

if ~isfield(p,'landfrac')
  [salti, landfrac] = usgs_deg10_dem(p.rlat, p.rlon);
  p.landfrac = landfrac;
end

  lala = find(p.rlat < -90); p.rlat(lala) = -89.9;
  lala = find(p.rlat > +90); p.rlat(lala) = +89.9;
  lala = find(p.rlon < -180); p.rlon(lala) = -179.9;
  lala = find(p.rlon > +180); p.rlon(lala) = +179.9;

iApodize = +1;  %% obs are not apodized, need to do that!
if iApodize > 0
  disp('apodizing the obs ....');
  addpath /asl/matlab2012/cris/unapod
  xobs = p.robs1;
  xham = box_to_ham(xobs);
  p.robs1 = xham;
end

h.vchan = instr_chans('cris');
h.ichan = (1:length(h.vchan))';

freqsNchansCRIS
figure(1); t0961 = rad2bt(ff(3),p.robs1(cind1(3),:));
           t1231 = rad2bt(ff(5),p.robs1(cind1(5),:));
figure(1); scatter_coast(p.rlon,p.rlat,20,t0961-t1231); title('bt961 - bt1231');

figure(2); scatter_coast(p.rlon,p.rlat,20,p.landfrac); title('landfrac');
iYes = input('Is this the file you want to process (-1/+1) ');
if iYes < 0
  disp('oops try again!');
  return
else
  disp('ok, proceeding with klayers and sarta')
end

fnamePaul = ['/asl/data/rtprod_cris/' ystr '/' mstr '/' dstr '/'];
fnamePaul= [fnamePaul 'sergio_cris_l1c_ecm.' ystr '.' mstr '.' dstr '.' gstr '.rtp'];
eeP = exist(fnamePaul);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if eeP == 0
  fprintf(1,' making %s \n',fnamePaul);

  p.pobs = zeros(size(p.solazi));
  p.upwell = ones(size(p.solazi));
  %p.irinst = AIRSinst*ones(1,nobs);
  %p.findex = grannum*ones(1,nobs);

  plot(p.rlon,p.rlat,'.')

  h.pfields=5; % (1=prof + 4=IRobs);
  h.nchan = length(h.ichan);
    
  %pa = {{'profiles','rtime','seconds since 1993'}};
  ha = {{'header','hdf file',filename}};

  [h,ha,p,pa] = rtpadd_ecmwf_data(h,ha,p,pa);     %%% add on ecm
  [h,ha,p,pa] = rtpadd_emis_DanZhou2(h,ha,p,pa);  
  %% add on emissivity, 99!!! points

  [nemis,efreq,seaemis]=cal_seaemis2(p.satzen,p.wspeed); 
  %% this is better, 19 points
  p.nemis = nemis;
  p.efreq = efreq;
  p.emis  = seaemis;
  p.nrho  = nemis;
  p.rfreq = efreq;
  p.rho   = (1-seaemis)/pi;

  for ii = 1 : length(p.stemp)
    xnemis(ii)   = nemis(ii);
    xefreq(:,ii) = efreq(:,ii);
    N = 1:p.nemis(ii);
    xemis(:,ii) = ...
      interp1(p.efreq(N,ii),p.emis(N,ii),xefreq(:,ii),[],'extrap');
  end

  figure(1)
  scatter_coast(p.rlon,p.rlat,10,xemis(6,:)); title('emissivity at 954 cm-1')

  p.nemis = xnemis;
  p.efreq = xefreq;
  p.emis  = xemis;
  p.nrho  = xnemis;
  p.rfreq = xefreq;
  p.rho   = (1-xemis)/pi;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  fip = mktemp('fx.ip.rtp');
  fop = mktemp('fx.op.rtp');
  frp = mktemp('fx.rp.rtp');
  fugh = mktemp('ugh.rtp');

  rtpwrite(fip,h,ha,p,pa);
  disp(' running klayers ....')

  klayerser = ['!' klayers ' fin=' fip ' fout=' fop ' >& ' fugh]; 
    eval(klayerser);
  disp(' running sarta ....')
  sartaer   = ['!' sarta   ' fin=' fop ' fout=' frp ' >& ' fugh]; 
    eval(sartaer);

  [h2,ha2,p2,pa2] = rtpread(frp);

  %h = h2;
  p.rcalc = p2.rcalc;

  rmer = ['!/bin/rm ' fip ' ' fop ' ' frp ' ' fugh]; eval(rmer);

  %%% now save the rtp file!

  fnamex = fnamePaul;
  rtpwrite(fnamex,h,ha,p,pa)
  fprintf(1,'saved off %s \n',fnamePaul);
else 
  fprintf(1,'well %s already exists, so did not do anyything! \n',fnamePaul);
end
