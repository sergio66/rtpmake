%% creates an rtp file for ONE granule
%% can be modified for more!

%[2011 08 01 1414]

%% see /home/schou/OLD/Iasi/mkday_full_new.m

addpath /asl/matlab/iasi/readers/
addpath /asl/matlab/aslutil/
addpath /asl/matlab/science/
addpath /asl/matlab/data_readers/
addpath /asl/matlab/h4toolsV201/
addpath /asl/matlab/iasi/utils
addpath /asl/matlab/iasi/uniform/
addpath /asl/matlab/h4tools/
addpath /asl/matlab/rtptools/
addpath /asl/matlab/gribtools/
addpath /asl/matlab/airs/readers/
addpath /asl/matlab/iasi/clear

indir = '/asl/data/IASI/L1C/';

% Fixed site matchup range {km}
site_range = 55.5;  % we 55.5 km for IASI (and AIRS?)

% Default CO2 ppmv
co2ppm_default = 385.0;

% No data value
nodata = -9999;

iCnt = 0;
yymmddgg = +1;
disp('Init : get the list of granules (loop) : prefer to do one at a time!')
while yymmddgg > 0
  yymmddgg = input('Enter [YY MM DD hhmm] : ');

  ystr = num2str(yymmddgg(1));
  mstr = num2str(yymmddgg(2),'%02d');
  dstr = num2str(yymmddgg(3),'%02d');
  gstr = num2str(yymmddgg(4),'%04d');

  year  = yymmddgg(1);
  month = yymmddgg(2);
  day   = yymmddgg(3);
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
  filename = [indir '/' ystr '/' mstr '/' dstr '/'];
  dir0 = filename;
  filename = [filename 'IASI_xxx_1C_M02_' ystr mstr dstr gstr '*Z'];

  thedir = dir(filename);
  if length(thedir) == 1
    fname = [dir0 thedir.name];
  else
    fprintf(1,'%s \n',filename);
    error('file does not exist');
    end

  [head,hattr,prof,pattr,s,isubset] = iasi_uniform_and_allfov_func(fname,1);

  iCnt = iCnt + 1;
  if iCnt == 1
    p = prof;
  else
    p.findex   = [p.findex prof.findex];
    p.atrack   = [p.atrack prof.atrack];
    p.xtrack   = [p.xtrack prof.xtrack];
    p.zobs     = [p.zobs prof.zobs];
    p.calflag  = [p.calflag prof.calflag];
    p.robs1    = [p.robs1 prof.robs1];
    p.rlat     = [p.rlat prof.rlat];
    p.rlon     = [p.rlon prof.rlon];
    p.rtime    = [p.rtime prof.rtime];
    p.scanang  = [p.scanang prof.scanang];
    p.satzen   = [p.satzen prof.satzen];
    p.satazi   = [p.satazi prof.satazi];
    p.solzen   = [p.solzen prof.solzen];
    p.solazi   = [p.solazi prof.solazi];
    p.salti    = [p.salti prof.salti];
    p.landfrac = [p.landfrac prof.landfrac];
    end
    
  yymmddgg = input('read another file Y/N = (+1/-1) ? ');
  end

disp('Processing ...........')

% Assign RTP attribute strings
hattr = { {'header' 'pltfid' 'IASI'}, ...
          {'header' 'instid' 'METOP2'} };

p.pobs = zeros(size(p.solazi));
p.upwell = ones(size(p.solazi));
%p.irinst = AIRSinst*ones(1,nobs);
%p.findex = grannum*ones(1,nobs);

h = head; 
ha = hattr;
pa = pattr;
clear head hattr pattr

figure(1)
plot(p.rlon,p.rlat,'.')
ix = input('Proceed to adding ECMWF/SARTA? (+1/-1) ');
if ix < 0
  error('lkjfs')
end

addpath /asl/matlab/science
addpath /asl/matlab/rtptoolsV201
h.vchan = (645:0.25:2760)';
h.ichan = 1 : length(h.vchan);
h.nchan = length(h.vchan);
h.pfields=5; % (1=prof + 4=IRobs);

[h,ha,p,pa] = rtpadd_ecmwf_data(h,ha,p,pa);     %%% add on ecmwf
[h,ha,p,pa] = rtpadd_emis_DanZhou2(h,ha,p,pa);  %% add on emissivity, 99! pts

[nemis,efreq,seaemis]=cal_seaemis2(p.satzen,p.wspeed);  %% better, 19 points
%p.nemis = nemis;
%p.efreq = efreq;
%p.emis  = seaemis;
%p.nrho  = nemis;
%p.rfreq = efreq;
%p.rho   = (1-seaemis)/pi;

for ii = 1 : length(p.stemp)
  xnemis(ii)   = nemis(ii);
  xefreq(:,ii) = efreq(:,ii);
  N = 1:p.nemis(ii);
  xemis(:,ii)  = interp1(p.efreq(N,ii),p.emis(N,ii),xefreq(:,ii),[],'extrap');
  end

addpath /home/sergio/MATLABCODE/PLOTTER/
figure(1)
scatter_coast(p.rlon,p.rlat,10,xemis(6,:)); title('emissivity at 954 cm-1')
p.nemis = xnemis;
p.efreq = xefreq;
p.emis  = xemis;
p.nrho  = xnemis;
p.rfreq = xefreq;
p.rho   = (1-xemis)/pi;

[h1,ha1,p1,h2,ha2,p2] = iasi_split_prof(h,p,ha);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ff    = [0.8224    0.8999    0.9611    1.1295    1.2313    2.6037    2.6164]*1000;
for ii = 1 : length(ff)
  junk = find(h.vchan >= ff(ii)); 
  cind1(ii) = junk(1);
end
  
t0820 = rad2bt(0822,p1.robs1(cind1(1),:));
t0960 = rad2bt(0961,p1.robs1(cind1(3),:));
t1231 = rad2bt(1231,p1.robs1(cind1(5),:));
addpath /home/sergio/MATLABCODE/PLOTTER/

figure(2)
  scatter_coast(p.rlon,p.rlat,30,t0960); title('obs BT(960 cm-1)')
  caxis([200 300]); colorbar
figure(3)
  scatter_coast(p.rlon,p.rlat,30,t0960-t1231); title('obs BT(960 cm-1)-BT(1231 cm-1)')
  caxis([-3 +3]); colorbar

klayers = '/asl/packages/klayers/Bin/klayers_airs';
sarta   = '/asl/packages/sartaV108/Bin/sarta_apr08_m140_wcon_nte';
sarta   = '/asl/packages/sartaV108/Bin/sarta_iasi_may09_iceaggr_waterdrop_biom_slabcloud_hg3_wcon_nte_swch4';
sarta = '/asl/packages/sartaV108/Bin/sarta_iasi_may09_wcon_nte_swch4_nh3';

%%%%%%%%%%%%%%%%%%%%%%%%%
fip  = mktemp('fx.ip.rtp');
fop  = mktemp('fx.op.rtp');
frp  = mktemp('fx.rp.rtp');
fugh = mktemp('ugh.rtp');

oldrtpwrite(fip,h1,ha1,p1,pa);
klayerser = ['!' klayers ' fin=' fip ' fout=' fop ' >& ' fugh]; eval(klayerser);
sartaer   = ['!' sarta   ' fin=' fop ' fout=' frp]; eval(sartaer);
[hx,hax,px,pax] = oldrtpread(frp);
p1.rcalc = px.rcalc;
ct0820 = rad2bt(0822,p1.rcalc(cind1(1),:));
ct0960 = rad2bt(0961,p1.rcalc(cind1(3),:));
ct1231 = rad2bt(1231,p1.rcalc(cind1(5),:));
figure(4)
  scatter_coast(p1.rlon,p1.rlat,10,(t0960-t1231)-(ct0960-ct1231)); 
  title('bias : obs-cal : BT(960 cm-1)-BT(1231 cm-1)')
  caxis([-3 +3]); colorbar

rmer = ['!/bin/rm ' fip ' ' fop ' ' frp ' ' fugh]; eval(rmer);

%%% now save the rtp file!
outdir = ['/asl/s1/sergio/IASI_DUST/' ystr '/' mstr '/' dstr '/'];
ee = exist(outdir);
if ee == 0
  mker = ['!/bin/mkdir -p ' outdir];
  eval(mker)
end

fnamex = [outdir 'iasi_allfov_' ystr '_' mstr '_' dstr '_' gstr '_pt1.rtp'];
oldrtpwrite(fnamex,h1,ha1,p1,pa)

%%%%%%%%%%%%

oldrtpwrite(fip,h2,ha2,p2,pa);
klayerser = ['!' klayers ' fin=' fip ' fout=' fop ' >& ' fugh]; eval(klayerser);
sartaer   = ['!' sarta   ' fin=' fop ' fout=' frp]; eval(sartaer);
[hx,hax,px,pax] = oldrtpread(frp);
p2.rcalc = px.rcalc;

rmer = ['!/bin/rm ' fip ' ' fop ' ' frp ' ' fugh]; eval(rmer);

%%% now save the rtp file!
fnamex = [outdir 'iasi_allfov_' ystr '_' mstr '_' dstr '_' gstr '_pt2.rtp'];
oldrtpwrite(fnamex,h2,ha2,p2,pa)