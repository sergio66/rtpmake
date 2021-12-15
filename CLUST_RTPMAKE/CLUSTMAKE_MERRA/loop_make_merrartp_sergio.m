%% creates an rtp file for ONE granule
%% can be modified for more!

klayers = '/asl/packages/klayers/BinV201/klayers_airs';

klayers = '/home/sergio/klayersV205/BinV201/klayers_airs_x_testmeCLOUDversion';
klayers = '/asl/packages/klayers/Bin/klayers_airs';
klayers = '/asl/packages/klayersV205//BinV201/klayers_airs';

sarta   = '/asl/packages/sartaV108/Bin/sarta_apr08_m140_wcon_nte';
%{
addpath /asl/matlab/science
addpath /asl/matlab/rtptoolsV201
addpath /asl/matlab/h4tools/
addpath /asl/matlab/rtptools/
addpath /asl/matlab/gribtools/
addpath /asl/matlab/airs/readers/
%}

addpath /asl/matlib/opendap
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

iCnt = 0;
lala = +1;
disp('Init : get the list of granules (loop) : prefer to do one at a time!')
while lala > 0
  iCnt = iCnt + 1;
  loop(iCnt).yymmddgg = input('Enter [YY MM DD GG] : ');
  lala = input('another? (-1 to quit) : ');
end

for ix = 1 : iCnt
  clear p h hattr pattr prof yymmddgg
  yymmddgg = loop(ix).yymmddgg;
  ystr = num2str(yymmddgg(1));
  mstr = num2str(yymmddgg(2),'%02d');
  dstr = num2str(yymmddgg(3),'%02d');
  gstr = num2str(yymmddgg(4),'%03d');

  fnamex = ['/strowdataN/s3/robinso2/PUYEHUE/'];
  fnamex = ['/asl/data/rtprod_airs/' ystr '/' mstr '/' dstr '/'];
  fnamex = [fnamex 'allfov_merra' ystr '_' mstr '_' dstr '_' gstr '.rtp'];
  ee = exist(fnamex);
  if ee > 0
    fprintf(1,' %s already exists ... processing next \n',fnamex);
  else
    fprintf(1,' making %s \n',fnamex);
  
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
    filename = ['/strowdataN/data/airs/Aqua_AIRS_Level1/AIRIBRAD.005/' ystr '/'];
    filename = ['/asl/data/airs/AIRIBRAD/' ystr '/'];
    filename = [filename num2str(days_so_far,'%03d') '/'];
    dir0 = filename;
    filename = [filename 'AIRS.' ystr '.' mstr '.' dstr '.' gstr];
    filename = [filename '.L1B.AIRS_Rad.v5*.hdf'];

    thedir = dir(filename);
    if length(thedir) == 1
      fname = [dir0 thedir.name];
    else
      fprintf(1,'%s \n',filename);
      disp('file does not exist');
      iYesNo = input('Get the datafile? (-1/+1) : ');
      if iYesNo > 0
        getstr = ['!getairs ' ystr mstr dstr ' 1 AIRIBRAD.005 ' num2str(gran)];
        eval(getstr);
        thedir = dir(filename);
        fname = [dir0 thedir.name];
        [a,b,c] = sdload(fname);
        scatter(a.Longitude(:),a.Latitude(:),10,a.scanang(:))
        iAHA = input('is this the gran you wanted???? (-1/+1) : ')
        if iAHA < 0
          return
        end
      else
        return
      end
    end

    % [a,b,c] = sdload(fname);
    % p.rlat = [p.rlat a.Latitude(:)'];
    % p.rlon = [p.rlon a.Longitude(:)'];
    % p.rtime = [p.rtime a.Time(:)'];
  
    [meantime, f, prof] = readl1b_all(fname);  
    p = prof;

    p.pobs = zeros(size(p.solazi));
    p.upwell = ones(size(p.solazi));
    %p.irinst = AIRSinst*ones(1,nobs);
    %p.findex = grannum*ones(1,nobs);

    plot(p.rlon,p.rlat,'.')

    pa = {{'profiles','rtime','seconds since 1993'}};
    ha = {{'header','hdf file',filename}};

    h.nchan = 2378;
    h.ichan = (1:2378)';
    h.vchan = f;
    h.pfields=5; % (1=prof + 4=IRobs);

    [h,ha,p,pa] = rtpadd_merra(h,ha,p,pa);     %%% add on ecmwf
    [h,ha,p,pa] = rtpadd_emis_wis(h,ha,p,pa);
  
%{
    [h,ha,p,pa] = rtpadd_emis_DanZhou2(h,ha,p,pa);  
    %% add on emissivity, 99!!! points

    [nemis,efreq,seaemis]=cal_seaemis2(p.satzen,p.wspeed); 
    %% this is better, 19 points
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
      xemis(:,ii) = ...
        interp1(p.efreq(N,ii),p.emis(N,ii),xefreq(:,ii),[],'extrap');
    end

    figure(1)
    %scatter_coast(p.rlon,p.rlat,10,xemis(6,:)); title('emissivity at 954 cm-1')
    scatter_coast(p.rlon,p.rlat,10,p.stemp); title('stemp')
    p.nemis = xnemis;
    p.efreq = xefreq;
    p.emis  = xemis;
    p.nrho  = xnemis;
    p.rfreq = xefreq;
    p.rho   = (1-xemis)/pi;
%}

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    cind1 = [ 532      758       903       1249      1291      2321      2333  ];
    ff    = [0.8224    0.8999    0.9611    1.1295    1.2313    2.6037    2.6164];
    ff = ff*1000;
    t0820 = rad2bt(0822,p.robs1(cind1(1),:));
    t0900 = rad2bt(0899,p.robs1(cind1(2),:));
    t0960 = rad2bt(0961,p.robs1(cind1(3),:));
    t1129 = rad2bt(1129,p.robs1(cind1(4),:));
    t1231 = rad2bt(1231,p.robs1(cind1(5),:));
    t2603 = rad2bt(2603,p.robs1(cind1(6),:));
    t2616 = rad2bt(2616,p.robs1(cind1(7),:));

    figure(2)
      scatter_coast(p.rlon,p.rlat,10,t0960); 
      title('obs BT(960 cm-1)')
    figure(3)
      scatter_coast(p.rlon,p.rlat,10,t0960-t1231); 
      title('obs BT(960 cm-1)-BT(1231 cm-1)')
    caxis([-3 +3]); colorbar

    fip = mktemp('fx.ip.rtp');
    fop = mktemp('fx.op.rtp');
    frp = mktemp('fx.rp.rtp');
    fugh = mktemp('ugh.rtp');

    rtpwrite(fip,h,ha,p,pa);
    klayerser = ['!' klayers ' fin=' fip ' fout=' fop ' >& ' fugh]; 
      eval(klayerser);
    sartaer   = ['!' sarta   ' fin=' fop ' fout=' frp]; 
      eval(sartaer);

    [h2,ha2,p2,pa2] = rtpread(frp);

    %h = h2;
    p.rcalc = p2.rcalc;

    rmer = ['!/bin/rm ' fip ' ' fop ' ' frp ' ' fugh]; eval(rmer);

    %%% now save the rtp file!
    rtpwrite(fnamex,h,ha,p,pa)

    ct0820 = rad2bt(0822,p.rcalc(cind1(1),:));
    ct0900 = rad2bt(0899,p.rcalc(cind1(2),:));
    ct0960 = rad2bt(0961,p.rcalc(cind1(3),:));
    ct1129 = rad2bt(1129,p.rcalc(cind1(4),:));
    ct1231 = rad2bt(1231,p.rcalc(cind1(5),:));
    ct2603 = rad2bt(2603,p.rcalc(cind1(6),:));
    ct2616 = rad2bt(2616,p.rcalc(cind1(7),:));

    figure(4)
      scatter_coast(p.rlon,p.rlat,40,(t2616-t1231)-(ct2616-ct1231)); 
      title('bias : obs-cal : BT(2616 cm-1)-BT(1231 cm-1)')
      caxis([-3 +3]); colorbar

    figure(4)
      scatter_coast(p.rlon,p.rlat,40,(t1129-t1231)); 
      title('bias : obs-cal : BT(1129 cm-1)-BT(1231 cm-1)')
      caxis([-3 +3]); colorbar

    figure(5); %% these are CO channels
    scatter_coast(p.rlon,p.rlat,40,(p.robs1(1865,:)-p.robs1(1867,:)))
    title('CO')

  end
end