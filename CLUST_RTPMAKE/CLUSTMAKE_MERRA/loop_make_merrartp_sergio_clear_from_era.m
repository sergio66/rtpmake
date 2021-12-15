% local running to test
% clustcmd -L loop_make_merrartp_sergio_clear_from_era.m 20020801:16:now
%
% otherwise when happy
% clustcmd -q medium -n 16 -p 2 loop_make_merrartp_sergio_clear_from_era.m 20020801:1:now
%
% or
% clustcmd -q long_contrib -n 16 -p 4 loop_make_merrartp_sergio_clear_from_era.m 20020801:now

%% creates an rtp file for ONE granule
%% can be modified for more!

klayers = '/asl/packages/klayers/BinV201/klayers_airs';

klayers = '/home/sergio/klayersV205/BinV201/klayers_airs_x_testmeCLOUDversion';
klayers = '/asl/packages/klayers/Bin/klayers_airs';
klayers = '/asl/packages/klayersV205//BinV201/klayers_airs';

sarta   = '/asl/packages/sartaV108/Bin/sarta_apr08_m140_wcon_nte';

addpath /asl/matlib/opendap
addpath /asl/matlib/science
%addpath /asl/matlib/rtptoolsV201
addpath /asl/matlib/h4tools/
addpath /asl/matlib/rtptools/
addpath /asl/matlib/gribtools/
%addpath /asl/matlib/airs/readers/

%addpath /home/shared/sergio/MATLABCODE
%addpath /asl/matlib/cris/readers/
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

yymmddgg = input('Enter [YY MM DD] : ');

%% [YY MM DD h m s] = datevec(JOB(1));
%% yymmddgg = [YY MM DD];

for ix = 0 : 23
  clear p h hattr pattr prof  hclr pclr
  ystr = num2str(yymmddgg(1));
  mstr = num2str(yymmddgg(2),'%02d');
  dstr = num2str(yymmddgg(3),'%02d');
  gstr = num2str(ix,'%02d');

  fnameDIR = ['/asl/data/rtprod_airs_0/' ystr '/' mstr '/' dstr '/'];
  fnamex   = [fnameDIR 'sergio_clr_merra.' ystr '.' mstr '.' dstr '.' gstr '.rtp'];
  fnameERA = [fnameDIR 'era.airs_l1bcm.' ystr '.' mstr '.' dstr '.' gstr '.v1.rtpZ'];

  ee = exist(fnamex);
  eeR = exist(fnameERA);

  if  ee > 0
    thedir = dir(fnamex);
    if thedir.bytes < 1000
      ee = 0;
    end
  end

  if  eeR > 0
    thedir = dir(fnameERA);
    if thedir.bytes < 1000
      eeR = 0;
    end
  end
  
  if ee > 0 & eeR > 0
    fprintf(1,' %s already exists ... processing next \n',fnamex);
  elseif ee > 0 & eeR == 0
    fprintf(1,' %s already exists ... but ERA does not ... processing next \n',fnamex);
  elseif ee == 0 & eeR == 0
    fprintf(1,' %s does not exist  ... and neither does ERA ... processing next \n',fnamex);
  else
    fprintf(1,' making %s \n',fnamex);

    year  = yymmddgg(1);
    month = yymmddgg(2);
    day   = yymmddgg(3);
    gran  = ix;

    %%% read in clr, effectively doing     [meantime, f, prof] = readl1b_all(fname);  

    [hclr,ha,pclr,pa] = rtpread(fnameERA);
    [hclr,ha,pclr,pa] = rtpgrow(hclr,ha,pclr,pa);
    iclear = find(pclr.iudef(1,:) == 1);

    prof.findex = pclr.findex(iclear);
    prof.atrack = pclr.atrack(iclear);
    prof.xtrack = pclr.xtrack(iclear);
    prof.zobs = pclr.zobs(iclear) * 1000;
    prof.calflag = pclr.calflag(:,iclear);
    prof.robs1 = pclr.robs1(:,iclear);
    prof.rlat = pclr.rlat(iclear);
    prof.rlon = pclr.rlon(iclear);
    prof.rtime = pclr.rtime(iclear);
    prof.scanang = pclr.scanang(iclear);
    prof.satzen = pclr.satzen(iclear);
%    prof.satazi = pclr.satazi(iclear);
    prof.solzen = pclr.solzen(iclear);
%    prof.solazi = pclr.solazi(iclear);
    prof.salti = pclr.salti(iclear);
    prof.landfrac = pclr.landfrac(iclear);

    p = prof;

    p.pobs = zeros(size(p.xtrack));
    p.upwell = ones(size(p.xtrack));
    %p.irinst = AIRSinst*ones(1,nobs);
    %p.findex = grannum*ones(1,nobs);

    plot(p.rlon,p.rlat,'.')

    pa = {{'profiles','rtime','seconds since 1993'}};
    ha = {{'header','era file',fnameERA}};

    h.nchan = 2378;
    h.ichan = (1:2378)';
    h.vchan = instr_chans('airs');
    h.vchan = hclr.vchan;;
    h.pfields=5; % (1=prof + 4=IRobs);

    [h,ha,p,pa] = rtpadd_merra(h,ha,p,pa);     %%% add on merra
    [h,ha,p,pa] = rtpadd_emis_wis(h,ha,p,pa);
  
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
    sartaer   = ['!' sarta   ' fin=' fop ' fout=' frp ' >& ' fugh]; 
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