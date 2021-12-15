%% creates an rtp file for ONE granule
%% can be modified for more!

klayers = '/asl/packages/klayersV205//BinV201/klayers_airs';

run_sarta.clear = -1;
run_sarta.cloud = +1;
run_sarta.cumsum = 9999;

iCodeX = 0;    %%% use Baran ice, and our water (done by Scott)
iCodeX = 1;    %%% use B.Baum ice, and our water (done by Sergio)
iCodeX = 2;    %%% use B.Baum ice, and MODIS variance water (done by Sergio)

iCodeX = 1;    %%% use B.Baum ice, and our water (done by Sergio); pretty much ran off 2012/09/20 AIRS with this, so do same for CrIS

if iCodeX == 0
  run_sarta.sartacloud_code = ...
    '/asl/packages/sartaV108/Bin/sarta_crisg4_nov09_iceaggr_waterdrop_desertdust_slabcloud_hg3_wcon_nte';
elseif iCodeX == 1
  run_sarta.sartacloud_code = ...
    '/asl/packages/sartaV108/Bin/sarta_crisg4_nov09_iceGHMbaum_waterdrop_desertdust_slabcloud_hg3_wcon_nte';
elseif iCodeX == 2
  run_sarta.sartacloud_code = ...
    '/asl/packages/sartaV108/Bin/sarta_crisg4_nov09_iceGHMbaum_waterdropMODIS_desertdust_slabcloud_hg3_wcon_nte';
end

addpath /home/sergio/MATLABCODE/matlib/clouds/sarta
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

yymmddgg = input('Enter [YY MM DD] : ');

%for ix = 0 : 143
for ix = 143 : -1 : 0
  clear p h hattr pattr prof  hclr pclr
  ystr = num2str(yymmddgg(1));
  mstr = num2str(yymmddgg(2),'%02d');
  dstr = num2str(yymmddgg(3),'%02d');
  gstr = num2str(ix,'%03d');

  fnameDIR = ['/asl/data/rtprod_cris/' ystr '/' mstr '/' dstr '/'];
  fnameECM = [fnameDIR 'ecm.cris_cspp_dev.' ystr '.' mstr '.' dstr '.' gstr '.Rv1.1d-Mv1.1c-1-g90d9ac4.rtpZ'];

  if iCodeX == 0
    fnamex   = [fnameDIR 'sergio_cris_cld_ecm.' ystr '.' mstr '.' dstr '.' gstr '.rtp'];
  elseif iCodeX == 1
    fnamex   = [fnameDIR 'sergio_cris_cld_ecm_iceGHMbaum.' ystr '.' mstr '.' dstr '.' gstr '.rtp'];
  elseif iCodeX == 2
    fnamex   = [fnameDIR 'sergio_cris_cld_ecm_iceGHMbaum_waterMODIS.' ystr '.' mstr '.' dstr '.' gstr '.rtp'];
  end

  ee = exist(fnamex);
  eeR = exist(fnameECM);

  if  ee > 0
    thedir = dir(fnamex);
    if thedir.bytes < 1000
      ee = 0;
    end
  end

  if  eeR > 0
    thedir = dir(fnameECM);
    if thedir.bytes < 1000
      eeR = 0;
    end
  end
  
  if ee > 0 & eeR > 0
    fprintf(1,' %s already exists ... processing next \n',fnamex);
  elseif ee > 0 & eeR == 0
    fprintf(1,' %s already exists ... but ECM does not ... processing next \n',fnamex);
  elseif ee == 0 & eeR == 0
    fprintf(1,' %s does not exist  ... and neither does ECM ... processing next \n',fnamex);
  else
    fprintf(1,' making %s \n',fnamex);

    year  = yymmddgg(1);
    month = yymmddgg(2);
    day   = yymmddgg(3);
    gran  = ix;

    %%% read in clr, effectively doing     [meantime, f, prof] = readl1b_all(fname);  

    [hclr,ha,prof,pa] = rtpread(fnameECM);
    [hclr,ha,prof,pa] = rtpgrow(hclr,ha,prof,pa);

    p = prof;
    h = hclr;

    p.pobs = zeros(size(p.xtrack));
    p.upwell = ones(size(p.xtrack));
    %p.irinst = AIRSinst*ones(1,nobs);
    %p.findex = grannum*ones(1,nobs);

    plot(p.rlon,p.rlat,'.')

    cldfields = {'CC','CIWC','CLWC'};
    [h,ha,p,pa] = rtpadd_ecmwf_data(h,ha,p,pa,cldfields); %%% add on ecm

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [p2] = driver_sarta_cloud_rtp(h,ha,p,pa,run_sarta);
    rtpwrite(fnamex,h,ha,p2,pa)

    % [h,ha,p2x,pa] = rtptrim_sartacloud(h,ha,p2,pa);
    % rtpwrite(fnamex,h,ha,p2x,pa)

  end
end