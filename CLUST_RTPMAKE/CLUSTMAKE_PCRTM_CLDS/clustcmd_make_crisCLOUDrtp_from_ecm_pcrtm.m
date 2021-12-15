% local running to test
% clustcmd -L clustcmd_make_crisCLOUDrtp_from_ecm_pcrtm.m 001:144
%
% otherwise when happy
% clustcmd -q medium -n 64 -p 2 clustcmd_make_crisCLOUDrtp_from_ecm_pcrtm.m 001:144
%
% or
% clustcmd -q long_contrib -n 64 -p 4 clustcmd_make_crisCLOUDrtp_from_ecm_pcrtm.m 001:144

%% copied from clustcmd_make_crisCLOUDrtp_from_ecm.m

klayers = '/asl/packages/klayersV205//BinV201/klayers_airs';

%%% set sarta exec

run_sarta.clear = -1;
run_sarta.cloud = -1;
run_sarta.cumsum = +9999;

codeX = 0; %% use default
codeX = 1; %% use new B Baum

code0 = '/asl/packages/sartaV108/BinV201/sarta_apr08_m140_iceaggr_waterdrop_desertdust_slabcloud_hg3_wcon_nte';
code1 = '/home/sergio/SARTA_CLOUDY/BinV201/sarta_apr08_m140x_iceGHMbaum_waterdrop_desertdust_slabcloud_hg3';

if codeX == 0 & run_sarta.cloud > 0
  icestr = '_pcrtm_and_sarta_baran_ice';
  run_sarta.sartacloud_code = code0;
elseif codeX == 1 & run_sarta.cloud > 0
  icestr = '_pcrtm_and_sarta_baum_ice';
  run_sarta.sartacloud_code = code1;
elseif run_sarta.cloud < 0
  icestr = '_pcrtm_only';
else
  error('codeX???')
end

% icestr = ['cloudy_airs_l1b_ecm' icestr '.'];
% icestr = ['cloudy_airs_l1b_ecm_fixed' icestr '.'];

addpath /home/sergio/MATLABCODE/matlib/clouds/sarta
addpath /home/sergio/MATLABCODE/matlib/clouds/pcrtm
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

%yymmddgg = input('Enter [YY MM DD] : ');
yymmddgg = [2012 09 20];

for ix = JOB
  clear p h hattr pattr prof  hclr pclr
  ystr = num2str(yymmddgg(1));
  mstr = num2str(yymmddgg(2),'%02d');
  dstr = num2str(yymmddgg(3),'%02d');
  gstr = num2str(ix,'%03d');

  fnameDIR = ['/asl/data/rtprod_cris/' ystr '/' mstr '/' dstr '/'];
  fnameECM = [fnameDIR 'ecm.cris_cspp_dev.' ystr '.' mstr '.' dstr '.' gstr '.Rv1.1d-Mv1.1c-1-g90d9ac4.rtpZ'];

  fnamex   = [fnameDIR 'sergio_pcrtm_cris_cld_ecm.' ystr '.' mstr '.' dstr '.' gstr '.rtp'];

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

    toucher = ['!/bin/touch ' fnamex]; eval(toucher);


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
    if ~isfield(p,'scanang')
      scanang = saconv(p.satzen, 805000*ones(size(p.satzen)));    
      p.scanang = scanang;
    end
 
    % I asked Breno to include cloud fields a few days ago
    % cldfields = {'CC','CIWC','CLWC'};
    % [h,ha,p,pa] = rtpadd_ecmwf_data(h,ha,p,pa,cldfields); %%% add on ecm

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    p.cfrac_ERA_ECM_seed = p.cfrac;
    [p2] = driver_pcrtm_cloud_rtp(h,ha,p,pa,run_sarta);
    rtpwrite(fnamex,h,ha,p2,pa)

    % [h,ha,p2x,pa] = rtptrim_sartacloud(h,ha,p2,pa);
    % rtpwrite(fnamex,h,ha,p2x,pa)

  end
end