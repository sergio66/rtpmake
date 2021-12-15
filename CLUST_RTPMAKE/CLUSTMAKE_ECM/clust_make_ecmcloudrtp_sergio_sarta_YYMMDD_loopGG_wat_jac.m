%% have to hardcode  YYMMDD before starting, loops over granules

% local running to test
% clustcmd -L clust_make_ecmcloudrtp_sergio_sarta_YYMMDD_loopGG_wat_jac.m 001:240
%
% otherwise when happy
% clustcmd -q medium -n 64 -p 4 clust_make_ecmcloudrtp_sergio_sarta_YYMMDD_loopGG_wat_jac.m 001:240
%
% or
% clustcmd -q long_contrib -n 64 -p 4 clust_make_ecmcloudrtp_sergio_sarta_YYMMDD_loopGG_wat_jac.m 001:240

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% March 11, 2011 is my good example
%yymmdd0  = input('Enter [YYYY MM DD] : ');
%iaGlist = input('Enter [GranStart GranList] : ');

yymmdd0  = [2012 09 20];   %% for SNO with CrIS
yymmdd0  = [2011 03 11];   %% for JGR paper
iaGlist  = JOB;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% creates an rtp file for ONE granule
%% can be modified for more!

klayers = '/asl/packages/klayers/Bin/klayers_airs';
sarta   = '/asl/packages/sartaV108/Bin/sarta_apr08_m140_wcon_nte';

klayers = '/asl/packages/klayers/Bin/klayers_airs';
sarta   = '/asl/packages/sartaV108/Bin/sarta_apr08_m140_wcon_nte';

klayers = '/asl/packages/klayersV205/BinV201/klayers_airs';
sarta   = '/asl/packages/sartaV108_PGEv6/Bin/sarta_airs_PGEv6_postNov2003';

addpath /asl/matlib/aslutil
addpath /asl/matlib/science
addpath /asl/matlib/rtptools
addpath /asl/matlib/h4tools/
addpath /asl/matlib/rtptools/
addpath /asl/matlib/gribtools/
addpath /home/sergio/MATLABCODE/matlib/clouds/sarta

theinds = (1 : 2378)';

%%%%%%%%%%%%%%%%%%%%%%%%%
%%% set sarta exec

run_sarta.clear = +1;
run_sarta.cloud = +1;
run_sarta.cumsum = 9999;

codeX = 0; %% use default with A. Baran params
codeX = 1; %% use new     with B. Baum, P. Yang params

code0 = '/asl/packages/sartaV108/BinV201/sarta_apr08_m140_iceaggr_waterdrop_desertdust_slabcloud_hg3_wcon_nte';
code1 = '/home/sergio/SARTA_CLOUDY/BinV201/sarta_apr08_m140x_iceGHMbaum_waterdrop_desertdust_slabcloud_hg3';

if codeX == 0
  icestr = '_sarta_baran_ice';
  run_sarta.sartacloud_code = code0;
elseif codeX == 1
  icestr = '_sarta_baum_ice';
  run_sarta.sartacloud_code = code1;
else
  error('codeX???')
end
icestr0 = ['cloudy_airs_l1b_ecm'          icestr '.'];
icestr  = ['cloudy_airs_l1b_ecm_wat_jac'  icestr '.'];
icestrP = ['cloudy_airs_l1b_ecm_fixed_pcrtm_only.'];

%%%%%%%%%%%%%%%%%%%%%%%%%

for ixx = 1 : length(iaGlist)
  ix = iaGlist(ixx);
  clear p h hattr pattr prof yymmddgg

  yymmddgg = [yymmdd0 ix];

  ystr = num2str(yymmddgg(1));
  mstr = num2str(yymmddgg(2),'%02d');
  dstr = num2str(yymmddgg(3),'%02d');
  gstr = num2str(yymmddgg(4),'%03d');

  fnameOUT = ['/asl/data/rtprod_airs/' ystr '/' mstr '/' dstr '/JAC/'];
  fnameOUT= [fnameOUT icestr ystr '.' mstr '.' dstr '.' gstr '.rtp'];

  fnameOUT0 = ['/asl/data/rtprod_airs/' ystr '/' mstr '/' dstr '/'];
  fnameOUT0= [fnameOUT0 icestr0 ystr '.' mstr '.' dstr '.' gstr '.rtp'];

  fnameOUTP = ['/asl/data/rtprod_airs/' ystr '/' mstr '/' dstr '/PCRTM/'];
  fnameOUTP= [fnameOUTP icestrP ystr '.' mstr '.' dstr '.' gstr '.rtp'];

  eeP = exist(fnameOUTP);
  ee0 = exist(fnameOUT0);
  ee  = exist(fnameOUT);

  if eeP > 0 & ee0 > 0 & ee == 0
    fprintf(1,' making %s \n',fnameOUT);
    toucher = ['!touch ' fnameOUT];
    eval(toucher)

    fip = mktemp('temp.ip.rtp');
    fop = mktemp('temp.op.rtp');
    frp = mktemp('temp.rp.rtp');
    ugh = mktemp('ugh');
  
    [h,ha,p0,pa] = rtpread(fnameOUT0);
    [h,ha,pP,pa] = rtpread(fnameOUTP);

    if sum(abs(p0.rtime-pP.rtime)) > eps
      fprintf(1,'JOB = %3i rtimes do not match \n',JOB)
      error('rtimes do not match');
    end

    p = pP; 
    p = rmfield(p,'rad_allsky_std'); p = rmfield(p,'rad_clrsky'); p = rmfield(p,'ncol');
    p = rmfield(p,'pcrtm_co2_used'); p = rmfield(p,'overlap'); 
    p = rmfield(p,'pcrtm_iceOD');    p = rmfield(p,'pcrtm_iceDME'); p = rmfield(p,'pcrtm_iceCTOP');
    p = rmfield(p,'pcrtm_lvlODice'); p = rmfield(p,'pcrtm_iceODX');
    p = rmfield(p,'pcrtm_waterOD');    p = rmfield(p,'pcrtm_waterDME'); p = rmfield(p,'pcrtm_waterCTOP');
    p = rmfield(p,'pcrtm_lvlODwater'); p = rmfield(p,'pcrtm_waterODX');
    p = rmfield(p,'rcalc_std');

    p.cfrac12 = p0.cfrac12;
    p.cfrac   = p0.cfrac;
    p.cpsize  = p0.cpsize;
    p.cngwat  = p0.cngwat;
    p.cprtop  = p0.cprtop;
    p.cprbot  = p0.cprbot;
    p.ctype   = p0.ctype;
    p.cfrac2   = p0.cfrac2;
    p.cpsize2  = p0.cpsize2;
    p.cngwat2  = p0.cngwat2;
    p.cprtop2  = p0.cprtop2;
    p.cprbot2  = p0.cprbot2;
    p.ctype2   = p0.ctype2;

    %% doing WATER jacobians
    iw = find(p.ctype  == 101); p.cngwat(iw)  = p0.cngwat(iw) * 1.1;
    iw = find(p.ctype2 == 101); p.cngwat2(iw) = p0.cngwat2(iw) * 1.1;

    [h,ha,p,pa] = rtpadd_emis_wis(h,ha,p,pa);

    rtpwrite(fip,h,ha,p,pa)
    klayerser = ['!' klayers ' fin=' fip ' fout=' fop ' >& ' ugh];
      eval(klayerser);
    try
      [headRX2 hattrR2 profRX2 pattrR2] = rtpread(fop);
    catch me
      me
      fprintf(1,'oops : error running klayers, look at error log %s \n',ugh);
      %keyboard
      error('woof! try again!')
    end

    sarta = run_sarta.sartacloud_code;
    sartaer = ['!' sarta ' fin=' fop ' fout=' frp ' >& ' ugh];
      eval(sartaer);
    try
     [headRX2 hattrR2 profRX2 pattrR2] = rtpread(frp);
    catch me
      me
      fprintf(1,'oops : error running sarta clear, look at error log %s \n',ugh);
      %keyboard
      error('woof! try again!')
    end

    [h2,ha2,p2,pa2] = rtpread(frp);
    rmer = ['!/bin/rm ' fip ' ' fop ' ' frp ' ' ugh]; eval(rmer);

    fnamex = fnameOUT;
    p2x = p;
    p2x.rcalc = p2.rcalc;
    rtpwrite(fnamex,h2,ha2,p2x,pa2)
  else
    disp('new file already there, or important files missing : ')
    fprintf(1,'ee = %3i ee0 = %3i eeP = %3i \n',ee,ee0,eeP);
  end
end
