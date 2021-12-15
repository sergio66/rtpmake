%% have to hardcode  YYMMDD before starting, loops over granules

% local running to test
% clustcmd -L clust_make_AIRSL2_90x135_clearrtp_YYMMDD_loopGG.m 001:240
%
% otherwise when happy
% clustcmd -q medium -n 64 clust_make_AIRSL2_90x135_clearrtp_YYMMDD_loopGG.m 001:240
%
% or
% clustcmd -q long_contrib -n 64 clust_make_AIRSL2_90x135_clearrtp_YYMMDD_loopGG.m 001:240

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% March 11, 2011 is my good example
%yymmdd0  = input('Enter [YYYY MM DD] : ');
%iaGlist = input('Enter [GranStart GranList] : ');

yymmdd0  = [2011 12 11];   %% for proposal
yymmdd0  = [2011 11 11];   %% for proposal
yymmdd0  = [2011 10 11];   %% for proposal
yymmdd0  = [2011 09 11];   %% for proposal

yymmdd0  = [2011 03 11];   %% orig

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
addpath /home/sergio/MATLABCODE/L2AIRS_TO_RTP/

theinds = (1 : 2378)';

%%%%%%%%%%%%%%%%%%%%%%%%%
%%% set sarta exec

run_sarta.cumsum = 9999;

codeX = 0; %% use default with A. Baran params
codeX = 1; %% use new     with B. Baum, P. Yang params

code0   = '/asl/packages/sartaV108/BinV201/sarta_apr08_m140_iceaggr_waterdrop_desertdust_slabcloud_hg3_wcon_nte';
code1   = '/home/sergio/SARTA_CLOUDY/BinV201/sarta_apr08_m140x_iceGHMbaum_waterdrop_desertdust_slabcloud_hg3';
codeClr = '/asl/packages/sartaV108/BinV201/sarta_apr08_m140_wcon_nte_nh3';   %%% clear

if codeX == 0
  icestr = '_sarta2_baran_ice';
  run_sarta.sartacloud_code = code0;
elseif codeX == 1
  icestr = '_sarta2_baum_ice';
  run_sarta.sartacloud_code = code1;
else
  error('codeX???')
end
run_sarta.ice_water_separator = 440;   %%%%%%%%% <<<<<<<<<<<<<<<<<<<<<<<<< separate ice and cloud
run_sarta.randomCpsize        = 0;     %%%%%%%%% <<<<<<<<<<<<<<<<<<<<<<<<< use PCRTM sizes

icestr = ['cloudy_airs_l1b_ecm' icestr '.'];

iL2CldType = -1;   %% Kahn
iL2CldType = +1;   %% Joel
%%%%%%%%%%%%%%%%%%%%%%%%%

for ixx = 1 : length(iaGlist)
  yy = yymmdd0(1);  ystr = num2str(yy);
  mm = yymmdd0(2);  mstr = num2str(mm,'%02d');
  dd = yymmdd0(3);  dstr = num2str(dd,'%02d');
  gg = JOB;
  gstr = num2str(JOB,'%03d');

  dirOUT = ['/asl/data/rtprod_airs/' ystr '/' mstr '/' dstr '/AIRS_L2_CLOSURE/'];
  %% airs_cloudNprofL2data.2011.03.11.039_v6.rtp
  if iL2CldType == -1
    fnameOUT= [dirOUT 'allfovSUP_clear_airs_cloudNprofL2data.' ystr '.' mstr '.' dstr '.' gstr '_v6.rtp'];
  elseif iL2CldType == +1
    fnameOUT= [dirOUT 'allfovSUP_clear_airs_cloudNprofL2data_Joel.' ystr '.' mstr '.' dstr '.' gstr '_v6.rtp'];
  end

  eeP = exist(fnameOUT);

  if eeP == 0
    fprintf(1,' making %s \n',fnameOUT);
    %toucher = ['!touch ' fnameOUT];
    %eval(toucher)
    if iL2CldType == -1
      combine_cc_l2ret_v6_cloudsKahn_sergio_90x135_SUP_pfofiles(yy,mm,dd,gg,dirOUT,...
                                                   '/asl/packages/klayersV205/BinV201/klayers_airs',codeClr,-1);
    elseif iL2CldType == +1
      combine_cc_l2ret_v6_cloudsJoel_sergio_90x135_SUP_profiles(yy,mm,dd,gg,dirOUT,...
                                                   '/asl/packages/klayersV205/BinV201/klayers_airs',codeClr,-1);
    end
  else
    fprintf(1,' %s already exists \n',fnameOUT)
  end
end
