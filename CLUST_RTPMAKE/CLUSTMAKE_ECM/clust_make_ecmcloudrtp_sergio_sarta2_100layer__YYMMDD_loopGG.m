%% have to hardcode  YYMMDD before starting, loops over granules
%%
%% same as clust_make_ecmcloudrtp_sergio_sarta_YYMMDD_loopGG.m except we SEPARATE ice N water clouds
%% using run_sarta.ice_water_separator = 440; 
%% same as clust_make_ecmcloudrtp_sergio_sarta2_YYMMDD_loopGG.m except we use 100 layer code

% local running to test
% clustcmd -L clust_make_ecmcloudrtp_sergio_sarta2_100layer__YYMMDD_loopGG.m 001:240
%
% otherwise when happy
% clustcmd -q medium -n 32 -p 4 clust_make_ecmcloudrtp_sergio_sarta2_100layer__YYMMDD_loopGG.m 001:240
%
% or
% clustcmd -q long_contrib -n 32 -p 4 clust_make_ecmcloudrtp_sergio_sarta2_100layer__YYMMDD_loopGG.m 001:240

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% March 11, 2011 is my good example
%yymmdd0  = input('Enter [YYYY MM DD] : ');
%iaGlist = input('Enter [GranStart GranList] : ');

yymmdd0  = [2011 03 11];
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

klayers = '/home/sergio/klayersV205/BinV201/klayers_airs_x_testmeCLOUDversion';
sarta   = '/home/sergio/SARTA_CLOUDY/BinV201/sarta_apr08_m140x_iceGHMbaum_waterdrop_desertdust_slabcloud_hg3_100layerNEW';

addpath /asl/matlib/aslutil
addpath /asl/matlib/science
addpath /asl/matlib/rtptools
addpath /asl/matlib/h4tools/
addpath /asl/matlib/rtptools/
addpath /asl/matlib/gribtools/
addpath /home/sergio/MATLABCODE/matlib/clouds/sarta
addpath /home/sergio/MATLABCODE

theinds = (1 : 2378)';

%%%%%%%%%%%%%%%%%%%%%%%%%
%%% set sarta exec

run_sarta.clear = +1;
run_sarta.cloud = +1;

codeX = 0; %% use default with A. Baran params
codeX = 1; %% use new     with B. Baum, P. Yang params

code0 = '/asl/packages/sartaV108/BinV201/sarta_apr08_m140_iceaggr_waterdrop_desertdust_slabcloud_hg3_wcon_nte';
code1 = '/home/sergio/SARTA_CLOUDY/BinV201/sarta_apr08_m140x_iceGHMbaum_waterdrop_desertdust_slabcloud_hg3';

code0 = 'BOO not done';
code1 = '/home/sergio/SARTA_CLOUDY/BinV201/sarta_apr08_m140x_iceGHMbaum_waterdrop_desertdust_slabcloud_hg3_100layerNEW';

if codeX == 0
  icestr = '_sarta2_baran_ice_100layercld';
  run_sarta.sartacloud_code = code0;
elseif codeX == 1
  icestr = '_sarta2_baum_ice_100layercld';
  run_sarta.sartacloud_code = code1;
else
  error('codeX???')
end

%%% run_sarta.cumsum = 9999;
%run_sarta.randomCpsize        = 0;     %%%%%%%%% <<<<<<<<<<<<<<<<<<<<<<<<< use PCRTM sizes
run_sarta.ncol                = 25;    %%%%%%%%% <<<<< number of subcolumns per pixel
run_sarta.cfrac               = 1;     %%%%%%%%% <<<<< set cldfrac = 1, instead of random
run_sarta.ice_water_separator = 440;   %%%%%%%%% <<<<<<<<<<<<<<<<<<<<<<<<< separate ice and water cloud, ala PCRTM
run_sarta.ice_water_separator = -1;    %%%%%%%%% <<<<<<<<<<<<<<<<<<<<<<<<< overlap ice and water cloud, ala 2 slabs

if run_sarta.ice_water_separator == -1
  icestr = [icestr '_overlapIW'];
else
  icestr = [icestr '_separateIW_at_mb=' num2str(run_sarta.ice_water_separator)];
end

icestr = ['cloudy_airs_l1b_ecm' icestr '.'];

%%%%%%%%%%%%%%%%%%%%%%%%%

for ixx = 1 : length(iaGlist)
  ix = iaGlist(ixx);
  clear p h hattr pattr prof yymmddgg

  yymmddgg = [yymmdd0 ix];

  ystr = num2str(yymmddgg(1));
  mstr = num2str(yymmddgg(2),'%02d');
  dstr = num2str(yymmddgg(3),'%02d');
  gstr = num2str(yymmddgg(4),'%03d');

  fnameOUT = ['/asl/data/rtprod_airs/' ystr '/' mstr '/' dstr '/'];
  fnameOUT= [fnameOUT icestr ystr '.' mstr '.' dstr '.' gstr '.rtp'];

  eeP = exist(fnameOUT);

  if eeP == 0
    fprintf(1,' making %s \n',fnameOUT);
    toucher = ['!/bin/touch ' fnameOUT];
    eval(toucher);
  
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
      return
    end

    [a,b,c] = sdload_quiet(fname);
    p.rlat = [a.Latitude(:)'];
    p.rlon = [a.Longitude(:)'];
    p.rtime = [a.Time(:)'];
  
    [meantime, f, prof] = readl1b_all(fname);  
    p = prof;

    p.pobs = zeros(size(p.solazi));
    p.upwell = ones(size(p.solazi));
    %p.irinst = AIRSinst*ones(1,nobs);
    %p.findex = grannum*ones(1,nobs);

    plot(p.rlon,p.rlat,'.')

    pa = {{'profiles','rtime','seconds since 1993'}};
    ha = {{'header','hdf file',filename}};

    h.pfields=5; % (1=prof + 4=IRobs);


    h.nchan = length(theinds);
    h.ichan = theinds;;
    h.vchan = f(h.ichan);;

h
p

    clrfields = {'SP','SKT','10U','10V','TCC','CI','T','Q','O3'};
    cldfields = {'SP','SKT','10U','10V','TCC','CI','T','Q','O3',...
                 'CC','CIWC','CLWC'};
    [h,ha,p,pa] = rtpadd_ecmwf_data(h,ha,p,pa,cldfields); %%% add on ecm
    p0 = p;
  
    %[h,ha,p,pa] = rtpadd_emis_DanZhou2(h,ha,p,pa);  
    %p = Prof_add_emis(p,yymmddgg(1),yymmddgg(2),yymmddgg(3));  %% broken crap by whoever
    %p = rtpadd_emis_DanZhou(h,ha,p,pa);   %% lso totally broken crap
    [h,ha,p,pa] = rtpadd_emis_wis(h,ha,p,pa);

    %figure(1)
    %scatter_coast(p.rlon,p.rlat,10,p.nemis); 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [p2] = driver_sarta_cloud100layer_rtp(h,ha,p,pa,run_sarta);

    fnamex = fnameOUT;
    [h,ha,p2x,pa] = rtptrim_sartacloud(h,ha,p2,pa);
    rtpwrite(fnamex,h,ha,p2x,pa)
  else
    fprintf(1,' %s already exists \n',fnameOUT)
  end
end
