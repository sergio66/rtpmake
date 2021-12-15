%% creates an rtp file for ONE granule
%% can be modified for more!

%% JOB 01 - 12 = 2002/01/01 to 2002/12/31 
%% JOB 13 - 24 = 2003/01/01 to 2003/12/31 
%% etc

% local running to test
% clustcmd -L clustcmd_loop_make_merrartp_sergio_cloud_from_ecm.m 001:240
%
% otherwise when happy
% clustcmd -q medium -n 64 -p 2 clustcmd_loop_make_merrartp_sergio_cloud_from_ecm.m 001:240
%
% or
% clustcmd -q long_contrib -n 64 -p 4 clustcmd_loop_make_merrartp_sergio_cloud_from_ecm.m 001:240

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

yymmdd0  = [2011 04 11];   %% to prove to Breno that OpenDap has problems over land; 
                           %% make sure iFix = -1 in rtp_add_merra
yymmdd0  = [2011 03 11];   %% for JGR paper

iaGlist  = JOB;

%% creates an rtp file for ONE granule
%% can be modified for more!

klayers = '/asl/packages/klayers/Bin/klayers_airs';
sarta   = '/asl/packages/sartaV108/Bin/sarta_apr08_m140_wcon_nte';

klayers = '/asl/packages/klayers/Bin/klayers_airs';
sarta   = '/asl/packages/sartaV108/Bin/sarta_apr08_m140_wcon_nte';

klayers = '/asl/packages/klayersV205/BinV201/klayers_airs';
sarta   = '/asl/packages/sartaV108_PGEv6/Bin/sarta_airs_PGEv6_postNov2003';

addpath /asl/matlib/opendap
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

%%% icestr = ['MERRACLD/x1BADTESTcloudy_airs_l1b_merra' icestr '.'];
%%% icestr = ['MERRACLD/x2cloudy_airs_l1b_merra' icestr '.'];     %% this is 2011/04/11 test, with NO fix

icestr = ['MERRACLD/fixJUNK_x2cloudy_airs_l1b_merra' icestr '.'];     %% this is 2011/04/11 test, with    fix

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
    toucher = ['!touch ' fnameOUT];
    eval(toucher)
  
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

    [a,b,c] = sdload(fname);
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

    clrfields = {'SP','SKT','10U','10V','TCC','CI','T','Q','O3'};
    cldfields = {'SP','SKT','10U','10V','TCC','CI','T','Q','O3',...
                 'CC','CIWC','CLWC'};
    [hERA,ha,pERA,pa] = rtpadd_era_data(h,ha,p,pa,cldfields); %%% add on era
    [hECM,ha,pECM,pa] = rtpadd_ecmwf_data(h,ha,p,pa,cldfields); %%% add on ecm

    %% [hMERRA,ha,pMERRA,pa] = rtpadd_merra_cld(h,ha,p,pa); %%% add on merra, original before Feb 2014
    [hMERRA,ha,pMERRA,pa] = rtpadd_merra2(h,ha,p,pa); %%% add on merra, original before Feb 2014

     %% this has 4 less fields : cc,cfrac, wsource and udef

    pMERRA.cfrac = pECM.cfrac;
    [mmmm,nnnn] = size(pMERRA.ptemp);
    for ii = 1 : length(pMERRA.stemp)
      plevsMERRA = log10(pMERRA.plevs(1:pMERRA.nlevs(ii),ii));
      plevsECM = log10(pECM.plevs(1:pECM.nlevs(ii),ii));
      ccECM    = pECM.cc(1:pECM.nlevs(ii),ii);
      xyz      = ones(mmmm,1)*0;
      ccMERRA  = interp1(plevsECM,ccECM,plevsMERRA,[],'extrap');
      ccMERRA(ccMERRA < 0) = 0;
      ccMERRA(ccMERRA > 1) = 1;
      xyz(1:length(ccMERRA)) = ccMERRA;
      pMERRA.cc(:,ii) = xyz;
    end
    p = pMERRA;
    h = hMERRA;
    p0 = p;
  
    %[h,ha,p,pa] = rtpadd_emis_DanZhou2(h,ha,p,pa);  
    %p = Prof_add_emis(p,yymmddgg(1),yymmddgg(2),yymmddgg(3));  %% broken crap by whoever
    %p = rtpadd_emis_DanZhou(h,ha,p,pa);   %% lso totally broken crap
    [h,ha,p,pa] = rtpadd_emis_wis(h,ha,p,pa);

%{
    %figure(1)
    %scatter_coast(p.rlon,p.rlat,10,p.nemis); 

    TzAVG = sum(pMERRA.ptemp)/42;
    scatter_coast(pECM.rlon,pECM.rlat,30,TzAVG)

    lala = find(pMERRA.landfrac > 1-eps);
    plot(lala,TzAVG(lala))
    plot(pMERRA.ptemp(:,2300),pMERRA.plevs(:,2300))
    [pMERRA.ptemp(:,2300),pMERRA.plevs(:,2300)]
       [pMERRA.nlevs(2300) pMERRA.stemp(2300) pMERRA.spres(2300)]

    [pMERRA.ptemp(:,lala(1392)),pMERRA.plevs(:,lala(1392))]
       [pMERRA.nlevs(lala(1392)) pMERRA.stemp(lala(1392)) pMERRA.spres(lala(1392))]

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    error('ooo')
%}

    [p2] = driver_sarta_cloud_rtp(h,ha,p,pa,run_sarta);

    fnamex = fnameOUT;
    p2x = p2;
    [h,ha,p2x,pa] = rtptrim_sartacloud(h,ha,p2,pa);
    rtpwrite(fnamex,h,ha,p2x,pa)

  else
    fprintf(1,'%s already exists \n',fnameOUT)
  end
end
