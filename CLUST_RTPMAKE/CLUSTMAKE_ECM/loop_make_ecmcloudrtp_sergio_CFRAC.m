%% creates an rtp file for ONE granule
%% can be modified for more!

%% code is same as loop_make_eracloudrtp_sergio.m except that
%% a) it only does 6 chans 
%% b) changes the cfrac from default level2slab_clouds.m to
%%    level2slab_clouds_cfrac.m
%%    (see /home/sergio/MATLABCODE/CLOUD_ECMWF_ERA/PACKAGE)

klayers = '/asl/packages/klayers/Bin/klayers_airs';
sarta   = '/asl/packages/sartaV108/Bin/sarta_apr08_m140_wcon_nte';

sarta = '/asl/packages/sartaV108/BinV201/';
sarta = [sarta 'sarta_apr08_m140_iceaggr_waterdrop_desertdust_slabcloud_hg3_wcon_nte'];

addpath /asl/matlab/science
addpath /asl/matlab/rtptoolsV201
addpath /asl/matlab/h4tools/
addpath /asl/matlab/rtptools/
addpath /asl/matlab/gribtools/
addpath /asl/matlab/airs/readers/
addpath /home/sergio/MATLABCODE/CLOUD_ECMWF_ERA/PACKAGE
addpath /home/sergio/MATLABCODE/PLOTTER
addpath /home/sergio/MATLABCODE

iaChanList = [1:2378];
iaChanList = [359 445 903 1291 1557 1614];  %% for MATLABCODE/RATES_CLOUD

%%%%%%%%%%%%%%%%%%%%%%%%%
rCNGWATiceMult = 1.0;
rCPSIZEiceMult = 1.0;
%%%%%%%%%%%%%%%%%%%%%%%%%

rCfracType = 0.5; %% use constant |x| for cfrac1,cfrac2, cfrac12
rCfracType = 0.99; %% use constant |x| for cfrac1,cfrac2, cfrac12
rCfracType = -1;  %% use wghted average      of sum(CC) in cfrac/cfrac2/cfrac12
                  %%                         almost rCfracType = 0, but uses CC (91xn)
rCfracType = -2;  %% use wghted average of sum(CC) (91xn), then SAME cfrac/cfrac2/cfrac12
rCfracType = 0;   %% default, use Scott's methods which uses TCC (1xn)

rCfracType = -2;  %% use wghted average of sum(CC) (91xn), then SAME cfrac/cfrac2/cfrac12
rCfracType = -1;  %% use wghted average      of sum(CC) in cfrac/cfrac2/cfrac12
                  %%                         almost rCfracType = 0, but uses CC (91xn)

%%%%%%%%%%%%%%%%%%%%%%%%%
rCNGWATiceMult = 5.0; rCPSIZEiceMult = 1.0;
rCfracType = 0;   %% default, use Scott's methods which uses TCC (1xn)
rCfracType = -1;  %% use wghted average      of ccfrac in cfrac/cfrac2/cfrac12
                  %%                         almost rCfracType = 0, but uses CC (91xn)
rCNGWATiceMult = 10.0; rCPSIZEiceMult = 0.5;
%%%%%%%%%%%%%%%%%%%%%%%%%

addpath /home/sergio/MATLABCODE

lala = +1;
disp('Init : get the list of granules (loop) : ');
disp('          (+1) : prefer to do one at a time')
disp('          (+2) : enter in [yy mm dd] and then a list of [gg] (+N for 1:N) ');
lala = input('Enter (1) or (2) : ');

iCnt = 0;
if lala == 1
  while lala > 0
    iOK = -1;
    while iOK < 0
      blah = input('Enter [YY MM DD GG] : ');
      if length(blah) ~= 4
        disp('    oops : need 4 elements [YY MM DD GG] retry : ');
      else
        iOK = +1;
      end
    end
    iCnt = iCnt + 1;
    loop(iCnt).yymmddgg = blah;
    lala = input('  another? (-1 to quit) : ');
  end
elseif lala == 2
  iOK = -1;
  while iOK < 0
    blah = input('Enter [YY MM DD] : ');
    if length(blah) ~= 3
      disp('    oops : need 3 elements [YY MM DD] retry : ');
    else
      iOK = +1;
    end
  end
  yymmdd = blah;
  lalax  = input('Enter in +N (for gg = 1 : N) or -1 (for gg list) : ');
  if lalax > 0
    for iCnt = 1 : abs(lalax)
      loop(iCnt).yymmddgg = [[yymmdd] iCnt];
    end
  else
    gglist = input('Enter gglist [g1 g2 ..... gN] : ');
    for iCnt = 1 : length(gglist)
      loop(iCnt).yymmddgg = [[yymmdd] gglist(iCnt)];
    end
  end
end

for ix = 1 : iCnt
  yymmddgg = loop(ix).yymmddgg; 
  fprintf(1,'%5i : %4i/%2i/%2i  gran = %3i \n',ix,yymmddgg);
end
  
for ix = 1 : iCnt
  disp(' ')
  clear p h hattr pattr prof yymmddgg
  yymmddgg = loop(ix).yymmddgg;
  ystr = num2str(yymmddgg(1));
  mstr = num2str(yymmddgg(2),'%02d');
  dstr = num2str(yymmddgg(3),'%02d');
  gstr = num2str(yymmddgg(4),'%03d');

  fnamex = ['/strowdataN/s3/robinso2/PUYEHUE/CLOUDS_6chans/'];

  if (abs(rCNGWATiceMult-1) < eps) & (abs(rCPSIZEiceMult-1) < eps) 
    %% DEFAULT : no ice amount or ice size multiplier
    if rCfracType == 0
      fnamex = [fnamex 'TYPE0/ecmcld_type0_allfov_'];
    elseif rCfracType == -1
     fnamex = [fnamex 'TYPE1/ecmcld_type1_allfov_'];
    elseif rCfracType == -2
      fnamex = [fnamex 'TYPE2/ecmcld_type2_allfov_'];
    elseif rCfracType > 0  & rCfracType <= 1
      fnamex = [fnamex 'TYPEX/ecmcld_typeX_allfov_'];
    end
  elseif (abs(rCNGWATiceMult-1) > eps) & (abs(rCPSIZEiceMult-1) < eps) 
    %% have put ice amount multiplier, but no ice size multiplier
    if rCfracType == 0
      fnamex = [fnamex 'TYPEI_rCfracType=0_x' num2str(rCNGWATiceMult)];
      fnamex = [fnamex '/ecmcld_type0_allfov_'];
    elseif rCfracType == -1
      fnamex = [fnamex 'TYPEI_rCfracType=-1_x' num2str(rCNGWATiceMult)];
      fnamex = [fnamex '/ecmcld_type0_allfov_'];
    elseif rCfracType == -2
      fnamex = [fnamex 'TYPEI_rCfracType=-2_x' num2str(rCNGWATiceMult)];
      fnamex = [fnamex '/ecmcld_type0_allfov_'];
    end
  elseif (abs(rCNGWATiceMult-1) > eps) & (abs(rCPSIZEiceMult-1) > eps) 
    %% have put ice amount multiplier and ice size multiplier
    if rCfracType == 0
      fnamex = [fnamex 'TYPEIS_rCfracType=0_x' num2str(rCNGWATiceMult)];
      fnamex = [fnamex '/ecmcld_type0_allfov_'];
    elseif rCfracType == -1
      fnamex = [fnamex 'TYPEIS_rCfracType=-1_x' num2str(rCNGWATiceMult)];
      fnamex = [fnamex '/ecmcld_type0_allfov_'];
    elseif rCfracType == -2
      fnamex = [fnamex 'TYPEIS_rCfracType=-2_x' num2str(rCNGWATiceMult)];
      fnamex = [fnamex '/ecmcld_type0_allfov_'];
    else
      error('dude!!l which dir to put stuff into????')
    end
  end

  fnamex = [fnamex ystr '_' mstr '_' dstr '_' gstr '.rtp'];
  ee = exist(fnamex);
  if ee > 0
    fprintf(1,' %s already exists ... processing next \n',fnamex);
  else
    fprintf(1,' making %s \n',fnamex);
    toucher = ['!touch ' fnamex]; eval(toucher)
  
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
    filename = ...
      ['/strowdataN/data/airs/Aqua_AIRS_Level1/AIRIBRAD.005/' ystr '/'];
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

    % [a,b,c] = sdload_quiet(fname);
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

    clrfields = {'SP','SKT','10U','10V','TCC','CI','T','Q','O3'};
    cldfields = {'SP','SKT','10U','10V','TCC','CI','T','Q','O3',...
                 'CC','CIWC','CLWC'};
    %% TCC is total cloud cover (1 x nprofs),     gets stuffed into p.cfrac
    %%  CC is cloud cover,       (nlevs x nprofs) field, used in 
    %%                           level2slab_clouds_CFRAC

    [h,ha,p,pa] = rtpadd_ecmwf_data(h,ha,p,pa,cldfields); %%% add on ecm

    %% now subset p as needed
    h.nchan = length(iaChanList);
    h.ichan = iaChanList';
    h.vchan = f(iaChanList);

    p.calflag = p.calflag(iaChanList,:);
    p.robs1   = p.robs1(iaChanList,:);

    h0 = h;
    p0 = p;   %% so if loop crashes, you have orig profile p0 to debug!

    % ccx = quick_cloud_frac(p0);
    % [hA,pA] = level2slab_clouds_CFRAC(h0,p0,0);
    % [hB,pB] = level2slab_clouds_CFRAC(h0,p0,-1);
    % [hC,pC] = level2slab_clouds_CFRAC(h0,p0,0.25);
    % [hD,pD] = level2slab_clouds_CFRAC(h0,p0,-2);

    %if (abs(rCNGWATiceMult-1) < eps)
    %  [h,p] = level2slab_clouds_CFRAC(h0,p0,rCfracType,1.0);
    %else
    %  [h,p] = level2slab_clouds_CFRAC(h0,p0,rCfracType,rCNGWATiceMult);
    %end
    [h,p] = level2slab_clouds_CFRAC(h0,p0,rCfracType,rCNGWATiceMult,rCPSIZEiceMult);

    ciwc = p.ciwc;
    clwc = p.clwc;

    %% add on emissivity, 99!!! points
    [h,ha,p,pa] = rtpadd_emis_DanZhou2(h,ha,p,pa);  

    %% this is better, 19 points
    % p.nemis = nemis;
    % p.efreq = efreq;
    % p.emis  = seaemis;
    % p.nrho  = nemis;
    % p.rfreq = efreq;
    % p.rho   = (1-seaemis)/pi;
    clear nemis efreq seaemis xnemis xefreq xemis
    [nemis,efreq,seaemis]=cal_seaemis2(p.satzen,p.wspeed); 

    for ii = 1 : length(p.stemp)
      xnemis(ii)   = nemis(ii);
      xefreq(:,ii) = efreq(:,ii);
      N = 1:p.nemis(ii);
      xemis(:,ii) = ...
        interp1(p.efreq(N,ii),p.emis(N,ii),xefreq(:,ii),[],'extrap');
    end

    % figure(1)
    % scatter_coast(p.rlon,p.rlat,10,xemis(6,:)); 
    % title('emissivity at 954 cm-1')

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

    oldrtpwrite(fip,h,ha,p,pa);
    disp(' running klayers ....')
    klayerser = ['!' klayers ' fin=' fip ' fout=' fop ' >& ' fugh]; 
      eval(klayerser);
    disp(' running sarta ....')
    sartaer   = ['!' sarta   ' fin=' fop ' fout=' frp]; 
      eval(sartaer);

    [h2,ha2,p2,pa2] = oldrtpread(frp);

    %h = h2;
    p.rcalc = p2.rcalc;

    %{
    dodo = find(h.ichan == 1291);
    if length(dodo) == 1
      dBT = 200 : 320;
      tobs_1231 = rad2bt(1231,p.robs1(dodo,:));
      tcal_1231 = rad2bt(1231,p.rcalc(dodo,:));
      nObs = hist(tobs_1231,dBT); 
      nCal = hist(tcal_1231,dBT); 
    figure(3)
      plot(dBT,nObs,dBT,nCal)
    figure(4)
      scatter_coast(p.rlon,p.rlat,10,tobs_1231-tcal_1231); 
      title('bias : obs-cal : BT(1231 cm-1)')
      caxis([-10 +10]); colorbar
    end
    %}

    rmer = ['!/bin/rm ' fip ' ' fop ' ' frp ' ' fugh]; eval(rmer);

    %%% now save the rtp file!
    oldrtpwrite(fnamex,h,ha,p,pa)
  end
end