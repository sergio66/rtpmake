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

theinds = (1 : 2378)';

yymmdd0  = input('Enter [YYYY MM DD] : ');
iaGlist = input('Enter [GranList] : ');

for ixx = 1 : length(iaGlist)
  ix = iaGlist(ixx);
  clear p h hattr pattr prof yymmddgg

  yymmddgg = [yymmdd0 ix];

  ystr = num2str(yymmddgg(1));
  mstr = num2str(yymmddgg(2),'%02d');
  dstr = num2str(yymmddgg(3),'%02d');
  gstr = num2str(yymmddgg(4),'%03d');

  fnameOUT = ['/asl/data/rtprod_airs/' ystr '/' mstr '/' dstr '/'];
  fnameOUT= [fnameOUT 'airs_l1b_era.' ystr '.' mstr '.' dstr '.' gstr '.rtp'];

  eeP = exist(fnameOUT);

  if eeP == 0
    fprintf(1,' making %s \n',fnameOUT);
  
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

h
p

    [h,ha,p,pa] = rtpadd_era_data(h,ha,p,pa);     %%% add on era
  
    [h,ha,p,pa] = rtpadd_emis_DanZhou2(h,ha,p,pa);  
    %% add on emissivity, 99!!! points

    [nemis,efreq,seaemis]=cal_seaemis2(p.satzen,p.wspeed); 
    %% this is better, 19 points
    p.nemis = nemis;
    p.efreq = efreq;
    p.emis  = seaemis;
    p.nrho  = nemis;
    p.rfreq = efreq;
    p.rho   = (1-seaemis)/pi;

    for ii = 1 : length(p.stemp)
      xnemis(ii)   = nemis(ii);
      xefreq(:,ii) = efreq(:,ii);
      N = 1:p.nemis(ii);
      xemis(:,ii) = ...
        interp1(p.efreq(N,ii),p.emis(N,ii),xefreq(:,ii),[],'extrap');
    end

    figure(1)
    %scatter_coast(p.rlon,p.rlat,10,xemis(6,:)); title('emissivity at 954 cm-1')
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

    rtpwrite(fip,h,ha,p,pa);
    disp(' running klayers ....')

    klayerser = ['!' klayers ' fin=' fip ' fout=' fop ' >& ' fugh]; 
      eval(klayerser);
    disp(' running sarta ....')
    sartaer   = ['!' sarta   ' fin=' fop ' fout=' frp ' >& ' fugh]; 
      eval(sartaer);

    [h2,ha2,p2,pa2] = rtpread(frp);

    %h = h2;
    p.rcalc = p2.rcalc;

    rmer = ['!/bin/rm ' fip ' ' fop ' ' frp ' ' fugh]; eval(rmer);

    %%% now save the rtp file!

    fnamex = fnameOUT;
    rtpwrite(fnamex,h,ha,p,pa)
  end
end
