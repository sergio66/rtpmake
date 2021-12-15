%% have to hardcode  YYMMDD before starting, loops over granules

% local running to test
% clustcmd -L clust_make_eraclearrtp_sergio_sarta_YYMMDD_loopGG.m 001:240
%
% otherwise when happy
% clustcmd -q medium -n 64 -p 4 clust_make_eraclearrtp_sergio_sarta_YYMMDD_loopGG.m 001:240
%
% or
% clustcmd -q long_contrib -n 64 -p 4 clust_make_eraclearrtp_sergio_sarta_YYMMDD_loopGG.m 001:240

%% run with
%% sbatch --array=N1-N2 --output='testslurm' sergio_matlab_jobB.sbatch
%% N1 = 1, N2 = number of files to be processed

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% March 11, 2011 is my good example
%yymmdd0  = input('Enter [YYYY MM DD] : ');
%iaGlist = input('Enter [GranStart GranList] : ');

yymmdd0  = [2013 08 26];   %% Larrabee wants this
yymmdd0  = [2013 08 27];   %% Larrabee wants this

yymmdd0  = [2013 06 09];   %% for Steve testing

yymmdd0 = [2016 01 018]    %% is fill_era messed up compared to rtpadd_era_data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%JOB = str2num(getenv('SLURM_ARRAY_TASK_ID'));
JOB = 180
iaGlist  = JOB;

%% creates an rtp file for ONE granule
%% can be modified for more!

%klayers = '/asl/packages/klayers/Bin/klayers_airs';
%sarta   = '/asl/packages/sartaV108/Bin/sarta_apr08_m140_wcon_nte';

%klayers = '/asl/packages/klayers/Bin/klayers_airs';
%sarta   = '/asl/packages/sartaV108/Bin/sarta_apr08_m140_wcon_nte';

klayers = '/asl/packages/klayersV205/BinV201/klayers_airs';
sarta   = '/asl/packages/sartaV108_PGEv6/Bin/sarta_airs_PGEv6_postNov2003';

addpath /asl/matlib/aslutil
addpath /asl/matlib/science
addpath /asl/matlib/rtptools
addpath /asl/matlib/h4tools/
addpath /asl/matlib/rtptools/
addpath /asl/matlib/gribtools/

theinds = (1 : 2378)';

run_sarta.klayers_code    = klayers;
run_sarta.sartaclear_code = sarta;

%%%%%%%%%%%%%%%%%%%%%%%%%
%%% set sarta exec

icestr = ['allfov_era_'];

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

    cldfields = {'SP','SKT','10U','10V','TCC','CI','T','Q','O3',...
                 'CC','CIWC','CLWC'};
    clrfields = {'SP','SKT','10U','10V','TCC','CI','T','Q','O3'};

    [h,ha,p,pa] = rtpadd_era_data(h,ha,p,pa,clrfields); %%% add on era using GRIB files which we dont have now
    [p,h] = fill_era(p,h);                              %%% add on era the NEW WAY
error('slgkglks')

    p0 = p;
  
    %[h,ha,p,pa] = rtpadd_emis_DanZhou2(h,ha,p,pa);  
    %p = Prof_add_emis(p,yymmddgg(1),yymmddgg(2),yymmddgg(3));  %% broken crap by whoever
    %p = rtpadd_emis_DanZhou(h,ha,p,pa);   %% lso totally broken crap
    [h,ha,p,pa] = rtpadd_emis_wis(h,ha,p,pa);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [p2] = driver_sarta_clear_rtp(h,ha,p,pa,run_sarta);

    fnamex = fnameOUT;
    p.rcalc = p2.rcalc;
    rtpwrite(fnamex,h,ha,p,pa)

    figure(1); scatter(p.rlon,p.rlat,20,rad2bt(1231,p.robs1(1291,:))); title('BT1231 obs')
      colorbar; cx = caxis;
    figure(2); scatter(p.rlon,p.rlat,20,rad2bt(1231,p.rcalc(1291,:))); title('BT1231 cal')
      caxis(cx); colorbar
  else
    fprintf(1,'%s already exists \n',fnameOUT)
  end
end
