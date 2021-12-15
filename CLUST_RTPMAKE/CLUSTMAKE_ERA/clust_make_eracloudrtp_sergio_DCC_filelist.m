%% run with
%% sbatch --array=N1-N2 --output='testslurm' sergio_matlab_jobB.sbatch
%% N1 = 1, N2 = number of files to be processed

addpath /asl/matlab2012/airs/readers
addpath /asl/matlib/aslutil
addpath /asl/matlib/science
addpath /asl/matlib/rtptools
addpath /asl/matlib/h4tools/
addpath /asl/matlib/rtptools/
addpath /asl/matlib/gribtools/
addpath /home/sergio/MATLABCODE/matlib/clouds/sarta
addpath /home/sergio/MATLABCODE

%addpath /home/strow/cress/Work/Rtp
% addpath /home/strow/Matlab/Grib     WARNING /home/strow/Matlab/Grib/rtpadd_grib_data.m DIFFERENT than /asl/matlib/gribtools/rtpadd_era_data.m
% addpath /home/sergio/MATLABCODE/CRIS_HiRes             %% for sergio_fill_ecmwf
% addpath /home/strow/Git/rtp_prod2/grib                  %% for fill_ecm
addpath /asl/packages/rtp_prod2/grib


JOB = str2num(getenv('SLURM_ARRAY_TASK_ID'));
%JOB = 6

filelist = 'dcc_days.txt';
thefilelist = load(filelist);
thefilelist = thefilelist(JOB,:);

yymmdd0  = thefilelist(1:3); %% YY MM DD

daysINmonth = [31 28 31 30 31 30 31 31 30 31 30 31];
if mod(yymmdd0(1),4) == 0
  daysINmonth(2) = 29;
end

doy = sum(daysINmonth(1:yymmdd0(2)-1)) + yymmdd0(3);

dccINfile = ['/asl/rtp_lustre/rtp_airxbcal_v5/' num2str(yymmdd0(1)) '/dcc/era_airxbcal_day' num2str(doy,'%03d') '_dcc.rtp'];
if ~exist(dccINfile)
  fprintf(1,'%s DNE \n',dccINfile);
  return
end

[h0,ha0,p0,pa0] = rtpread(dccINfile);

grans = unique(p0.findex);
for ix = 1 : length(grans)
  oo = find(p0.findex == grans(ix));
  numfovs(ix) = length(oo);

  bonk = p0.atrack(oo);
  atrack1(ix) = mode(bonk);     %% pull everything along this atrack
  zonk = find(bonk == atrack1(ix));
  rlat_atrack1(ix)   = mean(p0.rlat(oo(zonk)));  
  numtimes_atrack1(ix) = length(find(bonk == atrack1(ix)));

  bonk = p0.xtrack(oo);
  xtrack1(ix) = mode(bonk);     %% pull everything along this atrack
  zonk = find(bonk == xtrack1(ix));
  rlon_xtrack1(ix)   = mean(p0.rlon(oo(zonk)));  
  numtimes_xtrack1(ix) = length(find(bonk == xtrack1(ix)));

  if length(oo) > 100
      bonk = p0.atrack(oo);
      [m,values,count] = mode_order(bonk);      
      atrack2(ix) = values(2);
        zonk = find(bonk == atrack2(ix));
        rlat_atrack2(ix)   = mean(p0.rlat(oo(zonk)));  
        numtimes_atrack2(ix) = length(find(bonk == atrack2(ix)));
      atrack3(ix) = values(3);
        zonk = find(bonk == atrack3(ix));
        rlat_atrack3(ix)   = mean(p0.rlat(oo(zonk)));  
        numtimes_atrack3(ix) = length(find(bonk == atrack3(ix)));

      bonk = p0.xtrack(oo);
      [m,values,count] = mode_order(bonk);      
      xtrack2(ix) = values(2);
        zonk = find(bonk == xtrack2(ix));
        rlon_xtrack2(ix)   = mean(p0.rlon(oo(zonk)));  
        numtimes_xtrack2(ix) = length(find(bonk == xtrack2(ix)));
      xtrack3(ix) = values(3);
        zonk = find(bonk == xtrack3(ix));
        rlon_xtrack3(ix)   = mean(p0.rlon(oo(zonk)));  
        numtimes_xtrack3(ix) = length(find(bonk == xtrack3(ix)));

  else
    atrack2(ix) = atrack1(ix); atrack3(ix) = atrack1(ix); rlat_atrack2(ix) = rlat_atrack1(ix); rlat_atrack3(ix) = rlat_atrack1(ix);
    xtrack2(ix) = xtrack1(ix); xtrack3(ix) = xtrack1(ix); rlon_xtrack2(ix) = rlon_xtrack1(ix); rlon_xtrack3(ix) = rlon_xtrack1(ix);
  end
  
  figure(1); plot(p0.rlon(oo),p0.rlat(oo),'.'); title(num2str(length(oo)));
    ax = axis; line([ax(1) ax(2)],[rlat_atrack1(ix) rlat_atrack1(ix)],'color','r','linewidth',2)
    ax = axis; line([ax(1) ax(2)],[rlat_atrack2(ix) rlat_atrack2(ix)],'color','m')
    ax = axis; line([ax(1) ax(2)],[rlat_atrack3(ix) rlat_atrack3(ix)],'color','m')
    ax = axis; line([rlon_xtrack1(ix) rlon_xtrack1(ix)],[ax(3) ax(4)],'color','b','linewidth',2)
    ax = axis; line([rlon_xtrack2(ix) rlon_xtrack2(ix)],[ax(3) ax(4)],'color','c')
    ax = axis; line([rlon_xtrack3(ix) rlon_xtrack3(ix)],[ax(3) ax(4)],'color','c')
  figure(2); plot(p0.xtrack(oo),p0.atrack(oo),'.'); title(num2str(length(oo)));
    ax = axis; line([ax(1) ax(2)],[atrack1(ix) atrack1(ix)],'color','r','linewidth',2)
    ax = axis; line([ax(1) ax(2)],[atrack2(ix) atrack2(ix)],'color','m')
    ax = axis; line([ax(1) ax(2)],[atrack3(ix) atrack3(ix)],'color','m')    
    ax = axis; line([xtrack1(ix) xtrack1(ix)],[ax(3) ax(4)],'color','b','linewidth',2)
    ax = axis; line([xtrack2(ix) xtrack2(ix)],[ax(3) ax(4)],'color','c')
    ax = axis; line([xtrack3(ix) xtrack3(ix)],[ax(3) ax(4)],'color','c')        
  
  %disp('ret'); pause
  pause(0.1)
end
[grans; numfovs]
figure(1); scatter(grans,numfovs,60,numtimes_atrack1,'filled'); colorbar;
  xlabel('Gran number'); ylabel('numfovs'); title('Most common atrack count'); colormap jet
figure(2); scatter(grans,numfovs,60,numtimes_xtrack1,'filled'); colorbar;
  xlabel('Gran number'); ylabel('numfovs'); title('Most common xtrack count'); colormap jet
figure(3); plot(xtrack1,numtimes_xtrack1,'bo',atrack1,numtimes_atrack1,'rx')

[Y,I] = sort(numfovs,'descend');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% creates an rtp file for ONE granule
%% can be modified for more!

klayers = '/asl/packages/klayers/Bin/klayers_airs';
sarta   = '/asl/packages/sartaV108/Bin/sarta_apr08_m140_wcon_nte';

klayers = '/asl/packages/klayers/Bin/klayers_airs';
sarta   = '/asl/packages/sartaV108/Bin/sarta_apr08_m140_wcon_nte';

klayers = '/asl/packages/klayersV205/BinV201/klayers_airs';
sarta   = '/asl/packages/sartaV108_PGEv6/Bin/sarta_airs_PGEv6_postNov2003';


theinds = (1 : 2378)';

%%%%%%%%%%%%%%%%%%%%%%%%%
%%% set sarta exec

run_sarta.clear = +1;
run_sarta.cloud = +1;
run_sarta.cumsum = -1;    %% this is "closer" to MRO but since cliuds are at centroid, does not do too well with DCC
run_sarta.cumsum = 9999;  %% larrabee likes this, puts clouds high so does well for DCC

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

%icestr = ['NEWLANDFRAC/cloudy_airs_l1b_era' icestr '.'];
%icestr = ['cloudy_airs_l1b_era' icestr '.'];
icestr = ['cloudy_airs_l1b_era' icestr];

%%%%%%%%%%%%%%%%%%%%%%%%%

ystr = num2str(yymmdd0(1));
fdirOUT = ['/asl/s1/sergio/rtp/rtp_airibrad_v5/dcc/' ystr '/'];

if ~exist(fdirOUT)
  mker = ['!/bin/mkdir -p ' fdirOUT];
  fprintf(1,'mker = %s \n',mker);
  eval(mker)
end
  
fnameOUT= [fdirOUT icestr '_era_airxbcal_day' num2str(doy,'%03d') '_dcc.rtp'];

eeP = exist(fnameOUT);
if eeP > 0
  fprintf(1,'%s already exists \n',fnameOUT)
  return
end

if eeP == 0
  fprintf(1,' making %s \n',fnameOUT);
  toucher = ['!touch ' fnameOUT];
  eval(toucher)
end

%%%%%%%%%%%%%%%%%%%%%%%%%
hall = [];
pall = [];

hjunk = [];
pjunk = [];
for ixx = 1 : length(grans)
  ix = grans(ixx);
  clear p h hattr pattr prof yymmddgg

  yymmddgg = [yymmdd0 ix];

  ystr = num2str(yymmddgg(1));
  mstr = num2str(yymmddgg(2),'%02d');
  dstr = num2str(yymmddgg(3),'%02d');
  gstr = num2str(yymmddgg(4),'%03d');
  
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

  %[meantime, f, prof] = readl1b_all(fname);  %% has the old rtime 1993
  [meantime, f, prof] = xreadl1b_all_quiet(fname);  %% has the new rtime 1958
  p = prof;

  figure(3); clf; scatter(p.rlon,p.rlat,10,real(rad2bt(1231,p.robs1(1291,:)))); title('BT1231 obs'); colormap jet; colorbar
    
  p.pobs = zeros(size(p.solazi));
  p.upwell = ones(size(p.solazi));
  %p.irinst = AIRSinst*ones(1,nobs);
  %p.findex = grannum*ones(1,nobs);

  pa = {{'profiles','rtime','seconds since 1993'}};
  ha = {{'header','hdf file',filename}};

  h.pfields=5; % (1=prof + 4=IRobs);

  h.nchan = length(theinds);
  h.ichan = theinds;;
  h.vchan = f(h.ichan);;

  %%% this is NEW
  p.landfrac_fromL1B = p.landfrac;
  p.salti_fromL1B = p.salti;
  [salti, landfrac] = usgs_deg10_dem(p.rlat, p.rlon);
  p.landfrac = landfrac;
  p.salti    = salti;

  clrfields = {'SP','SKT','10U','10V','TCC','CI','T','Q','O3'};
  cldfields = {'SP','SKT','10U','10V','TCC','CI','T','Q','O3',...
               'CC','CIWC','CLWC'};

  %[h,ha,p,pa] = rtpadd_era_data(h,ha,p,pa,cldfields); %%% add on era
  [p,h] = fill_era(p,h);
    
  %[h,ha,p,pa] = rtpadd_emis_DanZhou2(h,ha,p,pa);
  %p = Prof_add_emis(p,yymmddgg(1),yymmddgg(2),yymmddgg(3));  %% broken crap by whoever
  %p = rtpadd_emis_DanZhou(h,ha,p,pa);   %% lso totally broken crap
  %[h,ha,p,pa] = rtpadd_emis_wis(h,ha,p,pa);
  %addpath /asl/rtp_prod2/emis/
  %addpath /asl/rtp_prod2/util/
  addpath /asl/packages/rtp_prod2/emis/
  addpath /asl/packages/rtp_prod2/util/    
  [p,pa] = rtp_add_emis(p,pa);

  %% now subset what we need
  oo = find(p0.findex == grans(ixx));  
  ixGeorge = (p0.atrack(oo)-1)*90 + p0.xtrack(oo);  
  ixSergio = (p.atrack-1)*90 + p.xtrack;
  [Y,iA,iB] = intersect(ixGeorge,ixSergio);
  fprintf(1,'  compare rlon/rlat/rtime George vs SergioSubset %8.6f %8.6f %8.6f\n',[sum(abs(p0.rlon(oo) - p.rlon(iB))) sum(abs(p0.rlat(oo) - p.rlat(iB))) sum(abs(p0.rtime(oo) - p.rtime(iB)))])

  %% looking EW at fixed point
  booA1 = find(p.atrack == atrack1(ixx));
  booA2 = find(p.atrack == atrack2(ixx));
  booA3 = find(p.atrack == atrack3(ixx));

  %% looking NS at fixed point
  booB1 = find(p.xtrack == xtrack1(ixx));
  booB2 = find(p.xtrack == xtrack2(ixx));
  booB3 = find(p.xtrack == xtrack3(ixx));
  
  iCA = union(union(booA1,booA2),booA3);
  iCB = union(union(booB1,booB2),booB3);
  
  iC = union(iCA,iCB);
  iC = iCA;    %% just keep looking EW
  iD = union(iB,iC);

  figure(1); clf
  figure(2); clf
  
  figure(1); plot(p0.rlon(oo),p0.rlat(oo),'b.',p.rlon(iD),p.rlat(iD),'ro'); title(num2str(length(oo)));
    ax = axis; line([ax(1) ax(2)],[rlat_atrack1(ixx) rlat_atrack1(ixx)],'color','r','linewidth',2)
    ax = axis; line([ax(1) ax(2)],[rlat_atrack2(ixx) rlat_atrack2(ixx)],'color','m')
    ax = axis; line([ax(1) ax(2)],[rlat_atrack3(ixx) rlat_atrack3(ixx)],'color','m')
    ax = axis; line([rlon_xtrack1(ixx) rlon_xtrack1(ixx)],[ax(3) ax(4)],'color','b','linewidth',2)
    ax = axis; line([rlon_xtrack2(ixx) rlon_xtrack2(ixx)],[ax(3) ax(4)],'color','c')
    ax = axis; line([rlon_xtrack3(ixx) rlon_xtrack3(ixx)],[ax(3) ax(4)],'color','c')
  figure(2); plot(p0.xtrack(oo),p0.atrack(oo),'b.',p.xtrack(iD),p.atrack(iD),'ro'); title(num2str(length(oo)));
    ax = axis; line([ax(1) ax(2)],[atrack1(ixx) atrack1(ixx)],'color','r','linewidth',2)
    ax = axis; line([ax(1) ax(2)],[atrack2(ixx) atrack2(ixx)],'color','m')
    ax = axis; line([ax(1) ax(2)],[atrack3(ixx) atrack3(ixx)],'color','m')    
    ax = axis; line([xtrack1(ixx) xtrack1(ixx)],[ax(3) ax(4)],'color','b','linewidth',2)
    ax = axis; line([xtrack2(ixx) xtrack2(ixx)],[ax(3) ax(4)],'color','c')
    ax = axis; line([xtrack3(ixx) xtrack3(ixx)],[ax(3) ax(4)],'color','c')
  %disp('ret'); pause
  pause(0.1);
  
  fprintf(1,'ix = %3i of %3i : AIRS granule = %3i ixGeorge = %5i ixSergio = %5i \n',ixx,length(grans),ix,length(iA),length(iD));
  if isfield(p,'robsqual')
    p = rmfield(p,'robsqual');
  end
  [h2,p2] = subset_rtp_allcloudfields(h,p,[],[],iB);
  h2.vchan = h0.vchan;
  
  [h,p] = subset_rtp_allcloudfields(h,p,[],[],iD);
  h.vchan = h0.vchan;
  
  if ixx == 1
    hall = h;
    pall = p;

    hjunk = h2;
    pjunk = p2;

  else
    [hjunk,pjunk] = cat_rtp(hjunk,pjunk,h2,p2);
    [hall,pall] = cat_rtp(hall,pall,h,p);    
  end
end

fprintf(1,'FINALJUNK ixGeorge = %5i ixSergio = %5i \n',length(p0.rlon),length(pjunk.rlon));
[sum(p0.rlon-pjunk.rlon) sum(p0.rlat-pjunk.rlat) sum(p0.rtime-pjunk.rtime)]

fprintf(1,'FINAL     ixGeorge = %5i ixSergio = %5i \n',length(p0.rlon),length(pall.rlon));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h = hall;
p = pall;

[p2] = driver_sarta_cloud_rtp(h,ha,p,pa,run_sarta);

fnamex = fnameOUT;
[h,ha,p2x,pa] = rtptrim_sartacloud(h,ha,p2,pa);
rtpwrite(fnamex,h,ha,p2x,pa)

