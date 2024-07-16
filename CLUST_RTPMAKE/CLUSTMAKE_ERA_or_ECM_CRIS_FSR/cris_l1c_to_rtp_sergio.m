function [hd0, hattr, pd0, pattr, tstr] = cris_l1c_to_rtp_sergio(yy,mm,dd,rgrans,NWPmodel,iSNPPorJ1orJ2,iInterp)

tstr = [];

%function [hd0, pd0] = cris_l1c_to_rtp(rdate, rgrans)
% cp -a /home/chepplew/projects/cris/cris_l1c_to_rtp.m .
% cris_l1c_to_rtp()
%
%
% Currently 1 or 2 arguments, rdate and optional granule numbers.
% Choose date to process: rdate='2018/08/01';
% Optional rgrans member [1..240] (max 4).
%    determines: clear subset or allsky scattering: prod = CLR,SCT
% Hardwired SRC{1} = NSR or FSR
% Hardwired SRC{2} = SNPP or J1 or J2
% Hardwired CCAST version    'v20d'

% based on ~/projects/airs/airs_l1c_to_rtp.m
% and /home/sbuczko1/git/rtp_prod2/cris/create_cris_ccast_hires_clear_day_rtp.m

addpath ~motteler/shome/cris/ccast/source/          % fixmyQC
addpath /asl/matlib/rtptools                        % set_attr
addpath /asl/matlib/time
addpath /asl/matlib/aslutil                         % int2bits
addpath /home/sbuczko1/git/rtp_prod2/chirp/util/uniform_clear/
addpath /home/sbuczko1/git/rtp_prod2/grib/          % fill_era
addpath /home/sbuczko1/git/rtp_prod2/emis
addpath /home/sbuczko1/git/rtp_prod2/util           % genscratchpath
addpath /home/sbuczko1/git/rtp_prod2/cris/util      % guard_ind
addpath /home/sbuczko1/git/rtp_prod2/cris/util/uniform_clear
addpath /home/chepplew/projects/chirp               % cat_rtp_clh
addpath /home/sergio/MATLABCODE

% ---------- Input Options and setup ----------------

if (nargin == 4)
  error('dont wanna do clear subset')
  % only date supplied so do clear subset
  prod = 'clr';
  prod = 'sct';
  NWPmodel = 'era';
  iSNPPorJ1orJ2 = 0;    %SNPP
  iInterp = +1;         %interp the analysis
elseif (nargin == 5)
  % date and granule number so do allsky (sct = scattering)
  prod = 'sct';
  iSNPPorJ1orJ2 = 0;    %SNPP
  iInterp = +1;         %interp the analysis
elseif (nargin == 6)
  prod = 'sct';
  iInterp = +1;         %interp the analysis
elseif (nargin == 7)
  prod = 'sct';
end
  
% Check option: subset or Scattering allsky RTP production.
if(~ismember(prod,{'clr','sct'})) 
  error('Invalid run type (clear or scattering)');
  return;
end

% Switch source depending on {NPP,J01}.CrIS,  & {NSR, FSR}.
% e.g. SRC = {'NPP','FSR'}
if iSNPPorJ1orJ2 == 0
  SRC = {'NPP','FSR'};
elseif iSNPPorJ1orJ2 == 1
  SRC = {'J01','FSR'};
elseif iSNPPorJ1orJ2 == 2
  SRC = {'J02','FSR'};
end

SRC = lower(SRC);
if(~ismember(SRC{1},{'npp','j01'})); error('Invalid CrIS'); return; end
if(~ismember(SRC{2},{'nsr','fsr'})); error('Invalid spec.res'); return; end

switch SRC{1}
  case 'npp'
    switch SRC{2}
      case 'fsr'
        d.home = '/asl/cris/ccast/sdr45_npp_HR/';
      case 'nsr'
        d.home = '/asl/cris/ccast/sdr45_npp_LR/';
    end
  case 'j01'
    switch SRC{2}
      case 'fsr'
        d.home = '/asl/cris/ccast/sdr45_j01_HR/';
      case 'nsr'
        d.home = '/asl/cris/ccast/sdr45_j01_LR/';
    end
end   

% Check valid date requested
rdate = [num2str(yy,'%04d') '/' num2str(mm,'%02d') '/' num2str(dd,'%02d')];
fprintf('rdate = %s \n', rdate);
try
   rD    = datenum(rdate,'yyyy/mm/dd');
   jday  = day(datetime(rD,'convertfrom','datenum'),'dayofyear');
   syear = year(datetime(rD,'convertfrom','datenum'));
catch
   error('Incorrect Date Format')
   return
end

% complete directory path
d.home = [d.home sprintf('/%4d/%03d/',syear,jday)];

% CCAST processing version and get file listing
daysSince2002_0 = change2days(2019,03,31,2002);
daysSince2002_x = change2days(yy,mm,dd,2002);
if daysSince2002_x > daysSince2002_0
  vers = 'v20d';  %% wot Chris H code has
else
  vers = 'v20a';  %% hmm for HALO  /asl/cris/ccast/sdr45_j01_HR//2019/115/  I only see v20a
end
fprintf(1,' daysSince2002_x = %4i so use %s \n',daysSince2002_x,vers);

if (~ismember(vers,{'v20d'})) & (~ismember(vers,{'v20a'}))
  error('invalid ccast version'); 
  return; 
end  

fpattern = ['CrIS_SDR_*' vers '.mat'];
%[d.home fpattern]
d.dir = dir([d.home fpattern]);

% Check SDRs exist
if(length(d.dir) <= 1) 
  d.home
  fpattern
  error('Insufficient SDR granules found'); 
return; end 
  
% reorder the listing into granule order (shouldn't be needed)
gnum = [];
for fn=1:length(d.dir) 
  junk  = strsplit(d.dir(fn).name,{'_','.'});
  gnum  = [gnum str2double(junk{7}(2:end))];
  gajunk = junk{6};
  thetime(fn,1) = str2num(gajunk(2:3));
  thetime(fn,2) = str2num(gajunk(4:5));
  thetime(fn,3) = thetime(fn,1)*10 + (round(thetime(fn,2)/6)+0);   %%% GONNA HAVE PROBLEMS with FIRST and LAST
end
[ia ib] = sort(gnum);

% Check requested granule numbers to process if allsky
if(strcmp(prod, 'sct'))
  iign = intersect(rgrans, gnum);
  [~,~,iign] = intersect(rgrans, thetime(:,3));
  tstr = strsplit(d.dir(iign).name,{'_','.'});
  tstr = tstr{6};
  fprintf(1,'rgrans %3i trstr = %s \n',rgrans,tstr);
  if(~all(thetime(iign,3) == rgrans))
    error('Cant find all requested granules to process')
    return; end
else if(strcmp(prod, 'clr'))
  iign = 1:length(d.dir);
  end
end
  
%  -------------- Initialize for Processing ------------------

% Number of guard channels
nguard = 2;

% seconds between 1 Jan 1958 and 1 Jan 2000
tdif = 15340 * 24 * 3600;

[sID, sTempPath] = genscratchpath();
fn_rtp1 = fullfile(sTempPath, ['cris_' sID '.rtp']);
fn_rtp2 = fullfile(sTempPath, ['cris_' sID '_kl.rtp']);
fn_rtp3 = fullfile(sTempPath, ['cris_' sID '_sar.rtp']);

fn_rtp1 = mktempS('fx.ip.rtp');
fn_rtp2 = mktempS('fx.op.rtp');
fn_rtp3 = mktempS('fx.rp.rtp');
fn_rtp4 = mktempS('fx.xp.rtp');

% assign executables and command strings (clr or sct)
klayers_bin      = '/asl/packages/klayersV205/BinV201/klayers_airs_wetwater';
sartaclr_bin.nsr = '/home/chepplew/gitLib/sarta/bin/crisg4_oct16';
sartasct_bin.nsr = '';
%sartaclr_bin.fsr = '/home/chepplew/gitLib/sarta/bin/cris_hrg4_p2019dec18';
sartaclr_bin.fsr = '/home/chepplew/gitLib/sarta/bin/crisg4_oct16_aug20';
sartasct_bin.fsr = '/home/chepplew/gitLib/sarta/bin/crisg4_hires_dec17_iceGHMbaum_wdrop_ddust_sc_hg3_new';

klayers_run = [klayers_bin ' fin=' fn_rtp1 ' fout=' fn_rtp2 ' >& ugh'];

switch SRC{2}
  case 'nsr'
     sarta_run = [sartaclr_bin.nsr ' fin=' fn_rtp2 ' fout=' fn_rtp3 ' >& ugh'];
     if(strcmp(prod,'sct'))
     % TBD
     end
  case 'fsr'
     sarta_run = [sartaclr_bin.fsr ' fin=' fn_rtp2 ' fout=' fn_rtp3 ' >& ugh'];
     if(strcmp(prod,'sct'))
     sarta_run = [sartasct_bin.fsr ' fin=' fn_rtp2 ' fout=' fn_rtp3 ' >& ugh'];
     end
end  

disp(['Using SARTA run cmd: ' sarta_run]);

% Options
cfg = struct;
cfg.model = 'era';
cfg.model = NWPmodel;

uniform_cfg = struct;
uniform_cfg.uniform_test_channel = 961;   % ??? 900, (961) or 123{1,2} cm-1
uniform_cfg.uniform_bt_threshold = 0.4;   % 0.4
uniform_cfg.cscanlines           = 45;
uniform_cfg.ascanlines           = 135;

clear_cfg.clear_test_channel = 961;
clear_cfg.clear_ocean_bt_threshold = 4;
clear_cfg.clear_land_bt_threshold = 7;

% used for concatenating RTP
isfirst = 1;

trace = struct;
trace.githash = 'na';
trace.RunDate = 'na';

%  ======= Main loop over Granules for the day ============
for fn=iign

  %disp(['fn: ' num2str(fn)])
  try
    fnamex = [d.dir(fn).folder '/' d.dir(fn).name]; fprintf(1,'file %3i %s \n',fn,fnamex);
    load([d.dir(fn).folder '/' d.dir(fn).name]);
  catch ME
    disp(ME.identifier)
    continue;
  end

  junk = strsplit(d.dir(fn).name,{'.','_'});
  gran.date = junk{5}(2:end);
  gran.num  = str2double(junk{7}(2:end));

  nedn_LW = nLW;
  nedn_MW = nMW;
  nedn_SW = nSW;

  % determine source instrument from attributes
  % <TBD>
  
  % Assign Header variables
  %h = struct;
  %h.pfields = 4;  % robs1, no calcs in file
  %h.ptype   = 0;
  %h.ngas    = 0;
  %% product_name_platform: "SS1330"
  %h.instid  = 800; % AIRS
  %h.pltfid  = -9999;

  % sanity check for ccast QC
  if exist ('L1a_err') ~= 1
    [d.dir(fn).folder '/' d.dir(fn).name]
    error('L1a_err flags missing in ccast SDR file')
  end

  % get total obs count (m = 9 ifovs, n = xtrack)
  [nfov, nfor, nscan] = size(geo.Latitude);
  nobs = nfov * nfor * nscan;

  %---------------
  % copy geo data
  %---------------
  p = struct;
  p.rlat = single(geo.Latitude(:)');
  p.rlon = single(geo.Longitude(:)');

  p.rlon = wrapTo180(p.rlon);

  %p.rtime = reshape(ones(9,1) * (geo.FORTime(:)' * 1e-6 - tdif), 1, nobs);
  p.rtime = reshape(ones(9,1) * (geo.FORTime(:)' * 1e-6 ), 1, nobs);
  p.satzen = single(geo.SatelliteZenithAngle(:)');
  p.satazi = single(geo.SatelliteAzimuthAngle(:)');
  p.solzen = single(geo.SolarZenithAngle(:)');
  p.solazi = single(geo.SolarAzimuthAngle(:)');
  % Incorrect
  %p.zobs = single(geo.Height(:)');
  % SatelliteRange is zobs for nadir
  temp = squeeze(geo.SatelliteRange(5,:,:));
  temp = (nanmean(temp(15,:),2) + nanmean(temp(16,:),2))/2;
  p.zobs = ones(1,nobs)*temp;
  clear temp;

addpath /home/sergio/MATLABCODE/TIME
[xyy,xmm,xdd,xhh] = tai2utcSergio(p.rtime);        %%% <<<<<<<<<<<<<<<<<<<<<<<<<<<<< for SdSM old time
time_so_far = (xyy-2000) + ((xmm-1)+1)/12;
co2ppm = 368 + 2.077*time_so_far;  %% 395.6933
p.co2ppm = co2ppm;
fprintf(1,'CLIMATOLOGY co2ppm for FIRST %4i/%2i/%2i = %8.6f ppmv\n',xyy(1),xmm(1),xdd(1),p.co2ppm(1));
fprintf(1,'CLIMATOLOGY co2ppm for LAST  %4i/%2i/%2i = %8.6f ppmv\n',xyy(end),xmm(end),xdd(end),p.co2ppm(end));

  iobs = 1:nobs;
  p.atrack = int32( 1 + floor((iobs-1)/(nfov*nfor)) );
  p.xtrack = int32( 1 + mod(floor((iobs-1)/9),30) );
  p.ifov = int32( 1 + mod(iobs-1,9) );
  %p.udef      = [];
  %p.iudef     = [];
  %--------------------
  % copy radiance data
  %--------------------
  sg = 2;       % number of src guard chans
  dg = nguard;  % number of dst guard chans

  % true channel set sizes
  nLW = length(vLW) - 2 * sg;
  nMW = length(vMW) - 2 * sg;
  nSW = length(vSW) - 2 * sg;

  % total number of output channels
  nout = nLW + nMW + nSW + 6 * dg;

  % initialize radiance output
  p.robs1 = ones(nout, nobs, 'single') * NaN;

  [si, di] = guard_ind(sg, dg, nLW);
  rtmp = reshape(rLW, length(vLW), [], nobs);
  p.robs1(di, :) = single(rtmp(si, :));

  [si, di] = guard_ind(sg, dg, nMW);
  di = nLW + 2 * dg + di;
  rtmp = reshape(rMW, length(vMW), nobs);
  p.robs1(di, :) = single(rtmp(si, :));

  [si, di] = guard_ind(sg, dg, nSW);
  di = nLW + nMW + 4 * dg + di;
  rtmp = reshape(rSW, length(vSW), nobs);
  p.robs1(di, :) = single(rtmp(si, :));

  % set to 1, for now
  p.robsqual = zeros(1, nobs, 'single');

  % observer pressure
  p.pobs = zeros(1,nobs,'single');

  % upwelling radiances
  p.upwell = ones(1,nobs,'int32');
 
  %--------------------
  % set the p udefs
  %--------------------
  p.udef  = zeros(20, nobs, 'single');
  p.iudef = zeros(10, nobs, 'int32');

  % iudef 3 is granule ID as an int32
  %t1 = str2double(cellstr(geo.Granule_ID(:,4:16)))';
  t2 = int32(ones(nobs,1) * gran.num);
  p.iudef(3,:) = t2(:)';

  % iudef 4 is ascending/descending flag
  t1 = geo.Asc_Desc_Flag';
  t2 = int32(ones(nfov*nfor,1) * t1);
  p.iudef(4,:) = t2(:)';

  % iudef 5 is orbit number
  t1 = geo.Orbit_Number';
  t2 = int32(ones(nfov*nfor,1) * t1);
  p.iudef(5,:) = t2(:)';

  %-------------------------------
  % trim output to a valid subset
  %-------------------------------
  % get good data index
  iok = find(reshape(ones(9,1) * ~L1a_err(:)', 1, nobs));
  [eLW, eMW, eSW] = fixmyQC(L1a_err, L1b_stat);
  etmp = eLW | eMW | eSW;
  iok = find(~etmp(:)');

  p.rlat   = p.rlat(:, iok);
  p.rlon   = p.rlon(:, iok);
  p.rtime  = p.rtime(:, iok);
  p.satzen = p.satzen(:, iok);
  p.satazi = p.satazi(:, iok);
  p.solzen = p.solzen(:, iok);
  p.solazi = p.solazi(:, iok);
  p.zobs   = p.zobs(:, iok);
  p.pobs   = p.pobs(:, iok);
  p.upwell = p.upwell(:, iok);
  p.atrack = p.atrack(:, iok);
  p.xtrack = p.xtrack(:, iok);
  p.ifov   = p.ifov(:, iok);
  p.robs1  = p.robs1(:, iok);
  p.robsqual = p.robsqual(:, iok);
  p.udef   = p.udef(:, iok);
  p.iudef  = p.iudef(:, iok);

  % Assign attribute strings
  pattr = struct;
  pattr={{'profiles' 'iudef(1,:)' 'Dust flag:[1=true,0=false,-1=land,-2=cloud,-3=bad data]'},...
      {'profiles' 'iudef(2,:)' 'Dust_score:[>380 (probable), N/A if Dust Flag < 0]'},...
      {'profiles' 'iudef(3,:)' 'SceneInhomogeneous:[128=inhomogeneous,64=homogeneous]'},...
      {'profiles' 'iudef(4,:)' 'scan_node_type [0=Ascending, 1=Descending]'},...
      {'profiles' 'udef(1,:)' 'sun_glint_distance:[km to sunglint,-9999=unknown,30000=no glint]'},...
      {'profiles' 'udef(2,:)' 'spectral_clear_indicator:[2=ocean clr,1=ocean n/clr,0=inc. data,-1=land n/clr,-2=land clr]'},...
      {'profiles' 'udef(3,:)' 'BT_diff_SO2:[<-6, likely volcanic input]'},...
      {'profiles' 'udef(4,:)' 'Inhomo850:[abs()>0.84 likely inhomogeneous'},...
      {'profiles' 'udef(5,:)' 'Rdiff_swindow'},...
      {'profiles' 'udef(6,:)' 'Rdiff_lwindow'}};
 
  %%pattr = set_attr(pattr, 'robs1', 'inpath');
  %%pattr = set_attr(pattr, 'rtime', 'TAI:1958');

  noise = [nLW(:,1,1); nMW(:,1,1); nSW(:,1,1)];

  %-------------------
  % set header values
  %-------------------
  h = struct;
  h.nchan = nout;
  h.ichan = cris_ichan(nguard, 2, nLW, nMW, nSW);
  h.vchan = cris_vchan(nguard, userLW, userMW, userSW);
  h.ichan = h.ichan';
  h.vchan = h.vchan';
  h.pfields = 4; % 4 = IR obs
  h.ptype   = 0;
  h.vcmax   = max(h.vchan);
  h.vcmin   = min(h.vchan);

 %% keyboard_nowindow
 %% plot(h.vchan,noise)

  [salti, landfrac] = usgs_deg10_dem(p.rlat, p.rlon);
  p.landfrac = landfrac;
  p.salti    = salti;

  %-----------------------
  % set header attributes
  %-----------------------
  hattr = {{'header', 'instid', 'CrIS'}, ...
           {'header', 'reader', 'ccast2rtp'}, ...
          };


  i900 = find(h.vchan >= 900,1);
  figure(1); scatter_coast(p.rlon,p.rlat,25,rad2bt(900,p.robs1(i900,:))); 
  if iSNPPorJ1orJ2 == 0
    title('BT900 CRIS NPP FSR obs'); 
  elseif iSNPPorJ1orJ2 == 1
    title('BT900 CRIS J1 FSR obs'); 
  elseif iSNPPorJ1orJ2 == 2
    title('BT900 CRIS J2 FSR obs'); 
  end
  
  colormap jet
  pause(0.1);

  %% /umbc/xfs2/strow/asl/s1/sbuczko1/git/rtp_prod2/grib/fill_ecmwf.m cn mess up p.plon (needs wrapTo180)
  addpath /home/sergio/MATLABCODE/RTPMAKE/CLUST_RTPMAKE/GRIB

  % Add in model data ******************************
  fprintf(1, '>>> Add model: %s...', cfg.model)
  switch cfg.model
   case 'ecmwf'
     if iInterp <= 0
       [p,h,pattr]  = fill_ecmwf(p,h,pattr);
     else
       new_8_ECMfiles_interp_analysis
     end

   case 'era'
     if iInterp <= 0
       [p,h,pattr]  = fill_era(p,h,pattr);
       %[p,h,pattr]  = fill_era_interp(p,h,pattr);
     else
       new_4_ERAfiles_interp_analysis
     end

   case 'merra'
     [p,h,pattr]  = fill_merra(p,h,pattr);
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % add landfrac then emissivity
  %[~,~,p,pattr] = rtpadd_usgs_10dem(h,hattr,p,pattr);
  %[p,pattr]     = rtp_add_emis(p,pattr);

  %[h,ha,p,pa] = rtpadd_emis_DanZhou2(h,ha,p,pa);
  %p = Prof_add_emis(p,yymmddgg(1),yymmddgg(2),yymmddgg(3));  %% broken crap by whoever
  %p = rtpadd_emis_DanZhou(h,ha,p,pa);   %% lso totally broken crap
  %[h,ha,p,pa] = rtpadd_emis_wis(h,ha,p,pa);
  %addpath /asl/rtp_prod2/emis/
  %addpath /asl/rtp_prod2/util/    
  %addpath /asl/packages/rtp_prod2/emis/
  %addpath /asl/packages/rtp_prod2/util/
  %addpath /asl/rtp_prod2/emis/
  %addpath /asl/rtp_prod2/util/
 
  addpath /home/sergio/MATLABCODE/matlib/rtp_prod2/emis/
  addpath /home/sergio/MATLABCODE/matlib/rtp_prod2/util/

  p.rlon = wrapTo180(p.rlon);
  [p,pattr] = rtp_add_emis(p,pattr);

  h.pfields = 5;  % robs, model

  % Save file and copy raw rtp data.
  %rtpwrite(fn_rtp1,h,hattr,p,pattr)

  if(strcmp(prod,'clr'))
    % ================ clear Subset ===================
    % 1. Get uniform scene subset
    disp('clear subset calcs')
  
    [iuni, amax_keep] = cris_find_uniform(h, p, uniform_cfg);
    %[dbtun, mbt, amax] = chirp_find_uniform(h, p, uniform_cfg);
  
    %iuniform = find(abs(dbtun) < 1.0);   % my version
    %iuniform = find(abs(mbt) < 1.0);      % slb version
    nuni = length(iuni);
    if 0 == nuni
      fprintf(2,['>> No uniform FOVs found for granule %d. ' ...
                'SKIPPING\n'],fn)
      continue;
    end
  
    fprintf(1, '>> Uniform obs found: %d/12150\n', nuni);
    pdu = rtp_sub_prof(p,iuni);
  
    rtpwrite(fn_rtp1,h,hattr,pdu,pattr)
  
    % Now run klayers
    unix(klayers_run);
  
    %[hd2,~,pd2,~] = rtpread(fn_rtp2);
    % Now run sarta
    unix(sarta_run);
  
    % Read in new rcalcs and insert into origin pdu field
    %stFileInfo = dir(fn_rtp3);
    [hd3,~,pd3,~] = rtpread(fn_rtp3);
    pdu.rclr  = pd3.rcalc;
    %p = rmfield(p,'rcalc');
    %clear p3;
    hd3.pfields = 7;        % includes model, obs & calc
    hd3.ptype   = 0;        % back to model levels (pre-klayers)
  
    % save a copy for comparison later
    %sav_dir = '/home/chepplew/data/rtp/chirp_AQ/clear/2018/';
    %sav_fn  = ['era_chirp_' gran.date '_' gran.num '_clear.rtp'];
    %rtpwrite([sav_dir sav_fn],hd3,hattr,pdu,pattr)
  
    %nuniform = length(pdu.rtime);
    %[iflagsc, bto1232, btc1232] = xfind_clear_loANDhires(h, p, 1:nobs);
    [iflagsc, bto1232, btc1232] = xfind_clear_hires(hd3, pdu, 1:nuni);
    %[iflagsc, bto1232, btc1232] = chirp_find_clear(hd3, pdu, clear_cfg);    % 
    %[iflagsc, btoc, btcc]       = chirp_find_clear(hd3, pdu, clear_cfg); % slb version
    
  
    iclear_all = find(iflagsc == 0);
    iclear_sea = find(iflagsc == 0 & pdu.landfrac == 0);
    iclear = iclear_sea;
    nclear = length(iclear);
    fprintf(1, '>>>> Total of %d clear & uniform obs passed test\n', nclear);
    if 0 == nclear
      fprintf(2,['>> No clear FOVs found for granule %d. ' ...
                'SKIPPING\n'],fn)
      continue;
    end
  
    p_clr = rtp_sub_prof(pdu, iclear);
  
    if isfirst
      pd0  = p_clr;
      hd0  = hd3;
    else
      % concatenate new random rtp data into running random rtp structure
      [hd0, pd0] = cat_rtp_clh(hd0, pd0, hd3, p_clr);
    end
  
  end      % end if(prod == 'clr')
  % ===================== allsky ==================
  if(strcmp(prod,'sct'))
    disp('allsky calcs')

   run_sarta.co2ppm = p.co2ppm;
   run_sarta.clear = +1;
   run_sarta.cloud = +1;
   run_sarta.cumsum = -1;    %% this is "closer" to MRO but since cliuds are at centroid, does not do too well with DCC
   run_sarta.cumsum = 9999;  %% larrabee likes this, puts clouds high so does well for DCC
   
   iSlabCld_CumSumStrowORGeorge = -1;
   if iSlabCld_CumSumStrowORGeorge > 0
     run_sarta.cumsum = 9999;  %% strow pick, cloud at PEAK of wgt fcn
   else
     run_sarta.cumsum = -1;  %% aumann pick, cloud at wgt mean of profile
   end
   run_sarta.sartaclear_code = sartaclr_bin.fsr;
   run_sarta.sartacloud_code = sartasct_bin.fsr;
   
   run_sarta
   
    if ~isfield(p,'scanang')
      p.scanang = saconv(p.satzen,p.zobs);
    end

    [p2] = driver_sarta_cloud_rtp(h,hattr,p,pattr,run_sarta);
    [h,hattr,p2x,pattr] = rtptrim_sartacloud(h,hattr,p2,pattr);
    hd0 = h;
    pd0 = p2x;
  end      % end if(prod == 'sct')
  % =====================================================

  isfirst = 0;
  hdfml('closeall')

  delete(fn_rtp1, fn_rtp2, fn_rtp3);

end

