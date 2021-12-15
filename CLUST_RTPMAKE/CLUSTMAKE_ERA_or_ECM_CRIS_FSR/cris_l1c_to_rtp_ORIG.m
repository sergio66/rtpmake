function [hd0, pd0] = cris_l1c_to_rtp(rdate, rgrans)
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

% ---------- Input Options and setup ----------------

if(nargin == 1)
  % only date supplied so do clear subset
  prod = 'clr';
else if(nargin == 2)
  % date and granule number so do allsky (sct = scattering)
  prod = 'sct';
  if(length(rgrans) > 4) error('too many granules for allsky');
    return; end
  end
end
  
% Check option: subset or Scattering allsky RTP production.
if(~ismember(prod,{'clr','sct'})) 
  error('Invalid run type (clear or scattering)');
  return;
end

% Switch source depending on {NPP,J01}.CrIS,  & {NSR, FSR}.
% e.g. SRC = {'NPP','FSR'}
SRC = {'NPP','FSR'}
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
whos rdate; disp(rdate); fprintf('\n');
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
vers = 'v20d';
if(~ismember(vers,{'v20d'})) error('invalid ccast version'); return; end  

fpattern = ['CrIS_SDR_*' vers '.mat'];
d.dir = dir([d.home fpattern]);

% Check SDRs exist
if(length(d.dir) <= 1) error('Insufficnet SDR granules found'); return; end 

% reorder the listing into granule order (shouldn't be needed)
gnum = [];
for fn=1:length(d.dir) 
  junk  = strsplit(d.dir(fn).name,{'_','.'});
  gnum  = [gnum str2double(junk{7}(2:end))];
end
[ia ib] = sort(gnum);

% Check requested granule numbers to process if allsky
if(strcmp(prod, 'sct'))
  iign = intersect(rgrans, gnum);
  if(~all(iign == rgrans))
    error('Can.t find all requested granules to process')
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

% assign executables and command strings (clr or sct)
klayers_bin      = '/asl/packages/klayersV205/BinV201/klayers_airs_wetwater';
sartaclr_bin.nsr = '/home/chepplew/gitLib/sarta/bin/crisg4_oct16';
sartasct_bin.nsr = '';
%sartaclr_bin.fsr = '/home/chepplew/gitLib/sarta/bin/cris_hrg4_p2019dec18';
sartaclr_bin.fsr = '/home/chepplew/gitLib/sarta/bin/crisg4_oct16_aug20';
sartasct_bin.fsr = '/home/chepplew/gitLib/sarta/bin/crisg4_hires_dec17_iceGHMbaum_wdrop_ddust_sc_hg3_new';

klayers_run = [klayers_bin ' fin=' fn_rtp1 ' fout=' fn_rtp2 ' > ' ...
               '~/logs/klayers/klout.txt'];
switch SRC{2}
  case 'nsr'
     sarta_run = [sartaclr_bin.nsr ' fin=' fn_rtp2 ' fout=' fn_rtp3 ...
       ' > ~/logs/sarta/sarta_out.txt'];
     if(strcmp(prod,'sct'))
     % TBD
     end
  case 'fsr'
     sarta_run = [sartaclr_bin.fsr ' fin=' fn_rtp2 ' fout=' fn_rtp3 ...
       ' > ~/logs/sarta/sarta_out.txt'];
     if(strcmp(prod,'sct'))
     sarta_run = [sartasct_bin.fsr ' fin=' fn_rtp2 ' fout=' fn_rtp3 ...
       ' > ~/logs/sarta/sarta_out.txt'];
     end
end  

disp(['Using SARTA run cmd: ' sarta_run]);

% Options
cfg = struct;
cfg.model = 'era';

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

  disp(['fn: ' num2str(fn)])
  try
    load([d.dir(fn).folder '/' d.dir(fn).name]);
  catch ME
    disp(ME.identifier)
    continue;
  end

  junk = strsplit(d.dir(fn).name,{'.','_'});
  gran.date = junk{5}(2:end);
  gran.num  = str2double(junk{7}(2:end));

  % determine source instrument from attributes
  % <TBD>
  
  % Assign Header variables
  head = struct;
  head.pfields = 4;  % robs1, no calcs in file
  head.ptype   = 0;
  head.ngas    = 0;
  % product_name_platform: "SS1330"
  head.instid  = 800; % AIRS
  head.pltfid  = -9999;

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
  prof = struct;
  prof.rlat = single(geo.Latitude(:)');
  prof.rlon = single(geo.Longitude(:)');
  %prof.rtime = reshape(ones(9,1) * (geo.FORTime(:)' * 1e-6 - tdif), 1, nobs);
  prof.rtime = reshape(ones(9,1) * (geo.FORTime(:)' * 1e-6 ), 1, nobs);
  prof.satzen = single(geo.SatelliteZenithAngle(:)');
  prof.satazi = single(geo.SatelliteAzimuthAngle(:)');
  prof.solzen = single(geo.SolarZenithAngle(:)');
  prof.solazi = single(geo.SolarAzimuthAngle(:)');
  % Incorrect
  %prof.zobs = single(geo.Height(:)');
  % SatelliteRange is zobs for nadir
  temp = squeeze(geo.SatelliteRange(5,:,:));
  temp = (nanmean(temp(15,:),2) + nanmean(temp(16,:),2))/2;
  prof.zobs = ones(1,nobs)*temp;
  clear temp;

  iobs = 1:nobs;
  prof.atrack = int32( 1 + floor((iobs-1)/(nfov*nfor)) );
  prof.xtrack = int32( 1 + mod(floor((iobs-1)/9),30) );
  prof.ifov = int32( 1 + mod(iobs-1,9) );
  %prof.udef      = [];
  %prof.iudef     = [];
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
  prof.robs1 = ones(nout, nobs, 'single') * NaN;

  [si, di] = guard_ind(sg, dg, nLW);
  rtmp = reshape(rLW, length(vLW), [], nobs);
  prof.robs1(di, :) = single(rtmp(si, :));

  [si, di] = guard_ind(sg, dg, nMW);
  di = nLW + 2 * dg + di;
  rtmp = reshape(rMW, length(vMW), nobs);
  prof.robs1(di, :) = single(rtmp(si, :));

  [si, di] = guard_ind(sg, dg, nSW);
  di = nLW + nMW + 4 * dg + di;
  rtmp = reshape(rSW, length(vSW), nobs);
  prof.robs1(di, :) = single(rtmp(si, :));

  % set to 1, for now
  prof.robsqual = zeros(1, nobs, 'single');

  % observer pressure
  prof.pobs = zeros(1,nobs,'single');

  % upwelling radiances
  prof.upwell = ones(1,nobs,'int32');
 
  %--------------------
  % set the prof udefs
  %--------------------
  prof.udef  = zeros(20, nobs, 'single');
  prof.iudef = zeros(10, nobs, 'int32');

  % iudef 3 is granule ID as an int32
  %t1 = str2double(cellstr(geo.Granule_ID(:,4:16)))';
  t2 = int32(ones(nobs,1) * gran.num);
  prof.iudef(3,:) = t2(:)';

  % iudef 4 is ascending/descending flag
  t1 = geo.Asc_Desc_Flag';
  t2 = int32(ones(nfov*nfor,1) * t1);
  prof.iudef(4,:) = t2(:)';

  % iudef 5 is orbit number
  t1 = geo.Orbit_Number';
  t2 = int32(ones(nfov*nfor,1) * t1);
  prof.iudef(5,:) = t2(:)';

  %-------------------------------
  % trim output to a valid subset
  %-------------------------------
  % get good data index
  iok = find(reshape(ones(9,1) * ~L1a_err(:)', 1, nobs));
  [eLW, eMW, eSW] = fixmyQC(L1a_err, L1b_stat);
  etmp = eLW | eMW | eSW;
  iok = find(~etmp(:)');

  prof.rlat   = prof.rlat(:, iok);
  prof.rlon   = prof.rlon(:, iok);
  prof.rtime  = prof.rtime(:, iok);
  prof.satzen = prof.satzen(:, iok);
  prof.satazi = prof.satazi(:, iok);
  prof.solzen = prof.solzen(:, iok);
  prof.solazi = prof.solazi(:, iok);
  prof.zobs   = prof.zobs(:, iok);
  prof.pobs   = prof.pobs(:, iok);
  prof.upwell = prof.upwell(:, iok);
  prof.atrack = prof.atrack(:, iok);
  prof.xtrack = prof.xtrack(:, iok);
  prof.ifov   = prof.ifov(:, iok);
  prof.robs1  = prof.robs1(:, iok);
  prof.robsqual = prof.robsqual(:, iok);
  prof.udef   = prof.udef(:, iok);
  prof.iudef  = prof.iudef(:, iok);

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

  %-------------------
  % set header values
  %-------------------
  head = struct;
  head.nchan = nout;
  head.ichan = cris_ichan(nguard, 2, nLW, nMW, nSW);
  head.vchan = cris_vchan(nguard, userLW, userMW, userSW);
  head.ichan = head.ichan';
  head.vchan = head.vchan';
  head.pfields = 4; % 4 = IR obs
  head.vcmax   = max(head.vchan);
  head.vcmin   = min(head.vchan);

  %-----------------------
  % set header attributes
  %-----------------------
  hattr = {{'header', 'instid', 'CrIS'}, ...
           {'header', 'reader', 'ccast2rtp'}, ...
          };

  % Add in model data ******************************
  fprintf(1, '>>> Add model: %s...', cfg.model)
  switch cfg.model
   case 'ecmwf'
    [prof,head,pattr]  = fill_ecmwf(prof,head,pattr);
   case 'era'
    [prof,head,pattr]  = fill_era(prof,head,pattr);
   case 'merra'
    [prof,head,pattr]  = fill_merra(prof,head,pattr);
  end

% add landfrac then emissivity
  [~,~,prof,pattr] = rtpadd_usgs_10dem(head,hattr,prof,pattr);
  [prof,pattr]     = rtp_add_emis(prof,pattr);

  head.pfields = 5;  % robs, model

  % Save file and copy raw rtp data.
  %rtpwrite(fn_rtp1,head,hattr,prof,pattr)


  if(strcmp(prod,'clr'))
  % ================ clear Subset ===================
  % 1. Get uniform scene subset
  disp('clear subset calcs')

  [iuni, amax_keep] = cris_find_uniform(head, prof, uniform_cfg);
  %[dbtun, mbt, amax] = chirp_find_uniform(head, prof, uniform_cfg);

  %iuniform = find(abs(dbtun) < 1.0);   % my version
  %iuniform = find(abs(mbt) < 1.0);      % slb version
  nuni = length(iuni);
  if 0 == nuni
    fprintf(2,['>> No uniform FOVs found for granule %d. ' ...
              'SKIPPING\n'],fn)
    continue;
  end

  fprintf(1, '>> Uniform obs found: %d/12150\n', nuni);
  pdu = rtp_sub_prof(prof,iuni);

  rtpwrite(fn_rtp1,head,hattr,pdu,pattr)

  % Now run klayers
  unix(klayers_run);

  %[hd2,~,pd2,~] = rtpread(fn_rtp2);
  % Now run sarta
  unix(sarta_run);

  % Read in new rcalcs and insert into origin pdu field
  %stFileInfo = dir(fn_rtp3);
  [hd3,~,pd3,~] = rtpread(fn_rtp3);
  pdu.rclr  = pd3.rcalc;
  %prof = rmfield(prof,'rcalc');
  %clear p3;
  hd3.pfields = 7;        % includes model, obs & calc
  hd3.ptype   = 0;        % back to model levels (pre-klayers)


  % save a copy for comparison later
  %sav_dir = '/home/chepplew/data/rtp/chirp_AQ/clear/2018/';
  %sav_fn  = ['era_chirp_' gran.date '_' gran.num '_clear.rtp'];
  %rtpwrite([sav_dir sav_fn],hd3,hattr,pdu,pattr)

  %nuniform = length(pdu.rtime);
  %[iflagsc, bto1232, btc1232] = xfind_clear_loANDhires(head, prof, 1:nobs);
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

  prof_clr = rtp_sub_prof(pdu, iclear);

  if isfirst
    pd0  = prof_clr;
    hd0  = hd3;
  else
    % concatenate new random rtp data into running random rtp structure
    [hd0, pd0] = cat_rtp_clh(hd0, pd0, hd3, prof_clr);
  end

  end      % end if(prod == 'clr')
  % ===================== allsky ==================
  if(strcmp(prod,'sct'))
  disp('allsky calcs')
  
  rtpwrite(fn_rtp1,head,hattr,prof,pattr)

  % Now run klayers
  unix(klayers_run);

  %[hd2,~,pd2,~] = rtpread(fn_rtp2);
  % Now run sarta
  unix(sarta_run);

  [hd3,~,pd3,~] = rtpread(fn_rtp3);
  pd3.rcld  = pd3.rcalc;
  pd3 = rmfield(pd3,'rcalc');
  hd3.pfields = 7;        % includes model, obs & calc
  hd3.ptype   = 0;        % back to model levels (pre-klayers)

  if isfirst
    pd0  = pd3;
    hd0  = hd3;
  else
    % concatenate new random rtp data into running random rtp structure
    [hd0, pd0] = cat_rtp_clh(hd0, pd0, hd3, pd3);
  end

  end      % end if(prod == 'sct')
  % =====================================================

  isfirst = 0;
  hdfml('closeall')

  delete(fn_rtp1, fn_rtp2, fn_rtp3);

end



%{
hdfml('listinfo')
hdfml('closeall')

[ix,iy] = seq_match(fcris_mr,hdc.vchan);
clf;plot(hdc.vchan, rad2bt(hdc.vchan, nanmean(pdc.robs1,2)),'.-')
hold on; plot(fcris_mr(iy), rad2bt(fcris_mr(ix), nanmean(pdc.rclr(iy,:),2)),'.-')

btc = rad2bt(fcris_mr(ix), nanmean(pdc.rclr(iy,:),2));
bto = rad2bt(hdc.vchan, nanmean(pdc.robs1,2));
btbias=bto(1:end-4) - btc(5:end);

%% 1679 channels are used but SARTA grid is not same as CHIRP L1C data.
[ix iy] = seq_match(hd3.vchan, head.vchan);
plot(head.vchan(iy), rad2bt(head.vchan(iy), nanmean(prof.robs1(iy,:),2)),'-')
   hold on;
plot(hd3.vchan(ix), rad2bt(hd3.vchan(ix), nanmean(pd3.rcalc(ix,:),2)),'-')

% Match CHIRP wnum with Sarta build:
xc=load('/home/chepplew/myLib/data/chirp_nu_official.mat');
xs=importdata('/home/chepplew/gitLib/ftc_dev/chanLists/chirp_1702_list_all');
sarta_wn = xs(:,2);
chirp_wn = xc.wnum;
[ix iy] = seq_match(chirp_wn, sort(sarta_wn));
head.ichan = iy;

%}
