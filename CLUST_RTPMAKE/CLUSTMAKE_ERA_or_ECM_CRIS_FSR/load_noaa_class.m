addpath /home/sergio/MATLABCODE/TIME

%disp(['fn: ' num2str(fn)])
try
  fnamex = [d.dir(fn).folder '/' d.dir(fn).name]; fprintf(1,'file %3i %s \n',fn,fnamex);
  x = read_netcdf_lls([d.dir(fn).folder '/' d.dir(fn).name]);
catch ME
  disp(ME.identifier)
end

junk = strsplit(d.dir(fn).name,{'.','_'});
gran.date = junk{4}(1:8);
gran.num  = str2double(junk{6}(2:end));

nedn_LW = x.nedn_lw;
nedn_MW = x.nedn_mw;
nedn_SW = x.nedn_sw;

% sanity check for ccast QC
%if exist ('L1a_err') ~= 1
%  [d.dir(fn).folder '/' d.dir(fn).name]
%  error('L1a_err flags missing in ccast SDR file')
%end

% get total obs count (m = 9 ifovs, n = xtrack)
[nfov, nfor, nscan] = size(x.lat);
nobs = nfov * nfor * nscan;

%---------------
% copy geo data
%---------------
p = struct;
p.rlat = single(x.lat(:)');
p.rlon = single(x.lon(:)');

p.rlon = wrapTo180(p.rlon);

%p.rtime = reshape(ones(9,1) * (geo.FORTime(:)' * 1e-6 - tdif), 1, nobs);
p.rtime = reshape(ones(9,1) * (offset1958_to_1993 + x.obs_time_tai93(:)'), 1, nobs);
p.satzen = single(x.sat_zen(:)');
p.satazi = single(x.sat_azi(:)');
p.solzen = single(x.sol_zen(:)');
p.solazi = single(x.sol_azi(:)');

% Incorrect
p.zobs = nanmean(single(x.sat_alt(:)')) * ones(size(p.rlon));

iobs = 1:nobs;
p.atrack = int32( 1 + floor((iobs-1)/(nfov*nfor)) );
p.xtrack = int32( 1 + mod(floor((iobs-1)/9),30) );
p.ifov = int32( 1 + mod(iobs-1,9) );

sg = 2;       % number of src guard chans
dg = nguard;  % number of dst guard chans

% true channel set sizes
vLW = x.wnum_lw;
vMW = x.wnum_mw;
vSW = x.wnum_sw;

rLW = x.rad_lw;
rMW = x.rad_mw;
rSW = x.rad_sw;

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
t1 = x.asc_flag';
t2 = int32(ones(nfov*nfor,1) * t1);
p.iudef(4,:) = t2(:)';

% iudef 5 is orbit number
%t1 = geo.Orbit_Number';
%t2 = int32(ones(nfov*nfor,1) * t1);
%p.iudef(5,:) = t2(:)';

%-------------------------------
% trim output to a valid subset
%-------------------------------
% get good data index
%iok = find(reshape(ones(9,1) * ~L1a_err(:)', 1, nobs));
%[eLW, eMW, eSW] = fixmyQC(L1a_err, L1b_stat);
%etmp = eLW | eMW | eSW;
%iok = find(~etmp(:)');
iok = 1:length(p.rlat);

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

if ~exist('userLWMWSW.mat')
  old = load('/asl/cris/ccast/sdr45_j01_HR/2023/001/CrIS_SDR_j01_s45_d20230101_t2312080_g233_v20d.mat','userLW','userMW','userSW');
  userLW = old.userLW;
  userMW = old.userMW;
  userSW = old.userSW;
  clear old
  save userLWMWSW.mat userLW userMW userSW
else
  load userLWMWSW.mat
end
