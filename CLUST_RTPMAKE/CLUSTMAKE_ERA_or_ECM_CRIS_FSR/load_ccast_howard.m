%disp(['fn: ' num2str(fn)])
try
  fnamex = [d.dir(fn).folder '/' d.dir(fn).name]; fprintf(1,'file %3i %s \n',fn,fnamex);
  load([d.dir(fn).folder '/' d.dir(fn).name]);
catch ME
  disp(ME.identifier)
end
junk = strsplit(d.dir(fn).name,{'.','_'});
gran.date = junk{5}(2:end);
gran.num  = str2double(junk{7}(2:end));

nedn_LW = nLW;
nedn_MW = nMW;
nedn_SW = nSW;

% determine source instrument from attributes
% <TBD>

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

iobs = 1:nobs;
p.atrack = int32( 1 + floor((iobs-1)/(nfov*nfor)) );
p.xtrack = int32( 1 + mod(floor((iobs-1)/9),30) );
p.ifov = int32( 1 + mod(iobs-1,9) );

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
