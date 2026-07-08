function proc_daily_to_emd(ryear, rzone)

%
% NB there are a total of 415 16-day collections from 2002_s001 to 2020_s416
%
%

addpath /asl/matlib/time
addpath /asl/matlib/aslutil
addpath /home/motteler/shome/chirp_test             % read_netcdf_h5
addpath /home/chepplew/projects/rates_anomalies/tiled % zone_lu.lis


% Check input parameters
if(~ismember(ryear,[2002:2020])); error('Invalid year'); return; end

if(~ismember(rzone, [1:64])); error('Invalid zone band [ 1:64]'); return; end

% preferred wnums for strat temperature 720 to 760 cm-1 (fairs)
load('/home/chepplew/myLib/data/airs_f.mat','fairs');

% Get PDF profiler
load('/home/chepplew/myLib/data/p.mat')
p = p(1:12:end);

% Construct zone string 1..32:S, 33..64:N
AZ = importdata('zone_lu.lis');
AZ = strtrim(AZ);
junk = sort(AZ);
for i=64:-1:33
  j = 65-i;
  BZ{j} = junk{i};
end
for i=1:32
  j = i+32;
  BZ{j} = junk{i};
end
zn_str = BZ{rzone};
clear AZ;

% reference longitude boundaries of tiles:
lonB = [-180:5:180];

% Get data
d.home = '/asl/isilon/airs/tile_test7/';
%d.dpar = [d.home set_str '/'];
d.dpar = dir([d.home sprintf('%4d_s*',ryear)])

% Get requested zone band tile names (e.g. YYYY_sNNN/N81p75/*.nc)
d.ppar = [];
for i = 1:length(d.dpar)
  d.ppar = [d.ppar; dir([d.dpar(i).folder '/' d.dpar(i).name '/' zn_str '/*.nc'])];
end
d.ppar

% Extract 16-day collection numbers from directory listing, used for output array index
clear ncollect;
for k = 1:length(d.dpar)
  junk = strsplit(d.dpar(k).name,'_s');
  ncollect(k) =  str2num(junk{2});
end
ncollect = unique(ncollect);

% Check numbers of tiles
[xlen, ~] = size(d.ppar);     % check = 72* {22 or 23}.
junk = xlen/72;
if(floor(xlen/72) ~= xlen/72)
  error('incorrect number zonal bins'); return;
end
slen = junk;
zlen = floor(xlen/slen); % 72;     % always 72

% choose TWP tile: tile_2003_s009_N00p00_E110p00.nc ifn = 23.
xk = 0; sk = 1; zk = 0;

nsam   = zeros(slen,zlen);     % slen or pad to 415
dnum   = zeros(slen,zlen);
pctr   = zeros(slen,zlen,2645,100,'single');
rahot  = cell(slen,zlen,2645);
bthot  = cell(slen,zlen,2645);
for ifn=1:length(d.ppar)

    if(zk > 0 & rem(zk,72) == 0) 
      fprintf(1,'.')
      sk = sk+1;
      zk = 0;
    end
    zk = zk+1;
    disp(['ifn: ' num2str(ifn) ' sk: ' num2str(sk) ' zk: ' num2str(zk) ' ' d.ppar(ifn).name]);

    tilename = [d.ppar(ifn).folder '/' d.ppar(ifn).name];
    tilesize = [d.ppar(ifn).bytes];

    [s,~] = read_netcdf_h5(tilename);
    inan = find(s.lat > 1E6);
    s.tai93([inan])    = [];
    s.lat([inan])      = [];
    s.lon([inan])      = [];
    s.rad(:,[inan])    = [];
    s.sol_zen([inan])  = [];
    s.asc_flag([inan]) = [];
    xdnum = tai2dnum(airs2tai(s.tai93));
%    all_dnum(ncollect(sk),zk) = nanmean(dnum);
    dnum(sk,zk) = nanmean(xdnum);
    ns  = size(s.rad, 2);
%    nsam(ncollect(sk),zk) = ns; 
    nsam(sk,zk) = ns; 
    pns = floor(ns/100);
    %
    bta = real(rad2bt(fairs, s.rad));        % can get -ve radiance

%   Get emd

[imf,resid,info] = emd(Ysp(ich,:),'Interpolation','pchip','MaxNumIMF',10);


end      % end for ifn loop 




%{
figure; plot(fairs, rad2bt(fairs, pctr(:,100)),'-');hold on;
   plot(fairs, rad2bt(fairs, nanmean(rahot,2)),'-')
   plot(fairs, nanmean(bthot,2),'-')

% PCA compression
xorig= squeeze(pctr(1,1,:,:));
[eigenv, scores] = pca(xorig);
mu = mean(xorig);
nComp=25;         % trial and error
Xhat = scores(:,1:nComp) * eigenv(:,1:nComp)';
Xhat = bsxfun(@plus, Xhat, mu); 

 clf;plot(fairs, rad2bt(fairs, squeeze(pctr(1,1,:,:))),'b-')
  hold on; plot(fairs, real(rad2bt(fairs, Xhat)),'g-')


%}
