function proc_tiles_to_daily(ryear, rzone)

% Data are supplied in 16-day collections (23 collections per year)
% with 64 zonal bands each zone has 72 longitude tiles.
%
%
% One task to: load all tiles for zonal band rzone.
%              for one 16-day set (collection)
%              for given year.
% Compute day radiance average for each tile and 
%             whole zonal band, for good channels.
% Save daily averages.
% task: INPUTS: ryear, [integer] 2002 to 2020.
%               rzone. [integer] 1 to 64 (1:S.P, 64:N.P.)
%
% There are 72 tiles/zonal band. 64 zonal bands/set.
%           23 sets/year. 18 years
% Round 1: for SSW detection. Analyse NP and SP zone bands.
%
% 2020Nov18 CLH: change to one calendar year gulp. 

addpath /asl/matlib/time
addpath /asl/matlib/aslutil
addpath /home/motteler/shome/chirp_test             % read_netcdf_h5

% Check input parameters
if(~ismember(ryear,[2002:2020])); error('Invalid year'); return; end

%if(~ismember(rset,[1:23])); error('invalid set number [1:23]'); return; end

%if(ryear == 2002 & rset > 8); error('2002 has sets 1..8 only'); return; end

if(~ismember(rzone, [1:64])); error('Invalid zone band [ 1:64]'); return; end

% construct set directory string
%nset = (ryear-2002)*23 -16 + rset;
%set_str = sprintf('%4i_s%03i', ryear, nset);
%disp(['Processing set: ' set_str])

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
zn_str = BZ{rzone}
clear AZ;

% preferred wnums for strat temperature 720 to 760 cm-1 (fairs)
load('/home/chepplew/myLib/data/airs_f.mat','fairs');
achs = find(fairs > 660 & fairs < 721);

% Use only 'good' channels (chanset)
load('/home/chepplew/myLib/data/airs_chanset_for_retrievals.mat');
xchs = intersect(achs, chanset);

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

[xlen, ~] = size(d.ppar);     % check = 72* {22 or 23} = 1584 or 1656 
junk = xlen/72;
if(floor(xlen/72) ~= xlen/72)
  error('incorrect number zonal bins'); return;
end
slen = junk;
zlen = floor(xlen/slen); % 72;     % always 72

% File listings are in order: set:1.zone 1:72, set:2.zone:1-72
% zone listing order: E00p00..E175p00, W005p00..W180p00
tsize = []; nsam = [];
for zk = 1:zlen
  %for sk = 1:slen
    bin(zk) = struct('tim',[],'lat',[],'lon',[],'rad',[],'nsam',[],'day',[]);
  %end
end

disp('Starting load tile data')
xk = 0; sk = 1; zk = 0;
for ifn = 1:length(d.ppar)
    if(zk > 0 & rem(zk,72) == 0) 
      fprintf(1,'.')
      sk = sk+1;
      zk = 0;
    end
    zk = zk+1;
%    disp(['ifn: ' num2str(ifn) ' sk: ' num2str(sk) ' zk: ' num2str(zk)]);
    
    tilename = [d.ppar(ifn).folder '/' d.ppar(ifn).name];
    tilesize = [d.ppar(ifn).bytes];

    [s,~] = read_netcdf_h5(tilename);
    inan = find(s.lat > 1E6);
    s.tai93([inan]) = [];
    s.lat([inan])   = [];
    s.lon([inan])   = [];
    s.rad(:,[inan]) = [];
    dnum = tai2dnum(airs2tai(s.tai93));
  %
    bin(zk).tim   = [bin(zk).tim; dnum];
    bin(zk).lon   = [bin(zk).lon; s.lon];
    bin(zk).lat   = [bin(zk).lat; s.lat];
    bin(zk).rad   = [bin(zk).rad; s.rad(xchs,:)'];
    bin(zk).nsam  = [bin(zk).nsam; length(s.lat)];  

end

% do daily averages (required to combine samples from other tiles)
clear tS; for zk=1:zlen tS(zk) = floor(min(cell2mat({bin(zk).tim}))); end
clear tE; for zk=1:zlen tE(zk) = ceil(max(cell2mat({bin(zk).tim}))); end
tlen = tE-tS;
clear rdays; for zk=1:zlen rdays{zk} = [tS(zk):1:tE(zk)-1]; end

for zk = 1:zlen
  junk = bin(zk).tim;
  for iday = 1:tE(zk)-tS(zk)
      itt = find(tS(zk)+iday-1 < junk & junk <= tS(zk)+iday);
      bin(zk).day(iday).rmn = nanmean(bin(zk).rad(itt,:));
      bin(zk).day(iday).itot = numel(itt);
  end
end

% reformat daily radiance means
clear radmn;
for ich = 1:length(xchs)
  for zk = 1:zlen
    junk=[]; for i=1:tlen(zk) junk = [junk bin(zk).day(i).rmn(ich)]; end
    radmn{ich,zk} = junk;
  end
end


%{

junk={bin([1:72]).nsam};
plot(datetime(bin(1).tim, 'convertfrom','datenum'),'.')
plot(datetime(bin(1).tim, 'convertfrom','datenum'), bin(1).rad(:,7),'.')
hold on;
plot(datetime(rdays{1},'convertfrom','datenum'), radmn{7,1},'-')

%}




% Save daily radiance means.
g.home = '/home/chepplew/data/rates_anomalies/tiling/daily_means/';
if(~exist(g.home)) mkdir(g.home); end
savfn = [num2str(ryear) '_lw40_' zn_str '_daymean.mat'];

% select variables to save
lonB = [-180:5:180];
zonB = BZ; clear BZ;

% save
disp(['Saving to file: ' g.home savfn]);
save([g.home savfn],'lonB','zonB','rdays','radmn','fairs','xchs','-v7.3');


% END

%{
zk = 1;
rjunk = []; for i=1:tlen rjunk = [rjunk; bin(zk).day(i).rmn]; end
tjunk = days; 
figure; plot(datetime(tjunk,'convertfrom','datenum'),rjunk,'.-')

%}



