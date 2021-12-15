function [rtv] = read_sergio_retrieve_tiles()

% read_sergio_retrieve_tiles.m
%
% INPUTS: selected oem fields {'skt','robs','rcal','gas_{1..9}'}
%    poemNew.{ stemp, robs1, rcalc. ...}
% OEM_initialiazation are what was used to do the first set of jacobian calcs
% _prof0 are in the rtp that was read in originally, typically this would be 
% the same as OEMinitialization, excpet for the cloud fields
%
% OUTPUT: rtv: structure with fields:
%      lat, lon, tim
%      select oem fields
%
%
% load Sergio's retrieved rates from tiled data
% filenames like: retrievalV2_cld_timestep_lonbin_01_latbin_04_JOB_021716_iDET_4_iStemp_ColWV_5_iCenterFov_-1_iCO2_Yes_No_Switch_-1.mat

%        chistats: [1x1 struct]
%          cind10: [572 794 960 1325 1520 2588 2600 274 1156 1861 2168]
%         cov_set: [2 1 1 0.2000 0.1400 0.1400 0.2000 0.1400 0.1400 0.2000 100000 5.1020e+06 5.1020e+06]
%               g: [1536x1 double]
%              ha: {{1x3 cell}  {1x3 cell}}
%         hoemNew: [1x1 struct]
%     iBadDCCOnly: 0
%      iCenterFov: -1
%            iDET: 4
%       iERAorECM: 16
%      iSondeList: -1
%    iStemp_ColWV: 5
%          iaBest: [1x412 double]
%       oem_param: [1x1 struct]
%              pa: {}
%         poemNew: [1x1 struct]
%        settings: [1x1 struct]
%           tempx: [1x1 struct]
%          themem: [1x42 double]

% [hoemNew,poemNew] are the basic rtp structure
% oem_param gives OEM settings,  and the channels used for retrieval.

%cd /home/chepplew/projects/rates_anomalies/tiled/

addpath /asl/matlib/time

% check version number in the following path strings:
d.home = ['/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/MakeAvgCldProfs2002_2020/' ...
         'LookAtTimeSeries/RTP_PROFSV2/Cld/CRODGERS_FAST_CLOUD_Retrievals/V4/'];

d.dir  = dir([d.home 'retrievalV4_cld_timestep*.mat']);

% files are loaded in latbin order first, then lonbin (ie 64 latbins for each lonbin)


% get lat and lon bin numbers of files
latN = [];  lonN = [];
for i = 1:length(d.dir)

  junk = strsplit(d.dir(i).name, {'_','.'});
  latN(i) = str2num(junk{7});
  lonN(i) = str2num(junk{5});

end
% check for missing tiles
all.lats = [1:64];
all.lons = [1:72];
ims = [];
for i = 1:64
  [ia ib ] = find(latN == i);
  if(length(ia) ~= 72)
    ji = find(diff(ib) > 64);
    disp(['WARN: missing tile at lon = ', num2str(i)]);
    ims = [ims i];
  end
end

% Get one RTP structure field names
% <TBD>

% set output variables.
lat   = [];
lon   = [];
tim   = [];
stemp = [];

% For retrieved stemp use: poemNew.stemp

for i=1:length(d.dir)
  load([d.dir(i).folder '/' d.dir(i).name]);

  if(i == 1) allnames = fieldnames(poemNew)'; 
    % Get start and end dates
    mtime = tai2dnum(poemNew.rtime);
    sdate = datetime(mtime(1),'convertFrom','datenum');
    edate = datetime(mtime(end),'convertFrom','datenum');
  end

  % Get geolocation
  lat(latN(i),lonN(i),:)   = poemNew.rlat;
  lon(latN(i),lonN(i),:)   = poemNew.rlon;
  tim(latN(i),lonN(i),:)   = tai2dnum(poemNew.rtime);

  % Get requested data fields
  stemp(latN(i),lonN(i),:) = poemNew.stemp;

  robs(latN(i),lonN(i),:,:)  = poemNew.robs1;
  rcal(latN(i),lonN(i),:,:)  = poemNew.rcalc;
   
  if(~mod(i,100)) fprintf(1,'.'); end
end

rtv.lat   = lat;
rtv.lon   = lon;
rtv.tim   = tim;
rtv.stemp = stemp;

% END
%{
% Pick selected tile: lonbin 27, latbin 47 (lat:+38.5  lon: -50.0 )
% nanmean(rtv.lat(47,27,:)) =  38.6824
% nanmean(rtv.lon(47,27,:)) = -48.4830
% ClearBin    : indX lon = 54  87.5009     indY lat = 24 -23.3750')
% nanmean(rtv.lat(24,54,:)) = -22.36.
% nanmean(rtv.lon(24,54,:)) = 87.24
% CloudBin    : indX lon = 67 152.5007     indY lat = 35   6.8743')
