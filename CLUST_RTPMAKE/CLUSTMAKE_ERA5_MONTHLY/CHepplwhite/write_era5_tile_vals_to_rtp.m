function [iok] =  write_era5_tile_vals_to_rtp(era, scen)

%
%
%
%
%

addpath /asl/matlib/h4tools

iok = 0;
% check expected fields are present

xx = load('/home/chepplew/myLib/data/airs_f.mat');

% Generate head structure
head = struct;
head.pfields = 1;  % no calcs, no obs, model in file
head.ptype   = 0;  % 1: layers, 0: levels
head.ngas    = 8;
head.gunit   = [21,21]'; %[21,21,21,21,21,21,21,21]';
head.glist   = [1, 3]';  % [1,2,3,4,5,6,9,12]';
head.vcmin   = -9999;
head.vcmax   = -9999;
head.pmin    = 0.100;
head.pmax    = 1.0371e+3;
head.nchan   = 2645;
head.vchan   = xx.fairs;
head.ichan   = [1:2645]';

% head attributes
hattr = struct;
hattr={ {'header' 'pltfid' 'ASL'}, ...
        {'header' 'instid' 'ERA5'}, ...
        {'header' 'pfields' '[1=profile,2=calc,4=obs (val=sum)]'}, ...
        {'header' 'ptype' '[0=levels,1=layers,2=AIRS layers]'}};
  
% generate prof structure
nlev  = length(era.plevs);
nsam  = length(era.sdnum);

prof = struct;
prof.nlevs  = ones(1,nsam) * 37;
prof.plat   = ones(1,nsam) * scen.lat;
prof.plon   = ones(1,nsam) * scen.lon;
prof.rtime  = dnum2tai(era.sdnum)';
prof.spres  = ones(1,nsam) + 1.037e5;;
prof.stemp  = scen.skt';
prof.ptemp  = scen.tz;
prof.gas_1  = scen.q;
prof.gas_3  = scen.o3;
prof.plevs  = repmat(era.plevs,1,nsam);
prof.salti  = zeros(1,nsam);
prof.clwc   = scen.clwc;
prof.ciwc   = scen.ciwc;
prof.cc     = scen.cc;
%
%prof.satzen = [];
%prof.solzen = [];
%prof.efreq  = [];
%prof.nemis  = [];
%prof.rho    = [];

% prof attributes
pattr = struct;
pattr={{'profiles' 'iudef(1,:)' 'Dust flag:[1=true,0=false,-1=land,-2=cloud,-3=bad data]'},...
       {'profiles' 'iudef(2,:)' 'Dust_score:[>380 (probable), N/A if Dust Flag < 0]'}};

% --------------


savdr = '/home/chepplew/data/rates_anomalies/tiled/era5_mon/rtp/';
if(~exist(savdr))
  disp(['Creating ' savdr]);
  mkdir(savdr);
end

fnrtp = ['era5_tile_center_lonbin'  ...
         sprintf('%02d_latbin%02d',era.lonbin,era.latbin) '.in.rtp'];


savfn = [savdr fnrtp];

rtpwrite(savfn, head, hattr, prof, pattr);

iok = 1;

% hdfml('listinfo')
% hdfml('closeall')
 
src_pwd = pwd;
cd(savdr)

klayers_exec = '/asl/packages/klayersV205/BinV201/klayers_airs_wetwater';
fn1 = fnrtp;
fn2 = ['era5_tile_center_lonbin'  ...
         sprintf('%02d_latbin%02d',era.lonbin,era.latbin) '.op.rtp'];


klayers_run = [klayers_exec ' fin=' fn1 ' fout=' fn2 ' > ' ...
               '/home/chepplew/logs/klayers/klout.txt'];

disp('Running klayers')
try
  unix(klayers_run);
catch  ME
  fprintf(1, '%s\n', ME.message);
  %continue;
end

cd(src_pwd)



%}
