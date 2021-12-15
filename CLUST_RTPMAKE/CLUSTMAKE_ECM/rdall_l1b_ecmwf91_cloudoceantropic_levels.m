function [ierr] = rdall_l1b_ecmwf(grannum,granfile,ecmwffile,rtpfile,idwant0);

% function [ierr] =rdall_l1b_ecmwf(grannum,granfile,ecmwffile,rtpfile,idwant0);
%
% Reads in L1B granule data for all FOVs and the corresponding
% nearest ECMWF profiles; FILTERS for over ocean tropic 
% and saves it to an output RTP file.
%
% Input:
%    grannum   = integer(1x1), granule number
%    granfile  = string, name of L1B granule file
%    ecmwffile = string, name of ECMWF file
%    rtpfile   = string, name of RTP output file to create
%    idwant0   = optional (nwant x 1) wanted channel IDs for output; all
%       channels are output if idwant0 is not specified
%
% Output:
%    ierr = integer(1x1), error detection flag [0=no error, 1=error] 
%

% Created: 15 January 2003, Scott Hannon
% Update: 11 April 2005, Scott Hannon - add optional argument idwant0
% Update: 09 May 2006, S.Hannon - change readecmwf91_nearest.m ECMWF reader
% Update: 22 Jan 2007, S. Machado uses readecmwf91_nearest_gasNcloud
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Edit this section only if needed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ierr=1;

% AIRS instrument code number
AIRSinst=800;

addpath /asl/matlab/airs/readers
addpath /asl/matlab/gribtools
addpath /asl/matlab/h4tools
addpath /asl/matlab/rtptools
addpath /asl/matlab/science


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check idwant0
idwant = (1:2378)';   %%default

iX = -1;
if iX > 0
  if (nargin == 5)
     d = size(idwant0);
     if (length(d) ~= 2 | min(d) ~= 1)
        disp('Error, idwant0 must be a (nwant x 1) array')
        return
     end
     nwant = length(idwant0);
     idwant = idwant0;
  else
     idwant0 = (1:2378)';
     idwant = idwant0;
  end
end

nwant = length(idwant);
whos nwant idwant0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make an RTP file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read the AIRS granule data file
fprintf(1,'granfile = %s \n',granfile)
[meantime, f, prof] = readl1b_all(granfile);
%disp(['Mean time is ' num2str(meantime) ' Z hours'])
nobs=length(prof.rlat);
nchan=length(f);

if (nobs ~= 12150)
   disp(['NOTE! readl1b_all returned ' int2str(nobs) ' FOVs, not 12150'])
end

fprintf(1,'ecmwffile = %s \n',ecmwffile)
% Read the ECMWF model and pull out nearest profile for each observation
addpath /home/sergio/MATLABCODE/CLOUD_ECMWF/PACKAGE
disp('moron 1')
[head, pprof] = ...
  readecmwf91_nearest_gasNcloud_levelsprof(ecmwffile, prof.rlat, prof.rlon);
head
pprof
disp('moron 2')

% Add channel info
head.nchan=nchan;
head.ichan=(1:nchan)';
head.vchan=f;  % approximate frequency
% Update pfields
head.pfields=5; % (1=prof + 4=IRobs);


% Fields found in both L1B and ECMWF
%       salti: [1x1350 double]
%    landfrac: [1x1350 double]
% For now, stick with the AIRS values
prof.udef      = zeros(2,nobs);
prof.udef(1,:) = prof.salti;
prof.udef(2,:) = prof.landfrac;
salti    = prof.salti;
landfrac = prof.landfrac;
%sumx = sum(prof.udef(2,:) - prof.landfrac);
%fprintf(1,'1 : sum(udef(2,:) - lf) = %8.6f \n',sumx)

% Append ECMWF profile info to observations
pfields1 = fieldnames(pprof);
npfields1 = length(pfields1);
for j=1:npfields1;
   fname = pfields1{j};
   eval(sprintf('[m,n] = size(pprof.%s);', fname));
   if n ~= nobs
      fname
      error('ECMWF field has an unexpected number of columns');
   end
   eval(['prof.' fname '=pprof.' fname ';'])
end
%% make sure we have the right surface values from L1B, not from ECMWF
prof.udef(1,:) = salti;
prof.udef(2,:) = landfrac;
prof.landfrac  = landfrac; 
clear pprof

%sumx = sum(prof.udef(2,:) - prof.landfrac);
%fprintf(1,'2 : sum(udef(2,:) - lf) = %8.6f \n',sumx)


% Add standard observation info
%
prof.upwell = ones(1,nobs); % radiance is upwelling
prof.pobs   = zeros(1,nobs);
prof.irinst = AIRSinst*ones(1,nobs);
prof.findex = grannum*ones(1,nobs);

% leave cfrac as is
% prof.cfrac=zeros(1,nobs);
%

% Plug in sea surface emissivity & reflectivity
[nemis,efreq,seaemis]=cal_seaemis2(prof.satzen,prof.wspeed);
prof.nemis = nemis;
prof.efreq = efreq;
prof.emis  = seaemis;
prof.nrho  = nemis;
prof.rfreq = efreq;
prof.rho   = (1-seaemis)/pi;

%sumx = sum(prof.udef(2,:) - prof.landfrac);
%fprintf(1,'3 : sum(udef(2,:) - lf) = %8.6f \n',sumx)

clear f seaemis efreq nemis

% attribute string for robs1 data
ii=max( find(granfile == '/') );
if (length(ii) == 0)
   ii=0;
end
junk=granfile((ii+1):length(granfile));
robs1_str=['airibrad file=' junk];

% attribute string for profile
ii=max( find(ecmwffile == '/') );
if (length(ii) == 0)
   ii=0;
end
junk=ecmwffile((ii+1):length(ecmwffile));
prof_str=['nearest ECMWF file=' junk];


% Assign RTP attribute strings
hattr={ {'header' 'profile' prof_str} };


pattr={ {'profiles' 'robs1' robs1_str}, ...
        {'profiles' 'udef' '1=L1B salti, 2=L1B landfrac'} };

clear robs1_str prof_str uniform_str

%sumx = sum(prof.udef(2,:) - prof.landfrac);
%fprintf(1,'4 : sum(udef(2,:) - lf) = %8.6f \n',sumx)

% Subset the channels if not all are to be output
if (nwant < 2378)
   [head, prof] = subset_rtp(head, prof, [], idwant, []);
end


%%% now do the subsetting
%sumx1 = sum(prof.udef(2,:) - prof.landfrac);
%fprintf(1,'5 : sum(udef(2,:) - lf) = %8.6f \n',sumx1)

il1b_oceantropic = find((prof.udef(2,:) <= 0.01) & (abs(prof.rlat) < 50));
fprintf(1,'number of OCEAN,SUBTROPIC == %5i \n',length(il1b_oceantropic))

il1b_oceantropic = find((prof.landfrac <= 0.01) & (abs(prof.rlat) < 50));
fprintf(1,'number of OCEAN,SUBTROPIC == %5i \n',length(il1b_oceantropic))

if (length(il1b_oceantropic) >= 1)
  disp('here we are 1');
  [head, prof] = subset_rtp(head, prof, [], idwant, il1b_oceantropic);

  % Write to an RTP file
  disp('here we are 2');
  rtpwrite(rtpfile, head, hattr, prof, pattr)
  fprintf(1,'wrote to rtpfile  = %s \n',rtpfile)

  ierr = 0;
else
  disp('sorry .. did not find any FOVS over ocean tropics');
  ierr = +1;
  end

%%% end of file %%%
