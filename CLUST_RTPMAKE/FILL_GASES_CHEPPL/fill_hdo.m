function [head, prof, pattr] = fill_hdo(head, prof, pattr)

%
% fill_hdo: apply a variable HDO depletion to prof.udef(20,:)
%           with units per.mil e.g. -600 p.mil 
%           values limited to -999 to 999.
% Notes: prof.udef MAXUDEF=20 float32 (prof.iudef: MAXIUDEF=10)
% Update pattr with definition for udef(20,:)
%
%
%

addpath /asl/matlib/rtptools                  % set_attr

disp('fill_hdo.m: adding HDO depletion')

% add hdo_depletion factor to udef(20)
nprofs = length(prof.rlat);
iun = 20;
hdo_depl = -600.0;  % -680.0;
udef20 = hdo_depl*ones(1,nprofs);
prof.udef(iun,1:nprofs) = udef20; 

% add attribute string
pattr = set_attr(pattr, 'udef(20,:)','HDO depletion. Units: per.mil');

% No change to head

% ------ end of function -----------------




