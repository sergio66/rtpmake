% trend_co2_fill.m
%
% Sanity check for long term trend of CO2 from fill_co2.,
%        compared with simple MLO trend values.
%
% Method: Take a day of 1% AIRS Obs, set time stamp to date early in 
% sequence, apply fill_co2, do the calc, incrememt time stamp in monthly
% steps. Repeat for 20 years. 
% Do stats for global, zonal and spot location (MLO).
%

addpath /asl/matlib/h4tools
addpath /asl/matlib/time
addpath /home/chepplew/gitLib/rtp_fill_gases

% set up latitude bands
latB  = [-90:15:90];
latC  = 0.5*(latB(1:end-1) + latB(2:end));
nlats = length(latB);

% MLO location:
mlo.lat = 19.48;   % deg North 
mlo.lon = -155.58  % West. 

% Test file (from existing)
fn0 = '/asl/rtp/airs/airs_l1c_v672/random/2002/era_airicrad_day276_random.rtp';

[head,hattr,prof,~] = rtpread(fn0);
nobs = length(prof.rlat);

% set month and year
ref_dvec  = datevec(tai2dnum(prof.rtime) );
new_dvec  = ref_dvec;

% Loop of months/years
k = 1;
for iyr = 2003:2022       % 2003:2022
  for imn = 1:12         % 1:12      was mod([9:20],12)+1]
    new_dvec(:,1)   = iyr*ones(nobs,1);
    new_dvec(:,2)   = imn*ones(nobs,1);
    prof.rtime      = dnum2tai(datenum(new_dvec) );
%
    [head, hattr, prof1] = fill_co2(head,hattr,prof);
    [head, hattr, prof1] = fill_ch4(head,hattr,prof1);
    [head, hattr, prof1] = fill_n2o(head,hattr,prof1);
% --------------------------------------------------------
% Calculate global, zonal average profiles and MLO profile
% --------------------------------------------------------
    gas2.gbl(k,:) = mean(prof1.gas_2,2);
    gas4.gbl(k,:) = mean(prof1.gas_4,2);
    gas6.gbl(k,:) = mean(prof1.gas_6,2);
    clear iiwnt gas2.zon;
    for i=1:nlats-1
      iiwnt{i} = find(prof.rlat > latB(i) & prof.rlat <= latB(i+1));
      gas2.zon(k,:,i) = mean(prof1.gas_2(:,iiwnt{i}),2);
      gas6.zon(k,:,i) = mean(prof1.gas_6(:,iiwnt{i}),2);
      gas4.zon(k,:,i) = mean(prof1.gas_4(:,iiwnt{i}),2);
    end
% ---------------------
% find proximity to MLO
% ---------------------
    clear iiwnt
    iiwnt = find(prof.rlat > mlo.lat-2 & prof.rlat < mlo.lat+2 ...
              & prof.rlon > mlo.lon-2 & prof.rlon < mlo.lon+2);
    if(isempty(iiwnt)) 
      error('No profiles near MLO');
      disp(['iyr, imn: ' num2str(iyr) ' ' num2str(imn)])
      return;
    end

    iiwnt = median(iiwnt);
    gas2.mlo(k,:) = prof1.gas_2(:,iiwnt);

% Save a copy of time stamp
    gas2.tai(k) = median(prof1.rtime);

    k = k + 1;
    disp(num2str(k));
  end      % end imn for-loop

end      % end iyr for-loop

% Save data:
sav.dir = '/home/chepplew/projects/data/validation/trend_bias/';
sav.fname = 'trend_monthly_rtp_fill_gases.mat'';


%{
% SFC is end of vector
plot(tai2dtime(gas2.tai),gas2.mlo(:,58),'.-')

%}

