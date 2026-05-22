function [xhd0,xpdmat] = get_richardson_number_levels(hd0,ha0,pd0,pa0);

%% assumes hd0,pd0 are for LEVELS profile
%% easiest to use this in eg clustbatch_make_eraORecm_cloudrtp_sergio_sarta_YYMMDD_loopGG.m
%% all salti are in meters
%% all PBLH are  in meters

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
%{
addpath0

[hd0,ha0,pd0,pa0] = rtpread('/home/sergio/git/matlabcode/REGR_PROFILES_SARTA/REGR49_PROFILES_for_kCARTA_breakouts_for_SARTA/regr49_1013_385ppm.ip.rtp');
hd0.pfields = 7;
hd0.nchan = 1;
hd0.ichan = 1291;
hd0.vchan = 1231.2;
pd0.robs1 = ones(size(pd0.stemp))*ttorad(1231,290);
pd0.rcalc = ones(size(pd0.stemp))*ttorad(1231,290);
[xhd0,xpdmat] = get_richardson_number_levels(hd0,ha0,pd0,pa0);
%}
  
  
if hd0.ptype > 0
  error('need hd0.ptype == 0 == LEVELS')
end

fip = mktempS('richardson.ip.rtp');
fop = mktempS('richardson.op.rtp');

i1231 = find(hd0.vchan >= 1231,1);
i1231 = hd0.ichan(i1231);
[x1hd0,x1pd0] = subset_rtp_allcloudfields(hd0,pd0,[],i1231,1:length(pd0.stemp));
rtpwrite(fip,x1hd0,ha0,x1pd0,pa0);

set_path_to_execs
klayerser = ['!' klayers ' fin=' fip ' fout=' fop ' >& ugh'];
eval(klayerser);

[x2hd0,ha0,x2pd0,pa0] = rtpread(fop);
xhd0.vchan = x2hd0.vchan;
xhd0.ichan = x2hd0.ichan;
xhd0.nchan = x2hd0.nchan;
xhd0.ngas  = 1;
xhd0.glist = hd0.glist(1);
xhd0.gunit = hd0.gunit(1);

rmer = ['!/bin/rm ' fip ' ' fop];
eval(rmer);

%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isfield(pd0,'u')
  disp('putting in rubbish u,v,w')
  pd0.w = zeros(size(pd0.ptemp));
  pd0.v = zeros(size(pd0.ptemp));
  pd0.u = ones(size(pd0.ptemp)) * 10;
  pd0.u10 = ones (size(pd0.ptemp)) * 8;
  pd0.v10 = ones (size(pd0.ptemp)) * 8;    
end

xpdmat.robs1   = x1pd0.robs1;
xpdmat.rcalc   = x1pd0.rcalc;

xpdmat.nlevs   = x1pd0.nlevs;
xpdmat.plevs   = x1pd0.plevs;
xpdmat.ptemp   = x1pd0.ptemp;
xpdmat.gas_1   = x1pd0.gas_1;

xpdmat.atrack  = x2pd0.atrack;
xpdmat.xtrack  = x2pd0.xtrack;

xpdmat.scanang = x2pd0.scanang;
xpdmat.solzen  = x2pd0.solzen;
xpdmat.rlon    = x2pd0.rlon;
xpdmat.rlat    = x2pd0.rlat;

xpdmat.landfrac = x2pd0.landfrac;
xpdmat.spres    = x2pd0.spres;
xpdmat.salti    = x2pd0.salti;
xpdmat.stemp    = x2pd0.stemp;

%% wind velocity (u,v,w) at 91 levels
xpdmat.u        = pd0.u;
xpdmat.v        = pd0.v;
xpdmat.w        = pd0.w;     

%% wind velocity (u,v,w) at 10 m
xpdmat.wspeed   = x2pd0.wspeed;
xpdmat.u10       = pd0.u10;
xpdmat.v10       = pd0.v10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% more /home/sergio/git/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/klayersV205_140levs/Doc/gas_units_code.txt
%         21   mass mixing ratio in (g/g) or (kg/kg), dry air
%              Grams of gas X per gram of "dry air"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% convert to potential temp

Rd_Cp = 0.286;      %% Rd/Cp
P0    = 1000;       %% mb
g     = 9.81;       %% m/s2

if hd0.glist(1) ~= 1
  error('need gasID = 1');
end
if hd0.gunit(1) ~= 21
  error('need gunit = 21 for gasID = 1');
end

xpdmat.Tvirtual = xpdmat.ptemp .* (1 + 0.61 * xpdmat.gas_1);   %% Tvirtual = T (1 + 0.61 r) where r is mix ratio in g/g or kg/kg
xpdmat.ptemp_pot = xpdmat.Tvirtual .* ((P0./xpdmat.plevs).^Rd_Cp);

nn = length(xpdmat.stemp);
for ii = 1 : nn
  mm = pd0.nlevs(ii);
  xpdmat.surf_Tvirtual(ii) = xpdmat.stemp(ii) .* (1 + 0.61 * xpdmat.gas_1(mm,ii));   %% Tvirtual = T (1 + 0.61 r) where r is mix ratio in g/g or kg/kg
  xpdmat.stemp_pot(ii)     = xpdmat.surf_Tvirtual(ii) .* ((P0./xpdmat.spres(ii)).^Rd_Cp);
end

set_Ri_critical

for ii = 1 : nn
  %% these are after klayers
  nlevs = x2pd0.nlevs(ii);
  plevs = x2pd0.plevs(1:nlevs,ii);
  zalts = x2pd0.palts(1:nlevs,ii);
  
  %% these are before klayers ie raw sonde or NWP
  NNlevs = pd0.nlevs(ii);
  PPlevs = pd0.plevs(1:NNlevs,ii);
  xpdmat.zalts(1:NNlevs,ii) = interp1(log(plevs),zalts,log(PPlevs),[],'extrap');
end  

if iVers_Ri == 1
  cp = 1004.7;   %J/kg/K
  xpdmat.surf_staticE = cp * xpdmat.surf_Tvirtual + g * xpdmat.salti;
  xpdmat.staticE      = cp * xpdmat.Tvirtual      + g * xpdmat.zalts;
end

% debug
% if iVers_Ri == 1
%   pcolor((xpdmat.staticE - cp*xpdmat.ptemp_pot)./xpdmat.staticE);
%   shading interp; colorbar
% 
%   plot(nanmean((xpdmat.staticE - cp*xpdmat.ptemp_pot)./xpdmat.staticE,2),nanmean(xpdmat.zalts,2)/1000)
%   axis([-1/2 +1/2 0 20])
%   title('frac error  staticE - cp*Tpot')
%   plotaxis2; xlabel('frac error'); ylabel('hgt km')
%   keyboard_nowindow  
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% convert to Ri =  g     (Tpot(z)-Tpot0)*(z-z0)
%%                 ---   ----------------------
%%                 Tpot0   horizspeed^2
%%
%% where Tpot0,zpot0 = potential temp at surface, surface altitude
%% units m/s2 K m / (K m2/s2) = m2/s2/(m2/s2) = []

% The bulk Richardson number (R_{ib}) can be negative near the
% surface, specifically in convective conditions where the surface is
% warmer than the air above, indicating unstable atmospheric conditions.
% It signifies that buoyancy forces dominate over wind shear in the
% lower boundary layer. For more details, visit skybrary.aero.

%   Dry Adiabatic Lapse Rate (DALR): Unsaturated air cools at a rate of
%   approximately 10K per 1 km (5.5 F per 10000 ft) of ascent

%   Moist Adiabatic Lapse Rate (MALR): Saturated air cools at a slower
%   rate because condensation releases latent heat, roughly 6 K per km
%   (3 F per 1000 ft)

MALR = 6;
DALR = 10;

%speed_sqr = (xpdmat.u).^2 + (xpdmat.v).^2;  WRONG before 5/20/26
speed_sqr = (xpdmat.u - xpdmat.u10).^2 + (xpdmat.v - xpdmat.v10).^2;
	       
[mm,nn] = size(xpdmat.gas_1);  %% 91 x 12150 for NWP, but for sondes can be any number  of levels, esp if this is regr49.iprtp or ecm83.ip.rtp

xpdmat.Ri    = nan(size(xpdmat.gas_1));
xpdmat.lapse = nan(size(xpdmat.gas_1));

for ii = 1 : nn
  %% these are after klayers
  nlevs = x2pd0.nlevs(ii);
  plevs = x2pd0.plevs(1:nlevs,ii);
  zalts = x2pd0.palts(1:nlevs,ii);

  %% these are before klayers ie raw sonde or NWP
  NNlevs = pd0.nlevs(ii);
  PPlevs = pd0.plevs(1:NNlevs,ii);
  %% xpdmat.zalts(1:NNlevs,ii) = interp1(log(plevs),zalts,log(PPlevs),[],'extrap');
  s2    = speed_sqr(1:NNlevs,ii);
  tp    = xpdmat.ptemp_pot(1:NNlevs,ii);
  tps   = xpdmat.stemp_pot(ii);
  
  if iVers_Ri == 0
    %% simpler one
    %% xpdmat.Ri(1:NNlevs,ii) = (tp - tps) .* (xpdmat.zalts(1:NNlevs,ii) - xpdmat.salti(ii));  %% Siedel has (tp-tps)*(z-zs)
    %% WRONG before 5/20/26 : Bulk Ri method of Vogelezang and Holtslag [1996], referenced in 
    %% Dian Siedel (2012), Climatology of the planetary boundary layer over the
    %% continentalUnited States and Europe, JOURNAL OF GEOPHYSICAL
    %% RESEARCH, VOL. 117, D17106, doi:10.1029/2012JD018143, 2012    
    xpdmat.Ri(1:NNlevs,ii) = (tp - tps) .* (xpdmat.zalts(1:NNlevs,ii));  %% Davy only has (tp-tps)*(z)
    xpdmat.Ri(1:NNlevs,ii) = g ./ tps ./s2 .* xpdmat.Ri(1:NNlevs,ii);
  elseif iVers_Ri == 1
    %% harder one, with static energy
    %% https://www.ecmwf.int/sites/default/files/elibrary/2017/17736-part-iv-physical-processes.pdf#section.3.10
    %% pg 50 of IFS DOCUMENTATION – Cy43r3 Operational implementation 11 July 2017, PART IV: PHYSICAL PROCESSES
    %%   3.10.1 Diagnostic boundary layer height, eqn 3.90    
    xpdmat.Ri(1:NNlevs,ii) = (xpdmat.staticE(1:NNlevs,ii) + xpdmat.surf_staticE(ii)) - g * (xpdmat.zalts(1:NNlevs,ii) + xpdmat.salti(ii));
    xpdmat.Ri(1:NNlevs,ii) = (xpdmat.staticE(1:NNlevs,ii) - xpdmat.surf_staticE(ii)) ./ xpdmat.Ri(1:NNlevs,ii);
    xpdmat.Ri(1:NNlevs,ii) = 2 * g ./ s2 .* xpdmat.Ri(1:NNlevs,ii) .* (xpdmat.zalts(1:NNlevs,ii) - xpdmat.salti(ii));

    %axpdmat.Ri(1:NNlevs,ii) = (tp - tps) .* (xpdmat.zalts(1:NNlevs,ii) - xpdmat.salti(ii));  %% Davy only has (tp-tps)*(z)
    %axpdmat.Ri(1:NNlevs,ii) = g ./ tps ./s2 .* axpdmat.Ri(1:NNlevs,ii);
    %plot(xpdmat.Ri(1:NNlevs,ii),xpdmat.zalts(1:NNlevs,ii)/1000,'bx-',axpdmat.Ri(1:NNlevs,ii),xpdmat.zalts(1:NNlevs,ii)/1000,'rx-'); axis([-40 120 0 10])    
    %keyboard_nowindow
  end

  numer = diff(xpdmat.ptemp(1:NNlevs,ii));      %% dT  [K]
  denom = diff(xpdmat.zalts(1:NNlevs,ii)/1000); %% dz [km]  	       
  xlapse = numer./denom;                 %% K/km
  xpdmat.lapse(1:NNlevs,ii) = -interp1(log(meanvaluebin(PPlevs)),xlapse,log(PPlevs),[],'extrap');   %% environment lapse rate

  %% stability
  xpdmat.stable(1:NNlevs,ii) = nan(NNlevs,1);
  boo = find(xpdmat.lapse(1:NNlevs,ii) < MALR);                               xpdmat.stable(boo,ii) = -1;  %% absolutely stable, Rising air is colder than its surroundings and sinks, regardless of moisture content. Typical in inversion layers.
  boo = find(xpdmat.lapse(1:NNlevs,ii) > DALR);                               xpdmat.stable(boo,ii) = +1;  %% absolutely unstable, Rising air is warmer than its surroundings and continues to rise, forming convective clouds (e.g., cumulus).
  boo = find(MALR <= xpdmat.lapse(1:NNlevs,ii) & xpdmat.lapse(1:NNlevs,ii) <= DALR); xpdmat.stable(boo,ii) = 0;   %% conditionally unstable, Stable if air is unsaturated, but unstable if forced to saturation.
  boo = find(abs(DALR - xpdmat.lapse(1:NNlevs,ii)) <= 0.01);                  xpdmat.stable(boo,ii) = -2;  %% neutral,  A lifted parcel stays at the new altitude

  %if iVers_Ri == 1
  %  plot(xpdmat.Ri(1:NNlevs,ii),xpdmat.zalts(1:NNlevs,ii)/1000); title([num2str(ii) ' ' num2str(xpdmat.salti(ii))] ); axis([0 1 0 5]); disp('ret to continue'); pause
  %end
  
  wah = xpdmat.Ri(1:NNlevs,ii);
  good = find(wah >= RiCritical);
  good = good(end);
  xpdmat.zPBLH_Ri(ii) = xpdmat.zalts(good,ii);
  xpdmat.pPBLH_Ri(ii) = interp1(xpdmat.zalts(1:NNlevs,ii),PPlevs,xpdmat.zPBLH_Ri(ii),[],'extrap');
end

xpdmat.zPBLH_Ri   = xpdmat.zPBLH_Ri/1000;
xpdmat.salti      = xpdmat.salti/1000;

iPlot = -1;
if iPlot > 0
  plot_richardson_PBLH_levels
end

