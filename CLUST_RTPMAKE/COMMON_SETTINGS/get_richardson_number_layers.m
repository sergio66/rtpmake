function [yhd0,ypdmat] = get_richardson_number_layers(h0,p0,xhd0,xpdmat)

%% assumes hd0,pd0 are for LAYERS profile, and
%%          xhd0,xpdmat came from corresponding NWP eg [xhd0,xpdmat] = get_richardson_number_levels(hd0,ha0,pd0,pa0);   so are in LEVELS formar

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{

g180 : Antartica and Southern Ocean mix
g200 : mostly ocean off California
g209 : mostly Nepal and India, some Arabian Sea
g210 : mostly Indian Ocean

load /home/sergio/nogit/sergio_temp_rtp_files/j1_ccast_hires/allfov/2024/11/13//interp_analysis_uvw_cloudy_airs_l1c_ecm_sarta_baum_ice.2024.11.13.200.mat   %% this gives levels_h0,levels_pdmat0
load /home/sergio/nogit/sergio_temp_rtp_files/singlefootprintretrievals_ccast_hires_jpss1/2024/11/13/interp_analysis_ecm_retr200_cris_-1_iDET_4_iStemp_ColWV_21_iCenterFov_-1_iCO2_Yes_No_Switch_-1_singlelayerclouds.mat

gran = 180;
loader = ['load /home/sergio/nogit/sergio_temp_rtp_files/j1_ccast_hires/allfov/2024/11/13//interp_analysis_uvw_cloudy_airs_l1c_ecm_sarta_baum_ice.2024.11.13.' num2str(gran) '.mat']; eval(loader)
loader = ['load /home/sergio/nogit/sergio_temp_rtp_files/singlefootprintretrievals_ccast_hires_jpss1/2024/11/13/'];
  loader = [loader '/interp_analysis_ecm_retr' num2str(gran) '_cris_-1_iDET_4_iStemp_ColWV_21_iCenterFov_-1_iCO2_Yes_No_Switch_-1_singlelayerclouds.mat']; eval(loader)

pecm = poemNew;
pecm.stemp = poemNew.stemp_OEMinitialization;
pecm.stemp = xpdmat.stemp;
pecm.ptemp = poemNew.ptemp_OEMinitialization;
pecm.gas_1 = poemNew.gas_1_OEMinitialization;

[yhd0,ypdmat] = get_richardson_number_layers(hoemNew,poemNew,xhd0,xpdmat)
[zhd0,zpdmat] = get_richardson_number_layers(hoemNew,pecm,xhd0,xpdmat)
%}
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
To determine the Planetary Boundary Layer Height PBLH or H
from the bulk Richardson number Ri_{B}, you must find the height
at which the cumulative Ri_{B} from the surface exceeds a critical
value, typically Ri_{Bc} \approx 0.25.

The method involves interpolating between model levels where
Ri_{B} crosses this threshold.Steps to Determine PBLH
HCalculate Ri_{B} Profiles:

Compute the bulk Richardson number from the surface z_{0} up to
several heights z using:Ri_{B}(z)=\frac{g}{\theta
_{vs}}\frac{(\theta _{v}(z)-\theta
_{vs})(z-z_{0})}{u(z)^{2}+v(z)^{2}}g: Gravity 9.81 \,
m/s^2\theta _{vs}: Virtual potential temperature at the
surface\theta_v(z): Virtual potential temperature at height
zu, v: Horizontal wind components

Define Critical Value Ri_{Bc}: Select a critical Richardson
number. While 0.25 is common, studies suggest using 0.24 for
strong stable layers, 0.31 for weak stable, and 0.39 for
unstable conditions for better accuracy.

Identify the Height: Find the lowest level z where Ri_B(z) \ge Ri_{Bc}

Interpolate: Linearly interpolate between the height below the
threshold and the height above it to find the exact height H

Context-Specific ConsiderationsUrban Areas: When dealing with urban
roughness sublayers, the formula is modified to account for
displacement height d and surface parameters.

Stable vs. Unstable: For stable boundary layers, Ri_{Bc} can range
from 0.5 to 3.0 depending on the specific study
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
if h0.ptype < 1
  error('need h0.ptype == 1,2 == LAYERS')
end

i1231 = find(h0.vchan >= 1231,1);
i1231 = h0.ichan(i1231);
[y2hd0,y2pd0] = subset_rtp_allcloudfields(h0,p0,[],i1231,1:length(p0.stemp));

yhd0.vchan = y2hd0.vchan;
yhd0.ichan = y2hd0.ichan;
yhd0.nchan = y2hd0.nchan;
yhd0.ngas  = 1;
yhd0.glist = h0.glist(1);
yhd0.gunit = h0.gunit(1);

%%%%%%%%%%%%%%%%%%%%%%%%%

ypdmat.robs1   = y2pd0.robs1;
ypdmat.rcalc   = y2pd0.rcalc;

ypdmat.plevs   = y2pd0.plevs;
ypdmat.nlevs   = y2pd0.nlevs;
ypdmat.ptemp   = y2pd0.ptemp;
ypdmat.gas_1   = y2pd0.gas_1;

ypdmat.atrack  = y2pd0.atrack;
ypdmat.xtrack  = y2pd0.xtrack;

ypdmat.scanang = y2pd0.scanang;
ypdmat.solzen  = y2pd0.solzen;
ypdmat.rlon    = y2pd0.rlon;
ypdmat.rlat    = y2pd0.rlat;

ypdmat.landfrac = y2pd0.landfrac;
ypdmat.spres    = y2pd0.spres;
ypdmat.salti    = y2pd0.salti;
ypdmat.stemp    = y2pd0.stemp;
ypdmat.wspeed   = y2pd0.wspeed;

%disp('warning ... using ECMWF stemp')
%ypdmat.stemp    = xpdmat.stemp;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ypdmat.u        = pd0.u;
%ypdmat.v        = pd0.v;
%ypdmat.w        = pd0.w;     
if ~isfield(ypdmat,'plays')
  ypdmat.plays = plevs2plays(ypdmat.plevs);
end

ypdmat.u = nan(size(ypdmat.ptemp));
ypdmat.v = nan(size(ypdmat.ptemp));
ypdmat.w = nan(size(ypdmat.ptemp));

if ~isfield(xpdmat,'nlevs')
  [mmjunk,nnjunk] = size(xpdmat.ptemp);
  xpdmat.nlevs = mmjunk * ones(size(xpdmat.ptemp));
end

for ii = 1 : length(ypdmat.stemp)
  n2lays = ypdmat.nlevs(ii)-1;
  p2lays = ypdmat.plays(1:n2lays,ii);
  
  n1levs = xpdmat.nlevs(ii);
  p1levs = xpdmat.plevs(1:n1levs,ii);
  u1levs = xpdmat.u(1:n1levs,ii);
  v1levs = xpdmat.v(1:n1levs,ii);
  w1levs = xpdmat.w(1:n1levs,ii);
  
  ypdmat.u(1:n2lays,ii) = interp1(log(p1levs),u1levs,log(p2lays),[],'extrap');
  ypdmat.v(1:n2lays,ii) = interp1(log(p1levs),v1levs,log(p2lays),[],'extrap');
  ypdmat.w(1:n2lays,ii) = interp1(log(p1levs),w1levs,log(p2lays),[],'extrap');  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if y2hd0.glist(1) ~= 1
  error('need gasID = 1');
end
if y2hd0.gunit(1) ~= 1
  error('need gunit = 1 for gasID = 1 ---> then we convert to gunit 21 g/g');
end

[ggLAY,ppmvLAY,ppmvAVG,ppmvMAX,pavgLAY,tavgLAY,ppmv500,ppmv75,ppmvSURF] = layers2gg(y2hd0,ypdmat,1:length(ypdmat.stemp),1);

[mmx,nnx] = size(ggLAY);
ypdmat.gg = nan(size(ypdmat.ptemp));
ypdmat.gg(1:mmx,:) = ggLAY;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% more /home/sergio/git/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/klayersV205_140levs/Doc/gas_units_code.txt
%         21   mass mixing ratio in (g/g) or (kg/kg), dry air
%              Grams of gas X per gram of "dry air"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% convert to potential temp

Rd_Cp = 0.286;      %% Rd/Cp
P0    = 1000;       %% mb
g     = 9.81;       %% m/s2

ypdmat.Tvirtual = ypdmat.ptemp .* (1 + 0.61 * ypdmat.gg);   %% Tvirtual = T (1 + 0.61 r) where r is mix ratio in g/g or kg/kg
ypdmat.ptemp_pot = ypdmat.Tvirtual .* ((P0./ypdmat.plays).^Rd_Cp);

nn = length(ypdmat.stemp);
for ii = 1 : nn
  mm = p0.nlevs(ii) - 1;
  plevs  = ypdmat.plevs(1:mm,ii);
  gg     = ypdmat.gg(1:mm,ii);
  ggSurf = interp1(log(plevs),gg,log(ypdmat.spres(ii)),[],'extrap');
  ypdmat.surf_Tvirtual(ii) = ypdmat.stemp(ii) .* (1 + 0.61 * ypdmat.gg(mm,ii));   %% Tvirtual = T (1 + 0.61 r) where r is mix ratio in g/g or kg/kg
  ypdmat.surf_Tvirtual(ii) = ypdmat.stemp(ii) .* (1 + 0.61 * ggSurf);             %% Tvirtual = T (1 + 0.61 r) where r is mix ratio in g/g or kg/kg  
  ypdmat.stemp_pot(ii)     = ypdmat.surf_Tvirtual(ii) .* ((P0./ypdmat.spres(ii)).^Rd_Cp);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% convert to Ri =  g     (Tpot(z)-Tpot0)*(z-z0)
%%                 ---   ----------------------
%%                 Tpot0   horizspeed^2
%%
%% where Tpot0,zpot0 = potential temp at surface, surface altitude
%% units m/s2 K m / (K m2/s2) = m2/s2/(m2/s2) = []

% Select a critical Richardson number. While 0.25 is common, studies
% suggest using 0.24 for strong stable layers, 0.31 for weak stable,
% and 0.39 for unstable conditions for better accuracy.

RiCritical = 0.27;

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

speed_sqr = (ypdmat.u).^2 + (ypdmat.v).^2;
	       
[mm,nn] = size(ypdmat.gas_1);  %% 101 x 12150 for klayers

ypdmat.Ri    = nan(size(ypdmat.gas_1));
ypdmat.zalts = nan(size(ypdmat.gas_1));

ypdmat.lapse = nan(size(ypdmat.gas_1));

for ii = 1 : nn
  %% these are after klayers
  nlevs = y2pd0.nlevs(ii);
  plevs = y2pd0.plevs(1:nlevs,ii);
  zalts = y2pd0.palts(1:nlevs,ii);

  nlays = y2pd0.nlevs(ii)-1;
  plays = y2pd0.plays(1:nlays,ii);
  zalts = interp1(log(plevs),zalts,log(plays),[],'extrap');
  ypdmat.zalts(:,ii) = nan;
  ypdmat.zalts(1:nlays,ii) = zalts;  
  
  s2    = speed_sqr(1:nlays,ii);
  tp    = ypdmat.ptemp_pot(1:nlays,ii);
  tps   = ypdmat.stemp_pot(ii);
  ypdmat.Ri(1:nlays,ii) = (tp - tps) .* (zalts - ypdmat.salti(ii));
  ypdmat.Ri(1:nlays,ii) = g ./ tps ./s2 .* ypdmat.Ri(1:nlays,ii);

  numer = diff(ypdmat.ptemp(1:nlays,ii));      %% dT  [K]
  denom = diff(zalts/1000);                    %% dz [km]  	       
  xlapse = numer./denom;                       %% K/km
  ypdmat.lapse(1:nlays,ii) = -interp1(log(meanvaluebin(plays)),xlapse,log(plays),[],'extrap');   %% environment lapse rate

  %% stability
  ypdmat.stable(1:nlays,ii) = nan(nlays,1);
  boo = find(ypdmat.lapse(1:nlays,ii) < MALR);                               ypdmat.stable(boo,ii) = -1;  %% absolutely stable, Rising air is colder than its surroundings and sinks, regardless of moisture content. Typical in inversion layers.
  boo = find(ypdmat.lapse(1:nlays,ii) > DALR);                               ypdmat.stable(boo,ii) = +1;  %% absolutely unstable, Rising air is warmer than its surroundings and continues to rise, forming convective clouds (e.g., cumulus).
  boo = find(MALR <= ypdmat.lapse(1:nlays,ii) & ypdmat.lapse(1:nlays,ii) <= DALR); ypdmat.stable(boo,ii) = 0;   %% conditionally unstable, Stable if air is unsaturated, but unstable if forced to saturation.
  boo = find(abs(DALR - ypdmat.lapse(1:nlays,ii)) <= 0.01);                  ypdmat.stable(boo,ii) = -2;  %% neutral,  A lifted parcel stays at the new altitude
  
  wah = ypdmat.Ri(1:nlays,ii);
  good = find(wah >= RiCritical);
  good = good(end);
  ypdmat.zPBLH_Ri_coarse(ii) = zalts(good);
  ypdmat.pPBLH_Ri_coarse(ii) = interp1(zalts,plays,ypdmat.zPBLH_Ri_coarse(ii),[],'extrap');

  levels_n    = xpdmat.nlevs(ii);
  levels_alts = xpdmat.zalts(1:levels_n,ii);
  levels_pres = xpdmat.plevs(1:levels_n,ii);
  wah = interp1(zalts,ypdmat.Ri(1:nlays,ii),levels_alts,[],'extrap');
  good = find(wah >= RiCritical);
  good = good(end);
  ypdmat.zPBLH_Ri(ii) = levels_alts(good);
  ypdmat.pPBLH_Ri(ii) = interp1(levels_alts,levels_pres,ypdmat.zPBLH_Ri(ii),[],'extrap');  
end

ocean = find(xpdmat.landfrac == 0);
land  = find(xpdmat.landfrac > 0.90);
if length(ocean)/length(xpdmat.landfrac) < 0.01
  ocean = land;
elseif length(land)/length(xpdmat.landfrac) < 0.01
  land = ocean;
end

figure(1); clf; scatter_coast(xpdmat.rlon,xpdmat.rlat,50,xpdmat.zPBLH_Ri);                   title('PBLH Ri from ECMWF');  cx = caxis;
figure(2); clf; scatter_coast(xpdmat.rlon,xpdmat.rlat,50,ypdmat.zPBLH_Ri);                   title('PBLH Ri from UMBC');   caxis(cx);
figure(3); clf; scatter_coast(xpdmat.rlon,xpdmat.rlat,50,xpdmat.zPBLH_Ri-ypdmat.zPBLH_Ri);   title('PBLH Ri ECMWF-UMBC'); colormap(usa2); caxis([-1 +1]*1500)
figure(4); clf; dz = -2000 : 100 : +2000; plot(dz,histc(xpdmat.zPBLH_Ri-ypdmat.zPBLH_Ri,dz));title('PBLH Ri ECMWF-UMBC'); grid;
  fprintf(1,'mean(NWP = UMBC) = %8.6f m, std(NWP - UMBC) = %8.6f \n',mean(xpdmat.zPBLH_Ri-ypdmat.zPBLH_Ri),std(xpdmat.zPBLH_Ri-ypdmat.zPBLH_Ri))

figure(5); clf; plot(-diff(xpdmat.zalts(:,6000))/1000,xpdmat.zalts(2:91,6000)/1000,'bo-',-diff(ypdmat.zalts(:,6000))/1000,ypdmat.zalts(2:101,6000)/1000,'rx-')
  axis([0 2 0 5]); ylabel('hgt [km]'); xlabel('laayer thickness or diff between levels [km]'); legend('NWP levels xpdmat','UMBC layers ypdmat');
figure(6); clf; plot(xpdmat.Ri(:,6000),xpdmat.zalts(:,6000)/1000,'bo-',ypdmat.Ri(:,6000),ypdmat.zalts(:,6000)/1000,'rx-')
  axis([-2 2 0 5]); ylabel('hgt [km]'); xlabel('Ri []'); plotaxis2; legend('NWP levels xpdmat','UMBC layers ypdmat');

figure(5); clf; plot(-diff(nanmean(xpdmat.zalts(:,ocean),2))/1000,nanmean(xpdmat.zalts(2:91,ocean),2)/1000,'bo-',-diff(nanmean(ypdmat.zalts(:,ocean),2))/1000,nanmean(ypdmat.zalts(2:101,ocean),2)/1000,'rx-')
  axis([0 2 0 5]); ylabel('hgt [km]'); xlabel('laayer thickness or diff between levels [km]'); legend('NWP levels xpdmat','UMBC layers ypdmat');
figure(6); clf; plot(nanmean(xpdmat.Ri(:,ocean),2),nanmean(xpdmat.zalts(:,ocean),2)/1000,'bo-',nanmean(ypdmat.Ri(:,ocean),2),nanmean(ypdmat.zalts(:,ocean),2)/1000,'rx-')
  axis([-2 2 0 5]); ylabel('hgt [km]'); xlabel('Ri []'); plotaxis2; legend('NWP levels xpdmat','UMBC layers ypdmat');

figure(7); clf; plot(nanmean(xpdmat.ptemp(:,ocean),2),nanmean(xpdmat.zalts(:,ocean),2)/1000,'bo-',nanmean(ypdmat.ptemp(:,ocean),2),nanmean(ypdmat.zalts(:,ocean),2)/1000,'rx-')
  axis([200 300 0 15]); ylabel('hgt [km]'); xlabel('T(z) [K]'); plotaxis2; legend('NWP levels xpdmat','UMBC layers ypdmat');
  axis([290 300 0 2])
figure(8); clf; semilogx(nanmean(xpdmat.gas_1(:,ocean),2),nanmean(xpdmat.zalts(:,ocean),2)/1000,'bo-',nanmean(ypdmat.gg(:,ocean),2),nanmean(ypdmat.zalts(:,ocean),2)/1000,'rx-')
  axis([0.001 0.1 0 15]); ylabel('hgt [km]'); xlabel('WV MR [g/g]'); plotaxis2; legend('NWP levels xpdmat','UMBC layers ypdmat');
  axis([0.01 0.02 0 2])

lev_speed = sqrt((xpdmat.u).^2 + (xpdmat.v).^2);
lay_speed = sqrt((ypdmat.u).^2 + (ypdmat.v).^2);
figure(9); clf; plot(nanmean(lev_speed(:,ocean),2),nanmean(xpdmat.zalts(:,ocean),2)/1000,'bo-',nanmean(lay_speed(:,ocean),2),nanmean(ypdmat.zalts(:,ocean),2)/1000,'rx-')
  axis([0 10 0 15]); ylabel('hgt [km]'); xlabel('windspeed [m/s]'); plotaxis2; legend('NWP levels xpdmat','UMBC layers ypdmat');
  axis([0 10 0 2])

  %{
figure(7); clf; semilogy(nanmean(xpdmat.ptemp(:,ocean),2),nanmean(xpdmat.plevs(:,ocean),2),'bo-',nanmean(ypdmat.ptemp(:,ocean),2),nanmean(ypdmat.plays(:,ocean),2),'rx-')
  axis([250 300 500 1050]); ylabel('hgt [km]'); xlabel('T(z) [K]'); plotaxis2; legend('NWP levels xpdmat','UMBC layers ypdmat'); set(gca,'ydir','reverse')
figure(8); clf; semilogx(nanmean(xpdmat.gas_1(:,ocean),2),nanmean(xpdmat.plevs(:,ocean),2),'bo-',nanmean(ypdmat.gg(:,ocean),2),nanmean(ypdmat.plays(:,ocean),2),'rx-')
  axis([0.001 0.1 500 1050]); ylabel('hgt [km]'); xlabel('WV MR [g/g]'); plotaxis2; legend('NWP levels xpdmat','UMBC layers ypdmat'); set(gca,'ydir','reverse')
%}
  
%iPlot = -1;
%if iPlot > 0
%  plot_richardson_PBLH
%end

