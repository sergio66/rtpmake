function [yhd0,ypdmat] = get_richardson_number_layers_v2(h0,p0,xhd0,xpdmat,iPlot)

%% assumes hd0,pd0 are for LAYERS profile, and
%%         xhd0,xpdmat came from corresponding NWP eg [xhd0,xpdmat] = get_richardson_number_levels(hd0,ha0,pd0,pa0);   so are in LEVELS formar
%% all salti are in meters
%% all PBLH are  in meters
%%
%% compared to get_richardson_number_levels_v3.m in the Bulk Richardson numerotor, this uses skt instead of t2m to set the Ri (tv(z) - tvs) instead of (tv(z) - tv2m)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
see_long_comment_about_PBLH_from_Ri_vs_Tpot.txt
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
VERSION 1

JPSS-1 2024/11/13
g180 : Antartica and Southern Ocean mix
g200 : mostly ocean off California
g209 : mostly Nepal and India, some Arabian Sea
g210 : mostly Indian Ocean

gran = 180;
%% airs_l1c ---> cris_fsr
loader = ['load /home/sergio/nogit/sergio_temp_rtp_files/j1_ccast_hires/allfov/2024/11/13//interp_analysis_uvw_cloudy_cris_fsr_ecm_sarta_baum_ice.2024.11.13.' num2str(gran) '.mat']; eval(loader)
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
VERSION 2

JPSS-1 2024/11/13
g180 : Antartica and Southern Ocean mix
g200 : mostly ocean off California
g209 : mostly Nepal and India, some Arabian Sea
g210 : mostly Indian Ocean

%% airs_l1c ---> cris_fsr
gran = 180;
wspeedfile = ['/home/sergio/nogit/sergio_temp_rtp_files/j1_ccast_hires/allfov/2024/11/13//interp_analysis_uvw_cloudy_cris_fsr_ecm_sarta_baum_ice.2024.11.13.' num2str(gran) '.mat'];
loader = ['load /home/sergio/nogit/sergio_temp_rtp_files/singlefootprintretrievals_ccast_hires_jpss1/2024/11/13/'];
  loader = [loader '/interp_analysis_ecm_retr' num2str(gran) '_cris_-1_iDET_4_iStemp_ColWV_21_iCenterFov_-1_iCO2_Yes_No_Switch_-1_singlelayerclouds.mat']; eval(loader)

addpath /umbc/rs/pi_sergio/WorkDirDec2025/matlabcode/PBL_Retrievals/HALO_BdryLayer/PBL_Hgt_from_poemNew/
[pumbc,pecm,pbulkRI] = driver_compute_PBLH_poemNew(hoemNew,poemNew,wspeedfile);
%}
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
Claude n'est pas Cote code
Here's a summary of what changed get_richardson_number_layers.m --> get_richardson_number_layers_v2.m :
Lines 426-448 (Ri threshold search):

length(good) > 10 → length(good) >= 1 — was requiring 10+ crossings which is way too strict
Added comment explaining why good(end) is correct for descending zalts
Refined the fine-grid search to also cap at 4km (was missing this before)
Fallback now sets proper NaN without the bogus salti + 10 placeholder

Lines 503-521 (post-compute_pblh adjustment):

Detects convective regime via dTvpot_surf < -5K
Convective (flag=2): falls back to zPBLH_Tpot — physically the best method for unstable BLs
Wind shear (flag=3): keeps the min|Ri-RiCritical| fallback from your original code
Good crossing (flag=0): no change

New PBLH_Ri_flag — lets you filter by method in downstream analysis:
%}
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
%if max(xpdmat.zPBLH_Ri) > 100
%  %% they are in km!!!!!!!    
%  xpdmat.zPBLH_Ri = xpdmat.zPBLH_Ri/1000;
%  xpdmat.salti    = xpdmat.salti/1000;      
%end

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

ypdmat.palts   = y2pd0.palts;
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
ypdmat.stemp    = y2pd0.stemp;   %% use input (probably retrieved) stemp

ypdmat.wspeed   = xpdmat.wspeed;
ypdmat.u10      = xpdmat.u10;
ypdmat.v10      = xpdmat.v10;

ypdmat.mmw = mmwater_rtp(h0,p0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iFixST_ECM = -1;
iFixST_ECM = +2;

if iFixST_ECM  == 1
  disp('warning ... using ECMWF stemp')
  ypdmat.stemp    = xpdmat.stemp;
elseif iFixST_ECM == 2
  boo = find(abs(p0.rlat) <= 60 & p0.landfrac == 0);
  if length(boo) > 0
    disp('   >>>> warning ... using ECMWF stemp in get_richardson_number_layers_v2.m for ocean/rlat < 60 <<<<<<')
      ypdmat.stemp(boo)    = xpdmat.stemp(boo);
    disp('   >>>> warning ... using ECMWF stemp in get_richardson_number_layers_v2.m for ocean/rlat < 60 <<<<<<')
  end      
end      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ypdmat.u        = pd0.u;
%ypdmat.v        = pd0.v;
%ypdmat.w        = pd0.w;     

if ~isfield(ypdmat,'plays')
  ypdmat.plays = plevs2plays(ypdmat.plevs);
end

%[size(ypdmat.plevs) size(ypdmat.palts)]
for ii = 1 : length(p0.stemp)
  n2lays = ypdmat.nlevs(ii)-1;
  p2lays = ypdmat.plays(1:n2lays,ii);
  t2temp = ypdmat.ptemp(1:n2lays,ii);
  
  nlays = p0.nlevs(ii)-1;
  numer = p0.spres(ii) - p0.plevs(nlays,ii);
  denom = p0.plevs(nlays+1,ii) - p0.plevs(nlays,ii);  
  ypdmat.fracL(ii) = numer/denom;
  junk = p0.spres(ii) - p0.plevs(nlays,ii);
  junk = junk/log(p0.spres(ii)/p0.plevs(nlays,ii));
  
  ypdmat.oldPlaysL(ii)   = ypdmat.plays(nlays,ii);
  ypdmat.newPlaysL(ii)   = junk;
  ypdmat.plays(nlays,ii) = junk;

  ypdmat.oldPtempL(ii)   = p0.ptemp(nlays,ii);
  ypdmat.ptemp(nlays,ii) = interp1(log(p2lays),t2temp,log(junk),[],'extrap');
  ypdmat.newPtempL(ii)   = ypdmat.ptemp(nlays,ii);
  
  ypdmat.oldPlevsLp1(ii) = p0.plevs(nlays+1,ii);
  ypdmat.oldPaltsLp1(ii) = p0.palts(nlays+1,ii);  
  ypdmat.plevs(nlays+1,ii) = p0.spres(ii);
  ypdmat.palts(nlays+1,ii) = p0.salti(ii)+0.1;

  ypdmat.gas_1(nlays,ii) = ypdmat.gas_1(nlays,ii) * ypdmat.fracL(ii);
end
%[size(ypdmat.plevs) size(ypdmat.palts)]

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
see ../COMMON_SETTINGS/testing_layers2gg_layers2sphum_conversions.m
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ggLAY1,ppmvLAY,ppmvAVG,ppmvMAX,pavgLAY,tavgLAY,ppmv500,ppmv75,ppmvSURF] = layers2gg(y2hd0,ypdmat,1:length(ypdmat.stemp),1);

dpjunk = zeros(101,length(p0.stemp));
dpjunk(1:100,:) = diff(p0.plevs,1)*100;  %% change mb to Pa
ggLAY2 = molecules_to_mixratio(ypdmat.gas_1,abs(dpjunk),ypdmat.nlevs-1,1);

dLjunk = zeros(101,length(p0.stemp));
dLjunk(1:100,:) = abs(diff(p0.palts,1));      %% in meteres
[ggLAY2,qqLAY2] = recover_q_from_forward(p0.gas_1,p0.plays*100,p0.ptemp,abs(dLjunk),p0.nlevs-1);  %% gg is mass mix ratio while qq is sp. humidity

[mmjunk,nnjunk] = size(ggLAY1);
if mmjunk < (max(p0.nlevs)-1)
  size(ggLAY1)
  (max(p0.nlevs)-1)
  error('sizess do not jive')
end  
ggLAY = ggLAY1;             %% mine
ggLAY = ggLAY2(1:mmjunk,:); %% claude

[mmx,nnx] = size(ggLAY);
ypdmat.gg = nan(size(ypdmat.ptemp));
ypdmat.gg(1:mmx,:) = ggLAY;

[rhLAY] = layeramt2RH(y2hd0,ypdmat);

[mmx,nnx] = size(ggLAY);
[mmxx,nnxx] = size(rhLAY);
ypdmat.rh = nan(size(ypdmat.ptemp));
ypdmat.rh(1:min(mmx,mmxx),:) = rhLAY(1:min(mmx,mmxx),:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% more /home/sergio/git/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/klayersV205_140levs/Doc/gas_units_code.txt
%         21   mass mixing ratio in (g/g) or (kg/kg), dry air
%              Grams of gas X per gram of "dry air"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% convert to potential temp

Rd_Cp = 0.286;      %% Rd/Cp
P0    = 1000;       %% mb
g     = 9.81;       %% m/s2

ypdmat.ptemp_pot = ypdmat.ptemp  .* ((P0./ypdmat.plevs).^Rd_Cp);            %% use air temp
ypdmat.Tvirtual  = ypdmat.ptemp .* (1 + 0.61 * ypdmat.gg);                  %% Tvirtual = T (1 + 0.61 r) where r is mix ratio in g/g or kg/kg
ypdmat.Tvirtual_potential = ypdmat.Tvirtual .* ((P0./ypdmat.plays).^Rd_Cp); %% uses virtual temp

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

set_Ri_critical

ypdmat.zalts = nan(size(y2pd0.plevs));
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
end  

%size(ypdmat.Tvirtual)
[mmm,nnn] = size(ypdmat.zalts);
ypdmat.staticE = nan(size(ypdmat.Tvirtual));
if iVers_Ri == 1
  %%%%% DIFFERENCE v2 vs v3  
  cp = 1004.7;   %J/kg/K
  ypdmat.surf_staticE     = cp * ypdmat.surf_Tvirtual     + g * ypdmat.salti;
  ypdmat.staticE(1:mmm,:) = cp * ypdmat.Tvirtual(1:mmm,:) + g * ypdmat.zalts(1:mmm,:);
end

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

speed_sqr = (ypdmat.u).^2 + (ypdmat.v).^2;                           %% WRONG before 5/20/26 but retry here, no diff from above
speed_sqr = (ypdmat.u - ypdmat.u10).^2 + (ypdmat.v - ypdmat.v10).^2; %% hmm this gives a 0.5 km bias in PBLH!
ypdmat.speed_sqr = speed_sqr;

[mm,nn] = size(ypdmat.gas_1);  %% 101 x 12150 for klayers

ypdmat.Ri           = nan(size(ypdmat.gas_1));
ypdmat.lapserate    = nan(size(ypdmat.gas_1));
ypdmat.PBLH_Ri_flag = zeros(1, nn);
%% PBLH_Ri_flag values:
%%   0 = good Ri threshold crossing found within 4km
%%   1 = no Ri crossing found (set before compute_pblh, may be overwritten)
%%   2 = convective BL (Tvpot_surf >> Tvpot_air): fell back to zPBLH_Tpot
%%   3 = wind shear case: Ri crosses >4km, capped using min|Ri-RiCritical|

for ii = 1 : nn
  %% these are after klayers
  nlevs = y2pd0.nlevs(ii);
  plevs = y2pd0.plevs(1:nlevs,ii);
  zalts = y2pd0.palts(1:nlevs,ii);

  nlays = y2pd0.nlevs(ii)-1;
  plays = y2pd0.plays(1:nlays,ii);
  zalts = ypdmat.zalts(1:nlays,ii);
  
  s2    = speed_sqr(1:nlays,ii);
  tp    = ypdmat.Tvirtual_potential(1:nlays,ii);
  tps   = ypdmat.stemp_pot(ii);

  if iVers_Ri == 0
    %% simmpler one
    %% ypdmat.Ri(1:nlays,ii) = (tp - tps) .* (zalts - ypdmat.salti(ii));  %% Siedel has (tp-tps)*(z-zs)
    %% WRONG before 5/20/26 : Bulk Ri method of Vogelezang and Holtslag [1996], referenced in 
    %% Dian Siedel (2012), Climatology of the planetary boundary layer over the
    %% continentalUnited States and Europe, JOURNAL OF GEOPHYSICAL
    %% RESEARCH, VOL. 117, D17106, doi:10.1029/2012JD018143, 2012    
    ypdmat.Ri(1:nlays,ii) = (tp - tps) .* (zalts);  %% Davy only has (tp-tps)*(z)
    ypdmat.Ri(1:nlays,ii) = g ./ tps ./s2 .* ypdmat.Ri(1:nlays,ii);
  elseif iVers_Ri == 1
    %% harder one, with static energy
    %% https://www.ecmwf.int/sites/default/files/elibrary/2017/17736-part-iv-physical-processes.pdf#section.3.10
    %% pg 50 of IFS DOCUMENTATION – Cy43r3 Operational implementation 11 July 2017, PART IV: PHYSICAL PROCESSES
    %%   3.10.1 Diagnostic boundary layer height, eqn 3.90

    %%%%% DIFFERENCE v2 vs v3      
    ypdmat.Ri(1:nlays,ii) = (ypdmat.staticE(1:nlays,ii) + ypdmat.surf_staticE(ii)) - g * (ypdmat.zalts(1:nlays,ii) + ypdmat.salti(ii));    
    ypdmat.Ri(1:nlays,ii) = (ypdmat.staticE(1:nlays,ii) - ypdmat.surf_staticE(ii)) ./ ypdmat.Ri(1:nlays,ii);
    ypdmat.Ri(1:nlays,ii) = 2 * g ./ s2 .* ypdmat.Ri(1:nlays,ii) .* (ypdmat.zalts(1:nlays,ii) - ypdmat.salti(ii));    
  end
  
  numer = diff(ypdmat.ptemp(1:nlays,ii));      %% dT  [K]
  denom = diff(zalts/1000);                    %% dz [km]  	       
  xlapse = numer./denom;                       %% K/km
  ypdmat.lapserate(1:nlays,ii) = -interp1(log(meanvaluebin(plays)),xlapse,log(plays),[],'extrap');   %% environment lapse rate

  %% stability
  ypdmat.stability(1:nlays,ii) = nan(nlays,1);
  boo = find(ypdmat.lapserate(1:nlays,ii) < MALR);                               ypdmat.stability(boo,ii) = -1;  %% absolutely stable, Rising air is colder than its surroundings and sinks, regardless of moisture content.
                                                                                                                 %% Typical in inversion layers.
  boo = find(ypdmat.lapserate(1:nlays,ii) > DALR);                               ypdmat.stability(boo,ii) = +1;  %% absolutely unstable, Rising air warmer than surroundings, continues to rise, forming convective clouds (e.g., cumulus).
  boo = find(MALR <= ypdmat.lapserate(1:nlays,ii) & ypdmat.lapserate(1:nlays,ii) <= DALR); ypdmat.stability(boo,ii) = 0;   %% conditionally unstable, Stable if air is unsaturated, but unstable if forced to saturation.
  boo = find(abs(DALR - ypdmat.lapserate(1:nlays,ii)) <= 0.01);                  ypdmat.stability(boo,ii) = -2;  %% neutral,  A lifted parcel stays at the new altitude

  levels_n    = xpdmat.nlevs(ii);
  levels_alts = xpdmat.zalts(1:levels_n,ii);
  levels_pres = xpdmat.plevs(1:levels_n,ii);

  dzalts = ypdmat.zalts(1:nlays,ii) - ypdmat.salti(ii);  
  wah = ypdmat.Ri(1:nlays,ii);

  %% NOTE: zalts is DESCENDING (TOA first, surface last)
  %% so layers near surface have LARGE indices
  %% "first crossing from surface upward" = good(end) in descending order
  good = find(wah >= RiCritical & dzalts/1000 <= 4);  

  if length(good) >= 1
    %% good(end) = lowest layer (closest to surface) where Ri >= RiCritical
    %% this is the correct "first crossing from surface upward" for descending zalts
    good_coarse = good(end);
    ypdmat.zPBLH_Ri_coarse(ii) = zalts(good_coarse);
    ypdmat.pPBLH_Ri_coarse(ii) = interp1(zalts,plays,ypdmat.zPBLH_Ri_coarse(ii),[],'extrap');

    %% refine on finer levels grid
    wah_fine = interp1(zalts,ypdmat.Ri(1:nlays,ii),levels_alts,[],'extrap');
    dz_fine  = levels_alts - ypdmat.salti(ii);
    good_fine = find(wah_fine >= RiCritical & dz_fine/1000 <= 4);
    if ~isempty(good_fine)
      good_fine = good_fine(end);   %% first crossing from surface in descending levels_alts
      ypdmat.zPBLH_Ri(ii) = levels_alts(good_fine);
    else
      ypdmat.zPBLH_Ri(ii) = ypdmat.zPBLH_Ri_coarse(ii);
    end
    ypdmat.pPBLH_Ri(ii) = interp1(levels_alts,levels_pres,ypdmat.zPBLH_Ri(ii),[],'extrap');
    ypdmat.PBLH_Ri_flag(ii) = 0;   %% good Ri crossing found

  else
    %% No Ri crossing found within 4km - convective/unstable profile
    ypdmat.zPBLH_Ri_coarse(ii) = NaN;
    ypdmat.pPBLH_Ri_coarse(ii) = NaN;
    ypdmat.zPBLH_Ri(ii)        = NaN;
    ypdmat.pPBLH_Ri(ii)        = NaN;
    ypdmat.PBLH_Ri_flag(ii)    = 1;   %% no Ri crossing - convective/unstable

  end

  moo = abs(plays -xpdmat.pPBLH_Ri(ii));
  moo = find(moo == min(moo),1);
  ypdmat.lapserate_at_pPBLH_Ri(ii) = ypdmat.lapserate(moo,ii);
  ypdmat.stability_at_pPBLH_Ri(ii) = ypdmat.stability(moo,ii);    
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  wah = interp1(zalts,ypdmat.gg(1:nlays,ii),levels_alts,[],'extrap');  %% find min gradient in z
  wah = gradient(wah,levels_alts);
  %good = find(wah == min(wah) & (levels_alts-ypdmat.salti(ii))/1000 <= 4.0,1);
  gah = find((levels_alts-ypdmat.salti(ii))/1000 <= 4);
  good = find(wah(gah) == min(wah(gah)),1);    
  ypdmat.zPBLH_gg(ii) = levels_alts(gah(good));
  ypdmat.pPBLH_gg(ii) = interp1(levels_alts,levels_pres,ypdmat.zPBLH_gg(ii),[],'extrap');
  
  wah = interp1(zalts,ypdmat.rh(1:nlays,ii),levels_alts,[],'extrap');  %% find min gradient in z
  wah = gradient(wah,levels_alts);
  %good = find(wah == min(wah) & (levels_alts-ypdmat.salti(ii))/1000 <= 4.0,1);
  gah = find((levels_alts-ypdmat.salti(ii))/1000 <= 4);
  good = find(wah(gah) == min(wah(gah)),1);    
  ypdmat.zPBLH_rh(ii) = levels_alts(gah(good));
  ypdmat.pPBLH_rh(ii) = interp1(levels_alts,levels_pres,ypdmat.zPBLH_rh(ii),[],'extrap');
  
  wah = interp1(zalts,ypdmat.Tvirtual(1:nlays,ii),levels_alts,[],'extrap');  %% find max gradient in z
  wah = gradient(wah,levels_alts);
  %good = find(wah == max(wah) & (levels_alts-ypdmat.salti(ii))/1000 <= 4.0,1);    
  gah = find((levels_alts-ypdmat.salti(ii))/1000 <= 4);
  good = find(wah(gah) == max(wah(gah)),1);
  ypdmat.zPBLH_Tvir(ii) = levels_alts(gah(good));
  ypdmat.pPBLH_Tvir(ii) = interp1(levels_alts,levels_pres,ypdmat.zPBLH_Tvir(ii),[],'extrap');

  wah = interp1(zalts,ypdmat.Tvirtual_potential(1:nlays,ii),levels_alts,[],'extrap');  %% find max gradient in z
  wah = gradient(wah,levels_alts);
  %good = find(wah == max(wah) & (levels_alts-ypdmat.salti(ii))/1000 <= 4.0,1);
  gah = find((levels_alts-ypdmat.salti(ii))/1000 <= 4);
  good = find(wah(gah) == max(wah(gah)),1);
  ypdmat.zPBLH_Tpot(ii) = levels_alts(gah(good));
  ypdmat.pPBLH_Tpot(ii) = interp1(levels_alts,levels_pres,ypdmat.zPBLH_Tpot(ii),[],'extrap');

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %% this is just a routine to compute PBLH from bulk_richardson
  pblhx = compute_pblh_Ri(ypdmat.zalts(1:nlays,ii),ypdmat.ptemp(1:nlays,ii),ypdmat.plays(1:nlays,ii),...
                         ypdmat.gg(1:nlays,ii),ypdmat.u(1:nlays,ii),ypdmat.v(1:nlays,ii),...
    		         ypdmat.stemp(ii),[ypdmat.u10(ii) ypdmat.v10(ii)],ypdmat.landfrac(ii));

  %%%%%%%%%%%%%%%%%%%%%%%%%

  ypdmat.zPBLH_Ri(ii) = pblhx;
  ypdmat.pPBLH_Ri(ii) = interp1(levels_alts,levels_pres,ypdmat.zPBLH_Ri(ii),[],'extrap');

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %%% now check if ypdmat.zPBLH_Ri <= 4 km above surface
  ypdmat.zPBLH_Ri_truecrossing_Rcrit(ii) = ypdmat.zPBLH_Ri(ii);
  ypdmat.pPBLH_Ri_truecrossing_Rcrit(ii) = ypdmat.pPBLH_Ri(ii);

  dTvpot_surf = ypdmat.Tvirtual_potential(nlays,ii) - ypdmat.stemp_pot(ii);
  bIsConvective = (dTvpot_surf < -5);   %% strongly unstable: surface >> air above

  if (ypdmat.zPBLH_Ri(ii) - ypdmat.salti(ii) <= 4000)
    %% Ri crossing within 4km above surface - good result, keep as-is
    ypdmat.PBLH_Ri_flag(ii) = 0;

  else
    %% Ri crossing > 4km or NaN - convective/unstable profile or wind shear issue
    %% see jpss gran 48, fov 12150, 2024/11/13 over Australia land
    %%   Ri starts out negative, starts swinging towards zero but at 3.8 km there
    %%   is a lot of wind shear and windspeed starts dropping, so Ri swings back
    %%   negative again. It only becomes >= 0.25 at about 8 km.

    dz = ypdmat.zalts(1:nlays,ii) - ypdmat.salti(ii);
    RR = ypdmat.Ri(1:nlays,ii);
    boo = find(dz <= 4000);

    if bIsConvective
      %% Convective BL: Ri is negative throughout, no meaningful threshold crossing
      %% Best fallback: use virtual potential temperature gradient method
      %% which is physically meaningful for both stable and unstable BLs
      ypdmat.zPBLH_Ri(ii) = ypdmat.zPBLH_Tpot(ii);   %% already in metres at this point
      ypdmat.pPBLH_Ri(ii) = ypdmat.pPBLH_Tpot(ii);
      ypdmat.PBLH_Ri_flag(ii) = 2;   %% convective: Ri failed, used Tpot gradient

    else
      %% Non-convective but Ri crosses threshold above 4km (wind shear case)
      %% Find where Ri is closest to RiCritical within 0-4km above surface
      RRX = abs(RR(boo) - RiCritical);
      moo = find(RRX == min(RRX));
      moo = moo(end);   %% take surface-side (end = lowest layer in descending zalts)
      ypdmat.zPBLH_Ri(ii) = ypdmat.zalts(boo(moo),ii);
      ypdmat.pPBLH_Ri(ii) = interp1(levels_alts,levels_pres,ypdmat.zPBLH_Ri(ii),[],'extrap');
      ypdmat.PBLH_Ri_flag(ii) = 3;   %% wind shear case: capped at 4km, min|Ri-Ric|
    end
  end
  
end

ypdmat.zPBLH_Tpot                  = ypdmat.zPBLH_Tpot/1000;
ypdmat.zPBLH_Tvir                  = ypdmat.zPBLH_Tvir/1000;
ypdmat.zPBLH_gg                    = ypdmat.zPBLH_gg/1000;
ypdmat.zPBLH_rh                    = ypdmat.zPBLH_rh/1000;
ypdmat.zPBLH_Ri_truecrossing_Rcrit = ypdmat.zPBLH_Ri_truecrossing_Rcrit/1000;  %% this is what you get from crossing Ri(z) at Ri_crtical
ypdmat.zPBLH_Ri                    = ypdmat.zPBLH_Ri/1000;       %% if ypdmat.zPBLH_Ri_truecrossing_Rcrit > 4, this is where Ri is minimum between 0-4 km
ypdmat.salti                       = ypdmat.salti/1000;

%figure(4); scatter_coast(ypdmat.rlon,ypdmat.rlat,50,ypdmat.zPBLH_Ri); caxis([0 4]); ax = axis; title('PBLH km')
%figure(5); scatter_coast(ypdmat.rlon,ypdmat.rlat,50,ypdmat.stemp); title('SKT')
%figure(6); scatter_coast(ypdmat.rlon,ypdmat.rlat,50,ypdmat.mmw);   title('mmw'); caxis([0 30])
%ix = 96; figure(7); scatter_coast(ypdmat.rlon,ypdmat.rlat,50,sqrt(ypdmat.u(ix,:).^2 + ypdmat.v(ix,:).^2)); axis(ax); title('Speed')
%keyboard_nowindow

if nargin == 4
  iPlot = -1;
end
if iPlot > 0
  plot_richardson_PBLH_layers
end
