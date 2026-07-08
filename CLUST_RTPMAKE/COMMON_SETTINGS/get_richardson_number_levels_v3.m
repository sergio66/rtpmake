function [xhd0,xpdmat] = get_richardson_number_levels_v3(hd0,ha0,pd0,pa0,iUseCdrag);

%% assumes hd0,pd0 are for LEVELS profile
%% easiest to use this in eg clustbatch_make_eraORecm_cloudrtp_sergio_sarta_YYMMDD_loopGG.m
%% all salti are in meters --> km at the end
%% all PBLH are  in meters --> km at the end
%%
%% compared to get_richardson_number_levels_v2.m in the Bulk Richardosn numerotor, this uses t2m instead of skt to set the Ri (tv(z) - tv2m) instead of (tv(z) - tvs)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
see_long_comment_about_PBLH_from_Ri_vs_Tpot.txt
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{

JPSS-1 2024/11/13
g180 : Antartica and Southern Ocean mix
g200 : mostly ocean off California
g209 : mostly Nepal and India, some Arabian Sea
g210 : mostly Indian Ocean

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
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin  < 5
  %% see COMMON_SETTINGS/ECMWF_PBLH_ustar*.pdf
  iUseCdrag = -1;  %% ERA5 does NOT put in Crag = b ugnd^2 term
  %% see COMMON_SETTINGS/ECMWF_PBLH_ustar*.pdf  
end

ZMAX = 4.0;   %% highest PBLH = RiCrssoing-SALTI
ZMAX = 7.0;   %% highest PBLH = RiCrssoing-SALTI  
    
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

clear xpdmat
xpdmat = struct;
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
if isfield(x1pd0,'t2m')
  disp('get_richardson_number_levels_v3.m : yay t2m available')
  xpdmat.t2m      = x1pd0.t2m;
else
  disp('get_richardson_number_levels_v3.m : t2m unavailable; use lowest level')
  for ii = 1 : length(x1pd0.stemp)
    nlevs = x1pd0.nlevs(ii);
    if x1pd0.plevs(1,ii) < x1pd0.plevs(nlevs,ii)
      xpdmat.t2m(ii)      = x1pd0.ptemp(nlevs,ii);
    else
      xpdmat.t2m(ii)      = x1pd0.ptemp(1,ii);
    end
  end
end  
if isfield(pd0,'pblh_nwp')
  disp('yay found pd0.pblh_nwp')
  xpdmat.pblh_nwp = pd0.pblh_nwp;
end

%% wind velocity (u,v,w) at 91 levels
xpdmat.u        = pd0.u;
xpdmat.v        = pd0.v;
xpdmat.w        = pd0.w;     

%% wind velocity (u,v,w) at 10 m
xpdmat.wspeed   = x2pd0.wspeed;
xpdmat.u10      = pd0.u10;
xpdmat.v10      = pd0.v10;

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

xpdmat.ptemp_pot = xpdmat.ptemp .* ((P0./xpdmat.plevs).^Rd_Cp);          %% use air temp
xpdmat.Tvirtual  = xpdmat.ptemp .* (1 + 0.61 * xpdmat.gas_1);            %% Tvirtual = T (1 + 0.61 r) where r is mix ratio in g/g or kg/kg
xpdmat.Tvirtual_potential = xpdmat.ptemp .* ((P0./xpdmat.plevs).^Rd_Cp); %% uses virtual temp

nn = length(xpdmat.stemp);
for ii = 1 : nn
  mm = pd0.nlevs(ii);
  xpdmat.surf_Tvirtual(ii) = xpdmat.stemp(ii) .* (1 + 0.61 * xpdmat.gas_1(mm,ii));   %% Tvirtual = T (1 + 0.61 r) where r is mix ratio in g/g or kg/kg
  xpdmat.stemp_pot(ii)     = xpdmat.surf_Tvirtual(ii) .* ((P0./xpdmat.spres(ii)).^Rd_Cp);

  xpdmat.t2m_Tvirtual(ii)  = xpdmat.t2m(ii) .* (1 + 0.61 * xpdmat.gas_1(mm,ii));   %% Tvirtual = T (1 + 0.61 r) where r is mix ratio in g/g or kg/kg
  xpdmat.t2m_pot(ii)       = xpdmat.t2m_Tvirtual(ii) .* ((P0./xpdmat.spres(ii)).^Rd_Cp);
  
end

xpdmat.rh = levels2RH(xpdmat.gas_1,xpdmat.ptemp,xpdmat.plevs);

set_Ri_critical  %% sets RiCritical and iVers_Ri

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

  xpdmat.surf_staticE = cp * xpdmat.t2m_Tvirtual  + g * xpdmat.salti;
  xpdmat.staticE      = cp * xpdmat.Tvirtual      + g * xpdmat.zalts;
end

% debug
% if iVers_Ri == 1
%   pcolor((xpdmat.staticE - cp*xpdmat.Tvirtual_potential)./xpdmat.staticE);
%   shading interp; colorbar
% 
%   plot(nanmean((xpdmat.staticE - cp*xpdmat.Tvirtual_potential)./xpdmat.staticE,2),nanmean(xpdmat.zalts,2)/1000)
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
xspeed_sqr = (xpdmat.u).^2 + (xpdmat.v).^2;  %% assume no velocity at 2m or 10 m, since you cannot get this from radiosonde observations
	       
[mm,nn] = size(xpdmat.gas_1);  %% 91 x 12150 for NWP, but for sondes can be any number  of levels, esp if this is regr49.iprtp or ecm83.ip.rtp

xpdmat.Ri        = nan(size(xpdmat.gas_1));
xpdmat.lapserate = nan(size(xpdmat.gas_1));

for ii = 1 : nn
  %% these are after klayers
  nlevs = x2pd0.nlevs(ii);
  plevs = x2pd0.plevs(1:nlevs,ii);
  zalts = x2pd0.palts(1:nlevs,ii);

  %%%%%%%%%%%%%%%%%%%%%%%%%
  if iUseCdrag < 0
    Rib = 0.0;
    ustar = 0.0;
    bdrag = 0.0;
  else  
    if pd0.landfrac(ii) < 0.1
      %% oocean
      Rib = 0.0;
      ustar = 0.0;
    else
      dz = 10.0;  %% 10 m
      Tav = 0.5*(pd0.t2m(ii)+pd0.stemp(ii));
      Usqr = (pd0.u10(ii).^2 + pd0.v10(ii).^2);
      Rib = g * dz * (pd0.t2m(ii)-pd0.stemp(ii))/Tav/Usqr;
  
      kappa = 0.4;
      z0 = 0.1;  %% roughness lenght = 1-4 m in urban areas, 0005 in soil, 0.1 m in crops, 0.25 m in high crops etc
      CdragN = (kappa/log(dz/z0))^2;
  
      if Rib < 0
        %% unstable
        b = 5.0;
        c = 5.0;
        d = 5.0;
        Cdrag = 1 + 3 * b * c * CdragN * sqrt((dz/z0)*abs(Rib));
        Cdrag = CdragN * (1 - 2*b*Rib/Cdrag);
      elseif Rib >= 0
        Cdrag = 1 + 2 * b * Rib/sqrt(1 + d * Rib);
        Cdrag = CdragN /Cdrag;
      end
      ustar = sqrt(Cdrag) * sqrt(Usqr);
      bdrag = 100.0;      
    end
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%
  
  %% these are before klayers ie raw sonde or NWP
  NNlevs = pd0.nlevs(ii);
  PPlevs = pd0.plevs(1:NNlevs,ii);
  %% xpdmat.zalts(1:NNlevs,ii) = interp1(log(plevs),zalts,log(PPlevs),[],'extrap');
  s2    = speed_sqr(1:NNlevs,ii) + bdrag*ustar;
  tp    = xpdmat.Tvirtual_potential(1:NNlevs,ii);
  tps   = xpdmat.stemp_pot(ii);
  tp2   = xpdmat.t2m_pot(ii);
  
  if iVers_Ri == 0
    %% simpler one
    %% xpdmat.Ri(1:NNlevs,ii) = (tp - tps) .* (xpdmat.zalts(1:NNlevs,ii) - xpdmat.salti(ii));  %% Siedel has (tp-tps)*(z-zs)
    %% WRONG before 5/20/26 : Bulk Ri method of Vogelezang and Holtslag [1996], referenced in 
    %% Dian Siedel (2012), Climatology of the planetary boundary layer over the
    %% continentalUnited States and Europe, JOURNAL OF GEOPHYSICAL
    %% RESEARCH, VOL. 117, D17106, doi:10.1029/2012JD018143, 2012    
    xpdmat.Ri(1:NNlevs,ii) = (tp - tps) .* (xpdmat.zalts(1:NNlevs,ii));                     %% Davy only has (tp-tps)*(z)
    xpdmat.Ri(1:NNlevs,ii) = (tp - tps) .* (xpdmat.zalts(1:NNlevs,ii) - xpdmat.salti(ii));  %% Davy only has (tp-tps)*(z) but let me put (z-zs)    
    xpdmat.Ri(1:NNlevs,ii) = g ./ tps ./s2 .* xpdmat.Ri(1:NNlevs,ii);

    xpdmat.Ri(1:NNlevs,ii) = (tp - tp2) .* (xpdmat.zalts(1:NNlevs,ii));                     %% Davy only has (tp-tp2)*(z)
    xpdmat.Ri(1:NNlevs,ii) = (tp - tp2) .* (xpdmat.zalts(1:NNlevs,ii) - xpdmat.salti(ii));  %% Davy only has (tp-tp2)*(z) but let me put (z-zs)    
    xpdmat.Ri(1:NNlevs,ii) = g ./ tp2 ./s2 .* xpdmat.Ri(1:NNlevs,ii);

  elseif iVers_Ri == 1
    %% harder one, with static energy
    %% https://www.ecmwf.int/sites/default/files/elibrary/2017/17736-part-iv-physical-processes.pdf#section.3.10
    %% pg 50 of IFS DOCUMENTATION – Cy43r3 Operational implementation 11 July 2017, PART IV: PHYSICAL PROCESSES
    %%   3.10.1 Diagnostic boundary layer height, eqn 3.90    
    xpdmat.Ri(1:NNlevs,ii) = (xpdmat.staticE(1:NNlevs,ii) + xpdmat.surf_staticE(ii)) - g * (xpdmat.zalts(1:NNlevs,ii) + xpdmat.salti(ii));
    xpdmat.Ri(1:NNlevs,ii) = (xpdmat.staticE(1:NNlevs,ii) - xpdmat.surf_staticE(ii)) ./ xpdmat.Ri(1:NNlevs,ii);
    xpdmat.Ri(1:NNlevs,ii) = 2 * g ./ s2 .* xpdmat.Ri(1:NNlevs,ii) .* (xpdmat.zalts(1:NNlevs,ii) - xpdmat.salti(ii));

  elseif iVers_Ri == 2
    %% simpler one
    %% xpdmat.Ri(1:NNlevs,ii) = (tp - tps) .* (xpdmat.zalts(1:NNlevs,ii) - xpdmat.salti(ii));  %% Siedel has (tp-tps)*(z-zs)
    %% WRONG before 5/20/26 : Bulk Ri method of Vogelezang and Holtslag [1996], referenced in 
    %% Dian Siedel (2012), Climatology of the planetary boundary layer over the
    %% continentalUnited States and Europe, JOURNAL OF GEOPHYSICAL
    %% RESEARCH, VOL. 117, D17106, doi:10.1029/2012JD018143, 2012
    error('iuiuiui')
    xs2   = xspeed_sqr(1:NNlevs,ii) + bdrag*ustar;    
    dz = diff(zalts);
    dz = dz(2:NNlevs);
    numer = diff(tp);
    denom = diff(xs2);
    xpdmat.Ri(2:NNlevs,ii) = g ./ tp(2:NNlevs,ii) .* dz .* numer ./denom;
  end		 

  numer = diff(xpdmat.ptemp(1:NNlevs,ii));      %% dT  [K]
  denom = diff(xpdmat.zalts(1:NNlevs,ii)/1000); %% dz [km]  	       
  xlapse = numer./denom;                 %% K/km
  xpdmat.lapserate(1:NNlevs,ii) = -interp1(log(meanvaluebin(PPlevs)),xlapse,log(PPlevs),[],'extrap');   %% environment lapse rate

  %% stability
  xpdmat.stability(1:NNlevs,ii) = nan(NNlevs,1);
  boo = find(xpdmat.lapserate(1:NNlevs,ii) < MALR);                               xpdmat.stability(boo,ii) = -1;  %% absolutely stable, Rising air is colder than its surroundings and sinks, regardless of moisture content.
                                                                                                                  %% Typical in inversion layers.
  boo = find(xpdmat.lapserate(1:NNlevs,ii) > DALR);                               xpdmat.stability(boo,ii) = +1;  %% absolutely unstable, Rising air warmer than surroundings, continues to rise, forming convective clouds (e.g., cumulus).
  boo = find(MALR <= xpdmat.lapserate(1:NNlevs,ii) & xpdmat.lapserate(1:NNlevs,ii) <= DALR); xpdmat.stability(boo,ii) = 0;   %% conditionally unstable, Stable if air is unsaturated, but unstable if forced to saturation.
  boo = find(abs(DALR - xpdmat.lapserate(1:NNlevs,ii)) <= 0.01);                  xpdmat.stability(boo,ii) = -2;  %% neutral,  A lifted parcel stays at the new altitude

  %if iVers_Ri == 1
  %  plot(xpdmat.Ri(1:NNlevs,ii),xpdmat.zalts(1:NNlevs,ii)/1000); title([num2str(ii) ' ' num2str(xpdmat.salti(ii))] ); axis([0 1 0 5]); disp('ret to continue'); pause
  %end

  dzalts = xpdmat.zalts(1:NNlevs,ii) - xpdmat.salti(ii);
  wah = xpdmat.Ri(1:NNlevs,ii);

  %% NOTE: zalts is DESCENDING (TOA first, surface last)
  %% so layers near surface have LARGE indices
  %% "first crossing from surface upward" = good(end) in descending order  
  %good = find(wah >= RiCritical);  
  good1 = find(wah >= RiCritical & dzalts/1000 <= ZMAX & dzalts/1000 >= 0.20);
  if length(good1) >= 1
    %% this is the "first(lowest) crossing from surface upward" for descending zalts        
    good = good1(end);   %% zalts(1) = 80 km, zalts(end) = MSL so go as low as possible when   good = find(wah >= RiCritical & dzalts/1000 <= ZMAX & dzalts/1000 >= 0.20);
    %% this is the "last (highest) crossing from surface upward" for descending zalts           
    %% good = good2(1);  %% zalts(1) = 80 km, zalts(end) = MSL so go as high as possible when  good = find(wah <= RiCritical & dzalts/1000 <= ZMAX & dzalts/1000 >= 0.20);
    xpdmat.z1PBLH_Ri(ii) = xpdmat.zalts(good,ii);
    xpdmat.p1PBLH_Ri(ii) = interp1(xpdmat.zalts(1:NNlevs,ii),PPlevs,xpdmat.z1PBLH_Ri(ii),[],'extrap');
  else
    xpdmat.z1PBLH_Ri(ii) = xpdmat.salti(ii) + 10;
    xpdmat.p1PBLH_Ri(ii) = 0.99999 * xpdmat.spres(ii);
  end
  
  good2 = find(wah <= RiCritical & dzalts/1000 <= ZMAX & dzalts/1000 >= 0.20);  
  if length(good2) >= 1
    %% this is the "first(lowest) crossing from surface upward" for descending zalts        
    %good = good1(end);   %% zalts(1) = 80 km, zalts(end) = MSL so go as low as possible when   good = find(wah >= RiCritical & dzalts/1000 <= ZMAX & dzalts/1000 >= 0.20);
    %% this is the "last (highest) crossing from surface upward" for descending zalts           
    good = good2(1);  %% zalts(1) = 80 km, zalts(end) = MSL so go as high as possible when  good = find(wah <= RiCritical & dzalts/1000 <= ZMAX & dzalts/1000 >= 0.20);
    xpdmat.z2PBLH_Ri(ii) = xpdmat.zalts(good,ii);
    xpdmat.p2PBLH_Ri(ii) = interp1(xpdmat.zalts(1:NNlevs,ii),PPlevs,xpdmat.z2PBLH_Ri(ii),[],'extrap');
  else
    xpdmat.z2PBLH_Ri(ii) = xpdmat.salti(ii) + 10;
    xpdmat.p2PBLH_Ri(ii) = 0.99999 * xpdmat.spres(ii);
  end

  xpdmat.zPBLH_Ri(ii) = 0.5*(xpdmat.z1PBLH_Ri(ii) + xpdmat.z2PBLH_Ri(ii));
  xpdmat.pPBLH_Ri(ii) = 0.5*(xpdmat.p1PBLH_Ri(ii) + xpdmat.p2PBLH_Ri(ii));  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %% this is just a routine to compute PBLH from bulk_richardson
  %[pblhx,Ribx,junkx] = compute_pblh_Ri(xpdmat.zalts(1:NNlevs,ii),xpdmat.ptemp(1:NNlevs,ii),xpdmat.plevs(1:NNlevs,ii),...
  %                                     xpdmat.gas_1(1:NNlevs,ii),xpdmat.u(1:NNlevs,ii),xpdmat.v(1:NNlevs,ii),...
  %                                     xpdmat.stemp(ii),[xpdmat.u10(ii) xpdmat.v10(ii)],xpdmat.landfrac(ii));
  
  %% notice use of t2m instead of skt!!!
  [pblhx,Ribx,junkx] = compute_pblh_Ri(xpdmat.zalts(1:NNlevs,ii),xpdmat.ptemp(1:NNlevs,ii),xpdmat.plevs(1:NNlevs,ii),...
                                       xpdmat.gas_1(1:NNlevs,ii),xpdmat.u(1:NNlevs,ii),xpdmat.v(1:NNlevs,ii),...
                                       xpdmat.t2m(ii),[xpdmat.u10(ii) xpdmat.v10(ii)],xpdmat.landfrac(ii));
  
  xpdmat.zPBLH_Ri_from_compute_pblh_Ri(ii) = pblhx;
  xpdmat.pPBLH_Ri_from_compute_pblh_Ri(ii) = interp1(xpdmat.zalts(1:NNlevs,ii),PPlevs,xpdmat.zPBLH_Ri_from_compute_pblh_Ri(ii),[],'extrap');  

  %%%%%%%%%%%%%%%%%%%%%%%%%

  moo = abs(PPlevs -xpdmat.pPBLH_Ri(ii));
  moo = find(moo == min(moo),1);
  xpdmat.lapserate_at_pPBLH_Ri(ii) = xpdmat.lapserate(moo,ii);
  xpdmat.stability_at_pPBLH_Ri(ii) = xpdmat.stability(moo,ii);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  levels_alts = xpdmat.zalts(1:NNlevs,ii);
  levels_pres = xpdmat.plevs(1:NNlevs,ii);  

  %wah = interp1(zalts,xpdmat.gas_1(1:NNlevs,ii),levels_alts,[],'extrap');  %% find min gradient in z, gas_1 in levels is in g/g
  wah = xpdmat.gas_1(1:NNlevs,ii);
  wah = gradient(wah,levels_alts);
  %good = find(wah == min(wah) & (levels_alts-xpdmat.salti(ii))/1000 <= ZMAX,1);
  gah = find((levels_alts-xpdmat.salti(ii))/1000 <= ZMAX);
  good = find(wah(gah) == min(wah(gah)),1);    
  xpdmat.zPBLH_gg(ii) = levels_alts(gah(good));
  xpdmat.pPBLH_gg(ii) = interp1(levels_alts,levels_pres,xpdmat.zPBLH_gg(ii),[],'extrap');
  
  %wah = interp1(zalts,xpdmat.rh(1:NNlevs,ii),levels_alts,[],'extrap');  %% find min gradient in z
  wah = xpdmat.rh(1:NNlevs,ii);
  wah = gradient(wah,levels_alts);
  %good = find(wah == min(wah) & (levels_alts-xpdmat.salti(ii))/1000 <= ZMAX,1);
  gah = find((levels_alts-xpdmat.salti(ii))/1000 <= ZMAX);
  good = find(wah(gah) == min(wah(gah)),1);    
  xpdmat.zPBLH_rh(ii) = levels_alts(gah(good));
  xpdmat.pPBLH_rh(ii) = interp1(levels_alts,levels_pres,xpdmat.zPBLH_rh(ii),[],'extrap');
  
  %wah = interp1(zalts,xpdmat.Tvirtual(1:NNlevs,ii),levels_alts,[],'extrap');  %% find max gradient in z
  wah = xpdmat.Tvirtual(1:NNlevs,ii);
  wah = gradient(wah,levels_alts);
  %good = find(wah == max(wah) & (levels_alts-xpdmat.salti(ii))/1000 <= ZMAX,1);    
  gah = find((levels_alts-xpdmat.salti(ii))/1000 <= ZMAX);
  good = find(wah(gah) == max(wah(gah)),1);
  xpdmat.zPBLH_Tvir(ii) = levels_alts(gah(good));
  xpdmat.pPBLH_Tvir(ii) = interp1(levels_alts,levels_pres,xpdmat.zPBLH_Tvir(ii),[],'extrap');

  %wah = interp1(zalts,xpdmat.Tvirtual_potential(1:NNlevs,ii),levels_alts,[],'extrap');  %% find max gradient in z
  wah = xpdmat.Tvirtual_potential(1:NNlevs,ii);
  wah = gradient(wah,levels_alts);
  %good = find(wah == max(wah) & (levels_alts-xpdmat.salti(ii))/1000 <= ZMAX,1);
  gah = find((levels_alts-xpdmat.salti(ii))/1000 <= ZMAX);
  good = find(wah(gah) == max(wah(gah)),1);
  xpdmat.zPBLH_Tpot(ii) = levels_alts(gah(good));
  xpdmat.pPBLH_Tpot(ii) = interp1(levels_alts,levels_pres,xpdmat.zPBLH_Tpot(ii),[],'extrap');

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if ii == 12150
    fprintf(1,'tps vs junk.tps = %2f %2f \n',tps,junkx.tps)
    figure(1); plot(xpdmat.Ri(:,ii),levels_alts/1000,'bx-',flipud(Ribx),levels_alts/1000,'ro-'); ylim([0 6]); plotaxis2; title('Bulk Ri []');
      xlim([-1 +1]*50);
      legend('this fcn','compute\_pblh\_Ri','location','best')
      figure(2); plot(tp, levels_alts/1000,'bx-', flipud(junkx.tp), levels_alts/1000,'ro-'); title('Tvir_pot [k]'); ylim([0 6])
      legend('this fcn','compute\_pblh\_Ri','location','best')    
      figure(3); plot(s2,levels_alts/1000,'bx-',flipud(junkx.speed_sqr),levels_alts/1000,'ro-'); title('Speed^2'); ylim([0 6])
      legend('this fcn','compute\_pblh\_Ri','location','best')
    figure(4); plot(xpdmat.u(1:NNlevs,ii), levels_alts/1000,'bx-', flipud(junkx.u), levels_alts/1000,'ro-',...
                    xpdmat.v(1:NNlevs,ii), levels_alts/1000,'cx-', flipud(junkx.v), levels_alts/1000,'mo-');         title('u,v [m/s]')
    legend('this fcn','compute\_pblh\_Ri','location','best')

    figure(5);  scatter_coast(xpdmat.rlon,xpdmat.rlat,20,xpdmat.zPBLH_Ri_from_compute_pblh_Ri/1000-xpdmat.salti/1000); colormap jet; caxis([0 6]); title('from compute pblh ri')
    %figure(5);  scatter_coast(xpdmat.rlon,xpdmat.rlat,20,xpdmat.zPBLH_Ri_from_compute_pblh_Ri); %% colormap jet; caxis([0 6]); title('from compute pblh ri')    
    
    figure(3); scatter_coast(xpdmat.rlon,xpdmat.rlat,20,xpdmat.zPBLH_Ri/1000-xpdmat.salti/1000);       title('SERGIO PBLH [km] '); caxis([0 6])
    figure(7); scatter_coast(xpdmat.rlon,xpdmat.rlat,20,xpdmat.z1PBLH_Ri/1000-xpdmat.salti/1000);      title('SERGIO PBLH1 [km] '); caxis([0 6])
    figure(8); scatter_coast(xpdmat.rlon,xpdmat.rlat,20,xpdmat.z2PBLH_Ri/1000-xpdmat.salti/1000);      title('SERGIO PBLH2 [km] '); caxis([0 6])
    if isfield(pd0,'pblh_nwp')
      figure(9); scatter_coast(xpdmat.rlon,xpdmat.rlat,20,pd0.pblh_nwp);                               title('ERA5 or ECM raw PBLH [km] '); caxis([0 6])              
    end

    figure(1); colormap jet
    figure(2); colormap jet
    figure(3); colormap jet
    figure(4); colormap jet
    figure(5); colormap jet
    figure(6); colormap jet
    figure(7); colormap jet
    figure(8); colormap jet
    figure(9); colormap jet
    
    pause(0.1);
    %%error('kljsljslfjasljflsaj in get_richardson_number_levels_v2.m')
    %disp('see the PNLH plots : ret to continue'); pause
    %[xpdmat.zPBLH_Ri_from_compute_pblh_Ri(ii) xpdmat.zPBLH_Ri(ii) xpdmat.z1PBLH_Ri(ii) xpdmat.z2PBLH_Ri(ii)]
    %[xpdmat.zPBLH_Ri_from_compute_pblh_Ri(01) xpdmat.zPBLH_Ri(01) xpdmat.z1PBLH_Ri(01) xpdmat.z2PBLH_Ri(01)]    
    %keyboard_nowindow
  end
  
  %%% now check if xpdmat.zPBLH_Ri <= 4 km above surface
  xpdmat.zPBLH_Ri_truecrossing_Rcrit(ii) = xpdmat.zPBLH_Ri(ii)/1000;
  xpdmat.pPBLH_Ri_truecrossing_Rcrit(ii) = xpdmat.pPBLH_Ri(ii);

  dTvpot_surf = xpdmat.Tvirtual_potential(NNlevs,ii) - xpdmat.stemp_pot(ii);
  bIsConvective = (dTvpot_surf < -5);   %% strongly unstable: surface >> air above

  xpdmat.PBLH_Ri_flag(ii) = -2;
   
  if (xpdmat.zPBLH_Ri(ii)/1000 - xpdmat.salti(ii)/1000 <= ZMAX)
    %% Ri crossing within 4km above surface - good result, keep as-is
    xpdmat.PBLH_Ri_flag(ii) = 0;

  else
    %% Ri crossing > 4km or NaN - convective/unstable profile or wind shear issue
    %% see jpss gran 48, fov 12150, 2024/11/13 over Australia land
    %%   Ri starts out negative, starts swinging towards zero but at 3.8 km there
    %%   is a lot of wind shear and windspeed starts dropping, so Ri swings back
    %%   negative again. It only becomes >= 0.25 at about 8 km.

    dz = xpdmat.zalts(1:NNlevs,ii) - xpdmat.salti(ii);
    RR = xpdmat.Ri(1:NNlevs,ii);
    boo = find(dz/1000 <= ZMAX);

    if bIsConvective
      %% Convective BL: Ri is negative throughout, no meaningful threshold crossing
      %% Best fallback: use virtual potential temperature gradient method
      %% which is physically meaningful for both stable and unstable BLs
      xpdmat.zPBLH_Ri(ii) = xpdmat.zPBLH_Tpot(ii);   %% already in metres at this point
      xpdmat.pPBLH_Ri(ii) = xpdmat.pPBLH_Tpot(ii);
      xpdmat.PBLH_Ri_flag(ii) = 2;   %% convective: Ri failed, used Tpot gradient

    else
      %% Non-convective but Ri crosses threshold above 4km (wind shear case)
      %% Find where Ri is closest to RiCritical within 0-4km above surface
      RRX = abs(RR(boo) - RiCritical);
      moo = find(RRX == min(RRX));
      moo = moo(end);   %% take surface-side (end = lowest layer in descending zalts)
      xpdmat.zPBLH_Ri(ii) = xpdmat.zalts(boo(moo),ii);
      xpdmat.pPBLH_Ri(ii) = interp1(levels_alts,levels_pres,xpdmat.zPBLH_Ri(ii),[],'extrap');
      xpdmat.PBLH_Ri_flag(ii) = 3;   %% wind shear case: capped at 4km, min|Ri-Ric|
    end
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %{
  if ii == 1 & iVers_Ri == 1
    axpdmat.Ri(1:NNlevs,ii) = (tp - tps) .* (xpdmat.zalts(1:NNlevs,ii) - xpdmat.salti(ii));  %% Davy only has (tp-tps)*(z)
    %axpdmat.Ri(1:NNlevs,ii) = g ./ tps ./s2 .* axpdmat.Ri(1:NNlevs,ii);
    axpdmat.Ri(1:NNlevs,ii) = g ./ tps ./(s2 + 100*100) .* axpdmat.Ri(1:NNlevs,ii);

    dz = diff(zalts);
    dz = dz(2:NNlevs);
    numer = diff(tp);
    denom = diff(xs2);
    xpdmat2.Ri(2:NNlevs,ii) = g ./ tp(2:NNlevs,ii) .* dz .* numer ./denom;
    
    plot(xpdmat2.Ri(1:NNlevs,ii),xpdmat.zalts(1:NNlevs,ii)/1000,'gx-',xpdmat.Ri(1:NNlevs,ii),xpdmat.zalts(1:NNlevs,ii)/1000,'bx-',axpdmat.Ri(1:NNlevs,ii),xpdmat.zalts(1:NNlevs,ii)/1000,'rx-'); axis([-40 120 0 10])

    plot(xpdmat.ptemp(1:NNlevs,ii),xpdmat.zalts(1:NNlevs,ii)/1000,'bx-',xpdmat.Tvirtual_potential(1:NNlevs,ii),xpdmat.zalts(1:NNlevs,ii)/1000,'bx-')
    keyboard_nowindow
  end
  %}
  
end

xpdmat.zPBLH_Ri   = xpdmat.zPBLH_Ri/1000;
xpdmat.salti      = xpdmat.salti/1000;

iPlot = -1;
if iPlot > 0
  plot_richardson_PBLH_levels
end

