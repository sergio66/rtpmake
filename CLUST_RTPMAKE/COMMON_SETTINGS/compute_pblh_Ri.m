function [pblh,Rib,junkstuff] = compute_pblh_Ri(z, T, p, w, u, v, skt, windvel10, landfrac, Rib_crit, b_us2, iUseCdrag)

% COMPUTE_PBLH  Estimate planetary boundary layer height via bulk Richardson number.
%
% INPUTS (all column vectors, ordered surface → model top):
%   z         - altitude [m]
%   T         - air temperature [K]
%   p         - pressure [Pa]
%   w         - water vapor mixing ratio [kg/kg]
%   u, v      - wind components [m/s]
%   skt       - surface temperature [K]  
%   windvel10 - wind speed at 10 m [m/s], has two components u,v
%   Rib_crit  - critical Richardson number (default: 0.25)
%   b_us2     - turbulent velocity scale term u*^2 [m^2/s^2] (default: 0.0)
%
% OUTPUT:
%   pblh     - planetary boundary layer height [m]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{

Rib = g/theta_vs (theta_v - theta_vs) (z - zs) / (uz − us)^2 + (vz−vs)^2+bu∗2

theta_v is virtual potential temperature
z_s is the surface height
u∗​,v∗​ are wind components
bu∗2​ is a turbulent velocity scale term (accounts for convective mixing)

The critical Ri_b threshold is typically 0.25

Key Characteristics of the IFS Scheme
--------------------------------------
Stable conditions: PBLH is found where Ri_b exceeds the critical value
going upward from the surface

Unstable/convective conditions: The scheme accounts for thermals and
convective overshooting via the bu∗2b u_*^2 bu∗2​ term

The TKE scheme: ERA5 uses a first-order turbulence closure with a
mass-flux scheme for convective boundary layers — the PBLH feeds into
the eddy-diffusivity/mass-flux (EDMF) framework

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 12
  %% see COMMON_SETTINGS/ECMWF_PBLH_ustar*.pdf    
  iUseCdrag = -1;
  %% see COMMON_SETTINGS/ECMWF_PBLH_ustar*.pdf  
end

set_Ri_critical  %% to get iVers_Ri

if nargin < 9, error('need z,P,T,q,u,v,skt,wspeed,landfrac'); end
if nargin < 10 || isempty(Rib_crit)
  Rib_crit = 0.25;
  Rib_crit = RiCritical;
end
if nargin < 11 || isempty(b_us2),    b_us2    = 0.0;  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

theta_v_excess = 0.0;

if length(windvel10) == 1
  if abs(windvel10) > 20
    windvel10 = 5;
  end
  windvel10(2) = 0.0;
end

[Y,I] = sort(z);  %% make sure index 1 === lowest altitude
z = z(I);
T = T(I);
p = p(I)*100;   %% mb to N/m2
w = w(I);
u = u(I);
v = v(I);

%% no need to sort wspeed, landfrac, skt

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% I have a ocean bias : if over ocean, PBLH is too low by about 0.25-0.5 km
if landfrac < 0.025    
  %% ocean -- could be ice though?????

  %% The ocean bias is a known issue with the basic bulk Ri
  %% method. Over the ocean, the surface layer is often well-mixed
  %% and near-neutral, so the wind shear denominator stays large and
  %% Ri builds up slowly, causing the method to miss the true
  %% inversion base. A few fixes:
  
  %% 1. Use a Slightly Higher Ri_crit Over Ocean
  %% Some implementations use Ri_crit = 0.3 over ocean vs 0.25 over land, which pushes PBLH up:
  %% tests showed d(Rcritical) of 0.02 = 0.06 km   
  Rib_crit_land  = 0.25;
  Rib_crit_ocean = 0.30;
  %%% Rib_crit_ocean = RiCritical;  %%% this should be correct?????
  Rib_crit = Rib_crit_ocean;
  
  %% 2. Add small virtual temperature excess to surface parcel  
  theta_v_excess = 0.5;    % K, typical range 0.3–1.0 K over ocean  
  %tvs = theta_v(1) + theta_v_excess;
  %% ERA5/IFS actually does something like this internally via the surface flux scheme.

  %% 3. Add a Turbulent Velocity Scale (most impactful)
  %% Over ocean, surface fluxes matter more. The bu∗2bu_*^2
  %% bu∗2​ term keeps Ri small in the mixed layer and pushes the
  %% diagnosed top upward:
  %% This alone often recovers 0.3–0.5 km over ocean — likely your main fix.

  %. Air-Sea Interactions (Oceanography & Meteorology)In meteorology
  % and ocean modeling, the drag coefficient parameterizes momentum
  % transfer at the air-sea interface. It is highly wind-speed dependent
  % and varies in distinct stages:
  %
  % Low Winds (< 5 m/s) : C_d often
  % decreases as wind speed increases over smooth water surfaces.Moderate
  %
  % Winds (5 - 30 m/s): \(C_{D}\) increases linearly or parabolically
  % with wind speed as ripples and larger gravity waves form, increasing
  % surface roughness.
  %
  % Extreme Winds / Hurricanes (> 30 m/s): Field and laboratory
  % observations show that \(C_{D}\) levels off (saturates) or even
  % decreases. This occurs because intense wave-breaking flattens the
  % sea surface and introduces thick layers of sea spray, essentially
  % causing "slip" and reducing aerodynamic friction.
  
  S10   = u(1)^2 + v(1)^2;  % u,v components at lowest level
  U10   = sqrt(S10);        % used this for surface windspeed over ocean

  U10   = norm(windvel10);
  %%%%%%%%%%%%%%%%%%%%%%%%%
  Cd    = 1.2e-3;           % typical ocean drag coefficient, constant

  % Cd depending on U10
  Cd(U10 < 11)  = 1.2e-3;   % low speed
  % Moderate to high winds
  Cd(U10 >= 11) = 1e-3 .* (0.49 + 0.065 .* U10(U10 >= 11));
  % Cap at very high winds
  Cd = min(Cd, 2.5e-3);
  %%%%%%%%%%%%%%%%%%%%%%%%%

  %%%%%%%%%%%%%%%%%%%%%%%%%
  u_st  = sqrt(Cd) * U10;   % friction velocity [m/s]
  b     = 100;              % Vogelezang & Holtslag 1996, constant
  b     = 100 + 5 .* U10;   % grows with wind speed  
  b_us2 = b * u_st^2;
  %%%%%%%%%%%%%%%%%%%%%%%%%

  % u(1)
  % v(1)
  % windvel10
  % U10
end

if iUseCdrag < 0
  b_us2 = 0.0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Virtual potential temperature ---
theta_v = compute_theta_v(T, p, w);

% --- Surface reference values ---
% --- Remember, everything sorted with ascending z ---
g  = 9.81;
zs = z(1);

%%%%%%%%%%%%%%%%%%%%%%%%%
%% tvs  = theta_v(1);
%% tvs  = theta_v(1) + 1.0*theta_v_excess;  %% 0.2 km too low for ECM
%% tvs  = theta_v(1) + 2.5*theta_v_excess;  %% 0.2 km too high for ECM
%% tvs  = theta_v(1) + 1.25*theta_v_excess; %% 0.2 km too high for ECM

%theta_skin   = skt .* (p0 ./ p(1)) .^ Rd_cp;
%tvs   = theta_skin .* (1 + 0.61 .* q(1,:)); 
tvs   = compute_theta_v(skt,p(1),w(1));
%%%%%%%%%%%%%%%%%%%%%%%%%

us   = u(1);
vs   = v(1);

us = windvel10(1);
vs = windvel10(2);

%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Bulk Richardson number at each level ---
nz  = numel(z);
Rib = zeros(nz, 1);

if iVers_Ri == 0
  %% b_us2 = 0;
  for k = 2:nz
      dtheta = theta_v(k) - tvs;
      dz     = z(k) - zs;
      dwind2 = (u(k) - us)^2 + (v(k) - vs)^2 + b_us2;
      junk(k) = dwind2;
      dwind2 = max(dwind2, 1e-6);           % avoid division by zero
      Rib(k) = (g / tvs) * dtheta * dz / dwind2;
  end
  junkstuff.tps       = tvs;
  junkstuff.tp        = theta_v;  

elseif iVers_Ri == 1
  % Tv_surf  = T_surf * (1 + 0.61*q_surf);
  % DSEv_surf = cp*Tv_surf + g*z_surf;

  cp = 1004;      % J/(kg K)
  g  = 9.81;      % m/s^2
  
  Tv_surf  = skt * (1 + 0.61*w(1));
  DSEv_surf = cp*Tv_surf + g*zs;

  DSEv = compute_dse(T,w,z);

  junk = (u-us).^2 + (v-vs).^2;
  WS2 = max(junk, 0.01);
  Rib = (DSEv - DSEv_surf) ./ (cp * Tv_surf) .* (z - zs) ./ WS2;
  Rib = Rib * g;
  
  junk = junk';
  junkstuff.tps       = Tv_surf;
  junkstuff.tp        = DSEv;
end

junkstuff.speed_sqr = junk';
junkstuff.dz        = z-zs;
junkstuff.b_us2     = b_us2;
junkstuff.u         = u;
junkstuff.v         = v;

% --- Find first level where Rib >= Rib_crit ---
idx = find(Rib >= Rib_crit, 1, 'first');

%% plot(Rib,z/1000); ylim([0 4]); title('compute\_pblh\_Ri'); 

if isempty(idx)
  pblh = z(end);                        % PBL top above profile
  return
end
if idx == 1
  pblh = z(1);
  return
end

% --- Linear interpolation between bracketing levels ---
k0   = idx - 1;
k1   = idx;
frac = (Rib_crit - Rib(k0)) / (Rib(k1) - Rib(k0));
pblh = z(k0) + frac * (z(k1) - z(k0));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function theta_v = compute_theta_v(T, p, w)
% COMPUTE_THETA_V  Virtual potential temperature.
%
%   T  - temperature [K]
%   p  - pressure [Pa]
%   w  - mixing ratio [kg/kg]

p0    = 100000.0;                         % reference pressure [Pa]
Rd_cp = 287.04 / 1004.0;                 % R_d / c_p

q       = w ./ (1 + w);                  % mixing ratio → specific humidity
theta   = T .* (p0 ./ p) .^ Rd_cp;      % potential temperature
theta_v = theta .* (1 + 0.61 .* q);     % virtual potential temperature

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function DSE_v = compute_dse(T,q,z)
% COMPUTE dry static energy

% DSE_v = cp * Tv + g * z
% Where:
% cp = specific heat of dry air at constant pressure = 1004 J/(kg·K)
% Tv = virtual temperature (K) = T * (1 + 0.61*q)
% g = gravitational acceleration = 9.81 m/s²
% z = geopotential height (m)

cp = 1004;      % J/(kg K)
g  = 9.81;      % m/s^2

Tv = T .* (1 + 0.61*q);     % virtual temperature, T in K, q in kg/kg

DSE_v = cp * Tv + g * z;    % J/kg

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{

cp = 1004;
g  = 9.81;

Tv_surf  = T_surf * (1 + 0.61*q_surf);
DSEv_surf = cp*Tv_surf + g*z_surf;

Tv  = T .* (1 + 0.61*q);
DSEv = cp*Tv + g*z;

WS2 = max((u-u_surf).^2 + (v-v_surf).^2, 0.01);

Ri = (DSEv - DSEv_surf) ./ (cp * Tv_surf) .* (z - z_surf) ./ WS2;

Why use DSE_v instead of Tv_pot? A few reasons some groups prefer it:

It's more directly tied to the energy budget of the PBL (since DSE is
what's conserved in dry adiabatic processes, with g*z baked in
explicitly rather than implicitly via potential temperature's
reference pressure)
Avoids any ambiguity in choice of reference pressure (1000 hPa) used
in the potential temperature definition
Numerically can behave better at very low pressures /very high
altitudes, though for PBL work (lowest few km) this rarely matters

For your specific application — since you're already computing Tv_pot
and it's working reasonably well (the 0.2km/0.5km biases we
discussed), there's likely no compelling reason to switch to DSE_v
unless you're matching a specific paper's methodology or ERA5
documentation that uses it.


%%%%%%%%%%%%%%%%%%%%%%%%%

Can you check if there is a missing value of g!!!!???

Start from
  DSEv = cp*Tv + g*z
So
  DSEv(z) - DSEv_surf = cp*(Tv(z) - Tv_surf) + g*(z - z_surf)
Thus
  Ri = [DSEv(z) - DSEv_surf] / (cp*Tv_surf) * (z-z_surf) / WS^2
     = [cp*(Tv(z)-Tv_surf) + g*(z-z_surf)] / (cp*Tv_surf) * (z-z_surf) / WS^2

That's wrong — it has an extra (z-z_surf) factor multiplying a term
that already contains g*(z-z_surf), giving a (z-z_surf)² term that
shouldn't be there, plus it doesn't reduce to the standard Ri
formula.

Let me redo this properly. The standard Ri is:
  Ri = (g/Tv_surf) * (Tv(z) - Tv_surf) * (z-z_surf) / WS^2

The only correct way to write this using DSE_v is to recognize that
g*(Tv(z)-Tv_surf)*(z-z_surf) is NOT simply related to
DSEv(z)-DSEv_surf in a clean way — DSE_v mixes the temperature and
height terms together additively, so you can't cleanly factor out
just the buoyancy part from it without separating cp*Tv and g*z back
out again.

So my previous formula was dimensionally and physically wrong — I
incorrectly conflated the two. The correct approach is simpler: just
stick with the standard Tv_pot formula directly:

Ri = (g/Tv_surf) * (Tv_pot(z) - Tv_pot_surf) * (z-z_surf) / WS^2

DSE_v isn't a drop-in replacement for Tv_pot in the Ri formula — they're different quantities serving different diagnostic purposes. Thanks for pushing back, that formula needed correcting rather than just a missing g.

%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
ERA documnetation on PBL : wget https://www.ecmwf.int/sites/default/files/elibrary/2017/17736-part-iv-physical-processes.pdf

This is the IFS bulk Richardson number formula using virtual dry
%static energy — exactly the quantity we discussed earlier! ECMWF uses
%this rather than virtual potential temperature directly.

The equations are:

|ΔU|² = u²_hbl + v²_hbl                                    (wind shear, ECMWF sets near-surface wind to 0)

s_vn   = cp·T_n·(1 + ε·q_n) + g·z_n                        (virtual DSE at lowest model level n)

s_v,hbl = cp·T_hbl·(1 + ε·q_hbl) + g·h_bl                  (virtual DSE at boundary layer height)  -- eqn 3.90

Ri_b = h_bl · 2g·(s_v,hbl - s_vn) / [(s_v,hbl + s_vn - g·h_bl - g·z_n)·|ΔU|²]

Key clarifications from the text:

eps here is the virtual temperature coefficient (≈0.61), and T·(1+eps·q)
is exactly the virtual temperature Tv we've been using — so s_v =
cp·Tv + g·z, confirming the virtual dry static energy formula from
earlier.

No friction velocity term! Since friction velocity is not known from
%radiosonde data, surface frictional effects are ignored in the bulk
%shear computation, and winds at 2m are set to zero — confirming what
%your code comment already said, and what we found in the Seidel et
%al. search earlier. You're right not to have ustar — ECMWF doesn't
%use it either in this diagnostic.
Arch Linux Man Pages

This is the Vogelezang & Holtslag (1996) method, found by Seidel et
al. (2012) to be the most appropriate algorithm for radiosondes,
reanalysis, and climate model data, suitable for both convective and
stable boundary layers, identifying a nonnegative height in all
cases. Arch Linux Man Pages

Scan direction: the boundary layer height is found by a vertical scan
from the surface upwards, with linear interpolation if the crossing
falls between two levels.

Important subtlety in the denominator: it's not just cp·Tv_surf like I
wrote earlier — it's 
  (s_v,hbl + s_vn - g·h_bl - g·z_n)
which expands to 
  cp·(Tv_hbl + Tv_n) — 
i.e. it uses the average of Tv at the two
levels, not just the surface value. That's a meaningful difference
from the simpler formula and from what I gave you earlier.

Key things to flag about this implementation:

1. Confirms what you already had right — ECMWF ignores friction
velocity (ustar) entirely in this diagnostic, and the original
Vogelezang & Holtslag near-surface wind treatment for radiosondes
sets the 2m wind to zero — but note for model-level application
(which is what you're doing with full profiles, not radiosondes),
you'd typically use the actual lowest-level wind u_n, v_n, not
zero. I left this as the actual lowest-level wind with a comment
explaining the choice — let me know if you want the strict
radiosonde-zero version instead. Arch Linux Man Pages

2. The big difference from a simple Tv_pot Ri — the denominator uses
(s_v,hbl + s_vn), the sum/average of virtual DSE at both levels, not
just the surface value. This makes the formula self-referential —
h_bl appears on both sides — which is why the function scans
candidate levels rather than computing Ri directly level-by-level
like your original code.

3. Scan direction — the IFS scans from the surface upward and finds
the lowest level where Ri_b reaches the critical value of 0.25, then
linearly interpolates between levels for the exact height. I
implemented exactly that. Arch Linux Man Pages

4. Should fix your convective ocean PBLH=0 problem — Seidel et
al. (2012) showed this algorithm identifies a nonnegative height in
all cases and works for both convective and stable boundary layers —
the wind shear term |ΔU|² between levels (not just surface) is what
gives it this property, unlike the simple surface-referenced Ri. Arch
Linux Man Pages

One important note — your data is TOA-first (descending), so I added a
%sort step at the top of the example usage to flip to surface-up
%before calling the function.

All this is in 
  get_PBLH_VogelezangHoltslag.m
  
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
