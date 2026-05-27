function pblh = compute_pblh(z, T, p, w, u, v, skt, wspeed10, landfrac, Rib_crit, b_us2)

% COMPUTE_PBLH  Estimate planetary boundary layer height via bulk Richardson number.
%
% INPUTS (all column vectors, ordered surface → model top):
%   z        - altitude [m]
%   T        - air temperature [K]
%   p        - pressure [Pa]
%   w        - water vapor mixing ratio [kg/kg]
%   u, v     - wind components [m/s]
%   skt      - surface temperature [K]  
%   wspeed10 - wind speed at 10 m [m/s], has two components u,v
%   Rib_crit - critical Richardson number (default: 0.25)
%   b_us2    - turbulent velocity scale term u*^2 [m^2/s^2] (default: 0.0)
%
% OUTPUT:
%   pblh     - planetary boundary layer height [m]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{

Rib​= g/theta_vs (theta_v - theta_vs) (z - zs) / (uz − us)^2 + (vz−vs)^2+bu∗2

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

if nargin < 9, error('need z,P,T,q,u,v,skt,wspeed,landfrac'); end
if nargin < 10 || isempty(Rib_crit), Rib_crit = 0.25; end
if nargin < 11 || isempty(b_us2),    b_us2    = 0.0;  end

theta_v_excess = 0.0;

if length(wspeed10) == 1
  wspeed10(2) = 0.0;
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
  
  S10   = u(1)^2 + v(1)^2;  % u,v compomemnts at lowest level
  U10   = sqrt(S10);        % used this for surface windspeed over ocean

  U10   = norm(wspeed10);
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
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Virtual potential temperature ---
theta_v = compute_theta_v(T, p, w);

% --- Surface reference values ---
g    = 9.81;
zs   = z(1);

%%%%%%%%%%%%%%%%%%%%%%%%%
%% tvs  = theta_v(1);
%% tvs  = theta_v(1) + 1.0*theta_v_excess;  %% 0.2 km too low for ECM
%% tvs  = theta_v(1) + 2.5*theta_v_excess;  %% 0.2 km too high for ECM
%% tvs  = theta_v(1) + 1.25*theta_v_excess;  %% 0.2 km too high for ECM

%theta_skin   = skt .* (p0 ./ p(1)) .^ Rd_cp;
%tvs   = theta_skin .* (1 + 0.61 .* q(1,:)); 
tvs   = compute_theta_v(skt,p(1),w(1));
%%%%%%%%%%%%%%%%%%%%%%%%%

us   = u(1);
vs   = v(1);

us = wspeed10(1);
vs = wspeed10(2);

%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Bulk Richardson number at each level ---
nz  = numel(z);
Rib = zeros(nz, 1);

for k = 2:nz
    dtheta = theta_v(k) - tvs;
    dz     = z(k) - zs;
    dwind2 = (u(k) - us)^2 + (v(k) - vs)^2 + b_us2;
    dwind2 = max(dwind2, 1e-6);           % avoid division by zero
    Rib(k) = (g / tvs) * dtheta * dz / dwind2;
end

% --- Find first level where Rib >= Rib_crit ---
idx = find(Rib >= Rib_crit, 1, 'first');

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
