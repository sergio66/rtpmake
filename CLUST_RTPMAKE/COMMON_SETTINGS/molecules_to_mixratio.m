function [w_mmr,q_specific_hum] = molecules_to_mixratio(N_mol_cm2, dp, nlays, gid, g)

% Convert layer-averaged water vapor [molecules/cm2] back to mixing ratio (sp humidity)
% 
%
% INPUTS:
%   N_mol_cm2  - water vapor column density [molecules/cm^2]
%   dp         - layer pressure thickness [Pa]  (positive)
%   gid        - HITRAN gasID  
%   g          - gravity [m/s^2], default 9.80665
%
% OUTPUT:
%   w_mmr          - water vapor mixing ratio [kg/kg]
%   q_specific_hum - specific humidity [kg/kg]

%{
see ../COMMON_SETTINGS/testing_layers2gg_layers2sphum_conversions.m
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
  disp('[w_mmr,q_specific_hum] = molecules_to_mixratio(N_mol_cm2, dp, nlays, gid, g)');    
  error('need at least 3 argument : profileQ, profile dp, numlay per profiles   .... gid and gravity are optional')
elseif nargin == 3
  gid = 1;
  g = 9.80665;
elseif nargin == 4
  g = 9.80665;
end

if size(N_mol_cm2) ~= size(dp)
  error('need size(N_mol_cm2N_mol_cm2) == size(dp)')
end  
[mm,nn] = size(N_mol_cm2);
if nn ~= length(nlays)
  error('need [mm,nn] size(N_mol_cm2N_mol_cm2); nn == length(nlays)')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if gid == 1
  Mw    = 0.018015;        % molar mass water [kg/mol]
  mass_g = 18;
  mass_g = Mw * 1000;
elseif gid == 2
  mass_g = 44;
elseif gid == 3
  mass_g = 48;
elseif gid == 4
  mass_g = 44;
elseif gid == 5
  mass_g = 28;
elseif gid == 6
  mass_g = 16;
elseif gid == 9
  mass_g = 64;
elseif gid == 12
  mass_g = 63;
else
  gid
  error('can only handle gid 1 2 3 4 5 6 9 12 since only have mass_g of these gases ... need to (easily) add mass of this new gas')
end

Na = 6.02214076e23;   % Avogadro [/mol]
%Mw    = 0.018015;
Mw = mass_g/1000;

% molecules/cm^2 → molecules/m^2
N_mol_m2 = N_mol_cm2 * 1e4;

% specific humidity [kg/kg]
q_specific_hum = (N_mol_m2 * Mw * g) ./ (dp * Na);

% mixing ratio [kg/kg]
w_mmr = q_specific_hum ./ (1 - q_specific_hum);

for ii = 1 : nn
  q_specific_hum(nlays(ii)+1:mm,ii) = nan;
  w_mmr(nlays(ii)+1:mm,ii)          = nan;  
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function N_mol_cm2 = mixratio_to_molecules(w_mmr, dp, gid, g)
% Convert mixing ratio [kg/kg] to layer water vapor [molecules/cm^2]
%
% INPUTS:
%   w_mmr      - water vapor mixing ratio [kg/kg]
%   dp         - layer pressure thickness [Pa]
%   g          - gravity [m/s^2], default 9.80665
%
% OUTPUT:
%   N_mol_cm2  - water vapor column density [molecules/cm^2]

if nargin == 2
  gid = 1;
  g = 9.80665;
elseif nargin == 3
  g = 9.80665;  
end

if gid == 1
  Mw    = 0.018015;        % molar mass water [kg/mol]
  mass_g = 18;
  mass_g = Mw * 1000;
elseif gid == 2
  mass_g = 44;
elseif gid == 3
  mass_g = 48;
elseif gid == 4
  mass_g = 44;
elseif gid == 5
  mass_g = 28;
elseif gid == 6
  mass_g = 16;
elseif gid == 9
  mass_g = 64;
elseif gid == 12
  mass_g = 63;
else
  gid
  error('can only handle gid 1 2 3 4 5 6 9 12 since only have mass_g of these gases ... need to (easily) add mass of this new gas')
end

Na    = 6.02214076e23;
%Mw    = 0.018015;
Mw = mass_g/1000;

q_specific_hum = w_mmr ./ (1 + w_mmr);

N_mol_m2  = (q_specific_hum * Na * dp) ./ (Mw * g);
N_mol_cm2 = N_mol_m2 / 1e4;

end
