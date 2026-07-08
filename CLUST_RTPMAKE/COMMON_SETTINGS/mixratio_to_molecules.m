function N_mol_cm2 = mixratio_to_molecules(q_mmr, dp, gid, g)

% this is opposite of molecules_to_mixratio.m
% Convert mixing ratio [kg/kg] to layer water vapor [molecules/cm^2]
%
% INPUTS:
%   q_mmr      - water vapor mixing ratio [kg/kg]
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

if gid == 1
  q_specific_hum = q_mmr ./ (1 + q_mmr);
  N_mol_m2  = (q_specific_hum * Na * dp) ./ (Mw * g);
  N_mol_cm2 = N_mol_m2 / 1e4;  
else
  q_specific_hum = q_mmr * nan;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
