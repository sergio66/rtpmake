function [w,q] = recover_q_from_forward(Q_mol_cm2, P_avg, T_avg, L, nlays, gid)

% Recover specific humidity from your forward integration
%
% in
%   Q_mol_cm2  - water vapor column [molecules/cm^2]
%   P_avg      - layer mean pressure [Pa]
%   T_avg      - layer mean temperature [K]
%   L          - layer geometric thickness [m]
%   nlays      - number of layers in each profile
%
% out
%  w,q           - mass mix ratio and specific humidity

if nargin == 5
  gid  = 1;
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

R   = 8.31446;     % universal gas constant [J/mol/K]
Na  = 6.02214076e23;

% Total molecules/cm^2 in layer
N_total_cm2 = (P_avg .* L .* Na) ./ (R .* T_avg .* 1e4);

% Volume mixing ratio
X = Q_mol_cm2 ./ N_total_cm2;

% Convert VMR to specific humidity
q = vmr_to_specific_humidity(X);
if gid ~= 1
  q = q * mass_g/18;  %% vmr_to_specific_humidity assumes WV, mass 18
end  

w = q./(1-q);

[mm,nn] = size(P_avg);
for ii = 1 : nn
  q(nlays(ii)+1:mm,ii) = nan;
  w(nlays(ii)+1:mm,ii) = nan;  
end

end
