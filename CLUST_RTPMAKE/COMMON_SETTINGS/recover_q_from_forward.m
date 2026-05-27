function [w,q] = recover_q_from_forward(Q_mol_cm2, P_avg, T_avg, L, nlays)

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

R   = 8.31446;     % universal gas constant [J/mol/K]
Na  = 6.02214076e23;

% Total molecules/cm^2 in layer
N_total_cm2 = (P_avg .* L .* Na) ./ (R .* T_avg .* 1e4);

% Volume mixing ratio
X = Q_mol_cm2 ./ N_total_cm2;

% Convert VMR to specific humidity
q = vmr_to_specific_humidity(X);

w = q./(1-q);

[mm,nn] = size(P_avg);
for ii = 1 : nn
  q(nlays(ii)+1:mm,ii) = nan;
  w(nlays(ii)+1:mm,ii) = nan;  
end

end
