function X = specific_humidity_to_vmr(q)
% Convert specific humidity to volume mixing ratio
% q   - specific humidity [kg/kg]
% X   - volume mixing ratio [mol/mol]

Mw  = 0.018015;
Md  = 0.028964;

X = (q .* Md) ./ (Mw + q .* (Md - Mw));

end
