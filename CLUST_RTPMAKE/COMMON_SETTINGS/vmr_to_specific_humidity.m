function q = vmr_to_specific_humidity(X,iDryOrWet)

% Convert volume mixing ratio to specific humidity
% X   - volume mixing ratio [mol/mol, not ppmv]
% q   - specific humidity [kg/kg]
% iDryOrWet = +1,+2 if using "dry" or accounting for "moist"
  
Mw  = 0.018015;    % molar mass water [kg/mol]
Md  = 0.028964;    % molar mass dry air [kg/mol]

% Pavg is moist air pressure — it includes the water vapor partial
% pressure e. So when you compute X = NH2O/NtotalX = N_{H_2O} / N_{total}
% the denominator already includes water vapor molecules. That's correct for VMR.
% 
% But when converting VMR -> specific humidity, if you use M_{dry}
% for the mean molar mmass of the denominator gas, you're slightly wrong over warm
% moist ocean where water vapor is 2-4% of total pressure:

if nargin == 1
  iDryOrWet = 1;
end

if iDryOrWet == 1
  M_moist = X .* Mw + (1 - X) .* Md;
  q = (X .* Mw) ./ M_moist;
elseif iDryOrWet == 2
  q = (X .* Mw) ./ (Md + X .* (Mw - Md));
end

end

