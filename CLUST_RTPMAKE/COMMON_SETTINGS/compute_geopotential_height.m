function z = compute_geopotential_height(T, w, sp, a_coef, b_coef, z_surf)

% Compute geopotential height on full model levels
% a_coef, b_coef: ECMWF hybrid coefficients, length nlev+1 (half levels)
% sp:             surface pressure [Pa]
% z_surf:         surface geopotential height [m] (orography / g)
% T, w:           [nlev x 1], surface = last index (ECMWF top-down ordering)

Rd  = 287.04;
g   = 9.80665;
nlev = length(T);

% --- Specific humidity ---
q = w ./ (1 + w);

% --- Virtual temperature ---
Tv = T .* (1 + 0.609 .* q);

% --- Half-level pressures (nlev+1 levels) ---
ph = a_coef + b_coef .* sp;        % [nlev+1 x 1]

% --- Full-level pressures ---
pf = 0.5 .* (ph(1:end-1) + ph(2:end));

% --- Thickness of each layer in ln(p) ---
alpha = zeros(nlev, 1);
for k = 1:nlev
    if k == 1 && ph(1) < 1e-9      % top of atmosphere
        alpha(k) = log(2);
    else
        alpha(k) = 1 - (ph(k) / (ph(k+1) - ph(k))) * log(ph(k+1)/ph(k));
    end
end

% --- Integrate geopotential upward from surface ---
% ECMWF integrates top-down, so we go from bottom (nlev) upward
phi = zeros(nlev, 1);
phi_half_below = z_surf * g;       % surface geopotential [m^2/s^2]

for k = nlev:-1:1
    dlog_p   = log(ph(k+1) / ph(k));
    phi_half_above = phi_half_below + Rd * Tv(k) * dlog_p;
    phi(k)   = phi_half_above + Rd * Tv(k) * alpha(k);
    phi_half_below = phi_half_above;
end

z = phi ./ g;   % geopotential height [m]
end
