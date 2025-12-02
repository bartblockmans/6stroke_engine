function params = convert_to_SI(params)

% Convert geometrical parameters to SI units
params.a = params.a * 1e-3; % [mm] -> [m]
params.B = params.B * 1e-3; % [mm] -> [m]
params.Rr = params.Rr * 1e-3; % [mm] -> [m]
params.Rp = params.Rp * 1e-3; % [mm] -> [m]
params.delta = params.delta * 1e-3; % [mm] -> [m]
params.L = params.L * 1e-3; % [mm] -> [m]
params.Vc = params.Vc * 1e-9; % [mm³] -> [m³]
