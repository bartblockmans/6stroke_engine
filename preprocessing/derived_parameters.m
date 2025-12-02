function params = derived_parameters(params)

% Gear tooth numbers
params.zp = 2 * params.Rp / params.mn;
params.zr = 2 * params.Rr / params.mn;

% Center distance
params.a = params.Rr - params.Rp;

% Transmission ratio
params.u_pc = - (params.Rr - params.Rp) / params.Rp;