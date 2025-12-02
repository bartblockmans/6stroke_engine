function params = ValveLiftForAnimation(params)

% Extra angle
deg_extra = 30;

% Phi_deg
phi_deg = linspace(0, 1080, 2161);

% Valve/port timing areas for reference/plots
A_int_m2  = valve_area_1080(phi_deg, params.IVO_A - deg_extra, params.IVC_A + deg_extra, params.D_int,  params.Lift_int);
A_exh_m2  = valve_area_1080(phi_deg, params.EVO_C1 - deg_extra, params.EVC_C1 + deg_extra, params.D_exh, params.Lift_exh) ...
          + valve_area_1080(phi_deg, params.EVO_F1 - deg_extra, params.EVC_F1 + deg_extra, params.D_exh, params.Lift_exh);

% Effective lifts from areas
Lift_int_m = zeros(size(phi_deg));
Lift_exh_m = zeros(size(phi_deg));
if params.D_int > 0
    Lift_int_m = A_int_m2 ./ (pi * params.D_int);
end
if params.D_exh > 0
    Lift_exh_m = A_exh_m2 ./ (pi * params.D_exh);
end

% Store
params.L_int_m = Lift_int_m;
params.L_exh_m = Lift_exh_m;

function A = valve_area_1080(phi_deg, open_deg, close_deg, Dv, Lmax)
% half-cosine lift between open/close across 1080 wrap
ph = mod(phi_deg,1080);
o  = mod(open_deg,1080);
c  = mod(close_deg,1080);

if o < c
    active = (ph >= o) & (ph <= c);
    dph    = c - o;
    tau    = (ph - o)/dph;
else
    active = (ph >= o) | (ph <= c);
    dph    = mod(c - o, 1080);
    tau    = zeros(size(ph));
    m1     = active & (ph >= o);
    m2     = active & (ph <= c);
    tau(m1)= (ph(m1) - o)/dph;
    tau(m2)= (ph(m2) + (1080 - o))/dph;
end

L = zeros(size(ph));
L(active) = Lmax * 0.5 .* (1 - cos(2*pi*tau(active)));
A = pi * Dv .* L;
end

end