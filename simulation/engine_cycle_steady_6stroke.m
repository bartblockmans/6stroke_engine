function result = engine_cycle_steady_6stroke(params)
%ENGINE_CYCLE_STEADY_6STROKE
% Integrate 6-stroke single-cylinder ICE model over multiple 1080° cycles
% until periodic steady state. State vector x = [p_cyl; m_cyl; mO2].
% -------------------------------------------------------------------------
% Bart Blockmans, 2024 - bart@blockmans.net
% -------------------------------------------------------------------------

N_rpm   = params.N_rpm;

% Sanity-check that all required geometric parameters are present (top-level)
geom_required = ["B","Rr","Rp","delta","L","Vc","psi0"];
ensure_fields(params, geom_required, 'engine_cycle_steady_6stroke geometry');

% Physical floors (can be overridden from params)
if ~isfield(params,'p_floor') || isempty(params.p_floor)
    params.p_floor = 5e3;         % [Pa] ~0.05 bar
end
if ~isfield(params,'m_floor') || isempty(params.m_floor)
    params.m_floor = 1e-9;        % [kg]
end
if ~isfield(params,'T_floor') || isempty(params.T_floor)
    params.T_floor = 200;         % [K]
end

%% Fuel per 1080° cycle and split
% 6-stroke cycle = 1080° = 3 crank revolutions.
% rev/s = N/60 => cycles/s = (N/60)/3 = N/180.
cycle_freq = N_rpm/180;           % [cycles/s]
period     = 1/cycle_freq;        % [s/cycle]
m_fuel_tot = params.mdot_fuel * period;     % [kg/1080°-cycle]
mf1 = params.fuel_split * m_fuel_tot;       % burn #1 fuel mass/cycle
mf2_cmd = (1 - params.fuel_split) * m_fuel_tot; % commanded burn #2 fuel/cycle

params.m_fuel_cycle_total = m_fuel_tot;
params.m_fuel1 = mf1;
params.m_fuel2_cmd = mf2_cmd;

%% Precompute geometry tables once (use global highest TDC for x_piston)
phi_tab = linspace(0,1080,2161);   % 0.5 deg grid
[V_mm3_tab, A_mm2_tab, ~, dVdth_mm3_per_deg_tab, ~, dxdth_mm_per_deg_tab, ~] = geometry_6stroke( ...
    phi_tab, params);

% Convert geometry outputs from mm-units to SI for the thermodynamics model
V_tab_m3      = V_mm3_tab * 1e-9;            % m^3
A_tab_m2      = A_mm2_tab * 1e-6;            % m^2
dVdphi_tab_m3 = dVdth_mm3_per_deg_tab * 1e-9; % m^3/deg
dxdphi_tab_m  = dxdth_mm_per_deg_tab * 1e-3;   % m/deg

% Cache geometry tables + monotone interpolants for reuse during ODE solve
params.geotab.phi     = phi_tab(:);
params.geotab.V       = V_tab_m3(:);
params.geotab.A       = A_tab_m2(:);
params.geotab.dVdphi  = dVdphi_tab_m3(:);
params.geotab.dxdphi  = dxdphi_tab_m(:);

params.geotab.V_s      = griddedInterpolant(params.geotab.phi, params.geotab.V,      'pchip');
params.geotab.A_s      = griddedInterpolant(params.geotab.phi, params.geotab.A,      'pchip');
params.geotab.dVdphi_s = griddedInterpolant(params.geotab.phi, params.geotab.dVdphi, 'pchip');
params.geotab.dxdphi_s = griddedInterpolant(params.geotab.phi, params.geotab.dxdphi, 'pchip');

% Compute reference displacement volume from geometry (max V - min V over cycle)
if ~isfield(params,'Vd_ref') || isempty(params.Vd_ref)
    params.Vd_ref = max(params.geotab.V) - min(params.geotab.V);  % [m^3]
end

% Initialize burn-2 oxygen scaling
if ~isfield(params,'lambda2') || isempty(params.lambda2)
    params.lambda2 = 1.0;
end
if ~isfield(params,'lambda2_relax') || isempty(params.lambda2_relax)
    params.lambda2_relax = 1.0;
end
params.lambda2_relax = max(0,min(1,params.lambda2_relax));  % clamp to [0,1]

%% Initial guess at phi = 0 (TDC1); near intake conditions
phi0 = params.phiSpan(1);
% Geometry from geometry_6stroke (converted to SI)
[V0, ~, ~] = geom_eval(phi0, params);
p0   = params.p_int;
T0   = params.T_int;
m0   = p0*V0/(params.R_mix*T0);
mO20 = m0 * params.YO2_air;   % initial O2 as if fresh air

x0 = [p0; m0; mO20];

%% Cycle loop with convergence on (p, m, mO2)
x_start = x0;
theta_sol_last = [];
x_sol_last     = [];
for cyc = 1:params.maxCycles
    % Integrate a single 1080° window
    rhs = @(ph, x) engine_rhs_6stroke(ph, x, params);
    opts = odeset('RelTol',1e-6,'AbsTol',1e-9);
    [phi_sol, x_sol] = ode15s(rhs, params.phiSpan, x_start, opts);

    % Hard-limit the solution to keep it in the physically admissible space
    x_sol = enforce_state_limits(x_sol, params);
    x_end = enforce_state_limits(x_sol(end,:).', params);

    % Update burn-2 oxygen cap (lambda2) using SOC2 snapshot
    try
        mO2_soc2 = interp1(phi_sol, x_sol(:,3), params.soc2, 'linear', 'extrap');
        m_soc2   = interp1(phi_sol, x_sol(:,2), params.soc2, 'linear', 'extrap');
        m_soc2   = max(m_soc2, params.m_floor);
        mO2_soc2 = min(max(mO2_soc2, 0), m_soc2);
        lambda2_new = mO2_soc2 / (params.nu_O2 * params.m_fuel2_cmd + 1e-12);
    catch
        lambda2_new = 1.0;
    end
    lambda2_new   = max(0.0, min(1.0, lambda2_new));
    params.lambda2 = params.lambda2 + params.lambda2_relax * (lambda2_new - params.lambda2);
    params.lambda2 = max(0.0, min(1.0, params.lambda2));

    % Scaled residual: p in bar, m in 0.1 g, mO2 in 0.01 g
    scale = [1e5; 1e-4; 1e-5];
    err = norm((x_end - x_start)./scale);
    fprintf('      Cycle %d: end-start scaled error = %.3e\n', cyc, err);

    if err < params.cycleTol
        fprintf('      Converged after %d cycles.\n', cyc);
        theta_sol_last = phi_sol;
        x_sol_last     = x_sol;
        break;
    end

    % Otherwise march to the next cycle with updated initial state
    x_start        = x_end;
    theta_sol_last = phi_sol;
    x_sol_last     = x_sol;
end

if isempty(theta_sol_last)
    warning('      Did not converge within maxCycles. Returning last cycle.');
    theta_sol_last = phi_sol;
    x_sol_last     = x_sol;
end

%% Interpolate onto uniform phi grid
phi_deg = params.phiGrid(:);
p_cyl   = interp1(theta_sol_last, x_sol_last(:,1), phi_deg, 'pchip');
m_cyl   = interp1(theta_sol_last, x_sol_last(:,2), phi_deg, 'pchip');
mO2     = interp1(theta_sol_last, x_sol_last(:,3), phi_deg, 'pchip');
state_interp = enforce_state_limits([p_cyl, m_cyl, mO2], params);
p_cyl = state_interp(:,1);
m_cyl = state_interp(:,2);
mO2   = state_interp(:,3);

% Temperature via EOS
[V_vec, ~, dV_dphi] = geom_eval(phi_deg, params);
T_cyl = p_cyl .* V_vec ./ (m_cyl * params.R_mix);
T_cyl = max(T_cyl, params.T_floor);

% Gas torque (two equivalent definitions)
T_ind       = gas_torque(phi_deg, p_cyl, params);        % via F * dx/dθ
T_ind_pdv   = gas_torque_pdv(phi_deg, p_cyl, params);    % via (p-p_crank) * dV/dθ

% IMEP over 1080° (use displacement volume from geometry)
p_dV    = p_cyl .* dV_dphi;
W_cycle = trapz(phi_deg, p_dV);                  % [J] = ∮ p dV over 1080°
IMEP_bar = (W_cycle / params.Vd_ref) / 1e5;      % [bar] vs geometry-derived displacement

% Mean indicated torque for 1080° = 6π rad
T_ind_mean = W_cycle / (6*pi);                   % [Nm]

% O2 mass at SOC2 (for info)
mO2_at_SOC2 = interp1(phi_deg, mO2, params.soc2, 'linear', 'extrap');

% Additional arrays for plotting and diagnostics
% ------------------------------------------------------------
% Geometry-derived piston kinematics and cumulative work
B_m     = params.B * 1e-3;                       % [m]
Ap_m2   = pi * (B_m^2) / 4;                      % [m^2]
V_min   = min(V_vec);
x_m     = max((V_vec - V_min) ./ Ap_m2, 0);      % [m] from highest TDC
x_mm    = x_m * 1e3;                            % [mm]
dx_dphi = params.geotab.dxdphi_s(mod(phi_deg,1080));  % [m/deg]
ddx_dphi2 = gradient(dx_dphi, mean(diff(phi_deg)));   % [m/deg^2]
W_cum   = cumtrapz(phi_deg, p_cyl .* dV_dphi);   % [J]

% Gas-side piston force
F_gas = (p_cyl - params.p_crank) * Ap_m2;        % [N]

% Heat terms
Aw_vec       = params.geotab.A_s(mod(phi_deg,1080));   % [m^2]
Qdot_wall_W  = Qdot_wall(phi_deg, p_cyl, T_cyl, Aw_vec, params);  % [W]
% Heat release per degree from the two Wiebes
Q1        = params.eta_comb1 * params.m_fuel1      * params.LHV;
Q2_cmd    = params.eta_comb2 * params.m_fuel2_cmd  * params.LHV;
dx1_dphi  = wiebe_dxb_dphi(phi_deg, params.soc1, params.dur1, params.m1, params.a1);
dx2_shape = wiebe_dxb_dphi(phi_deg, params.soc2, params.dur2, params.m2, params.a2);
dQ_dphi   = Q1 * dx1_dphi + (params.lambda2 * Q2_cmd) * dx2_shape;  % [J/deg]

% Burn fraction integrals (0→1 per Wiebe)
xburn1 = cumtrapz(phi_deg, dx1_dphi);
xburn2 = cumtrapz(phi_deg, dx2_shape);
if ~isempty(xburn1)
    xburn1 = xburn1 / max(xburn1(end), 1e-12);
end
if ~isempty(xburn2)
    xburn2 = xburn2 / max(xburn2(end), 1e-12);
end

% Time-domain utilities
deg_per_sec = 6 * params.N_rpm;
dt_per_deg  = 1 / deg_per_sec;                   % [s/deg]
t_sec       = (phi_deg - phi_deg(1)) * dt_per_deg;
omega_rad_s = params.N_rpm * 2*pi/60;
power_W     = T_ind * omega_rad_s;
imep_inst   = (p_cyl - params.p_crank) .* dV_dphi / max(params.Vd_ref, 1e-12); % [Pa]

% Valve/port timing areas for reference/plots
A_int_m2  = valve_area_1080(phi_deg, params.IVO_A, params.IVC_A, params.D_int,  params.Lift_int);
A_exh_m2  = valve_area_1080(phi_deg, params.EVO_C1, params.EVC_C1, params.D_exh, params.Lift_exh) ...
          + valve_area_1080(phi_deg, params.EVO_F1, params.EVC_F1, params.D_exh, params.Lift_exh);
A_scav_m2 = port_area(phi_deg, params.SP_open, params.SP_close, params.SP_perim, params.SP_hmax);

% Effective lifts from areas
Lift_int_m = zeros(size(phi_deg));
Lift_exh_m = zeros(size(phi_deg));
Lift_scav_m = zeros(size(phi_deg));
if params.D_int > 0
    Lift_int_m = A_int_m2 ./ (pi * params.D_int);
end
if params.D_exh > 0
    Lift_exh_m = A_exh_m2 ./ (pi * params.D_exh);
end
if params.SP_perim > 0
    Lift_scav_m = A_scav_m2 ./ params.SP_perim;
end

% Spark indicator (1 when the igniters fire, else 0) for plotting
spark = zeros(size(phi_deg));
flame = zeros(size(phi_deg));
phi_mod = mod(phi_deg, 1080);
if numel(phi_deg) > 1
    phi_step = max(abs(diff(phi_deg)));
else
    phi_step = 1;  % fall back to avoid zero-width smoothing
end
spark_events = [params.soc1, params.soc2];
spark_durations = [params.dur1, params.dur2];
for k = 1:numel(spark_events)
    soc = spark_events(k);
    dur = spark_durations(k);
    if ~isfinite(soc) || ~isfinite(dur) || dur <= 0
        continue;
    end
    soc_mod = mod(soc, 1080);
    delta = mod(phi_mod - soc_mod + 1080, 1080);  % [0,1080)
    in_window = delta <= dur;
    if ~any(in_window)
        continue;
    end
    % Smooth rise/fall portions, cap to duration and keep plateau at 1
    rise = max( min(0.1 * dur, 5), phi_step );
    rise = min(rise, dur);
    fall = rise;
    if rise + fall > dur
        rise = dur / 2;
        fall = dur - rise;
    end
    core = max(dur - rise - fall, 0);
    pulse = zeros(size(phi_deg));
    if rise > 0
        rise_mask = in_window & (delta < rise);
        tau = delta(rise_mask) / max(rise, eps);
        pulse(rise_mask) = 0.5 - 0.5*cos(pi * tau);  % smooth ramp 0→1
    end
    flat_mask = in_window & (delta >= rise) & (delta <= rise + core);
    pulse(flat_mask) = 1.0;
    if fall > 0
        fall_mask = in_window & (delta > rise + core) & (delta <= dur);
        tau_f = (delta(fall_mask) - (rise + core)) / max(fall, eps);
        pulse(fall_mask) = 0.5 - 0.5*cos(pi * (1 - tau_f));  % smooth 1→0
    end
    spark = max(spark, pulse);  % handle overlap between events gracefully
    
    flame_mask = in_window;
    if any(flame_mask)
        norm_delta = min(delta(flame_mask) ./ max(dur, eps), 1);
        flame_vals = 0.5 - 0.5*cos(pi * norm_delta);  % monotonic growth 0→1
        flame(flame_mask) = max(flame(flame_mask), flame_vals);
    end
end

%% Package
result.phi_deg           = phi_deg;
result.p_cyl             = p_cyl;
result.m_cyl             = m_cyl;
result.mO2               = mO2;
result.T_cyl             = T_cyl;
result.T_ind             = T_ind;
result.T_ind_pdv         = T_ind_pdv;
result.IMEP_bar          = IMEP_bar;
result.T_ind_mean        = T_ind_mean;
result.m_fuel_cycle_total= m_fuel_tot;
result.mO2_at_SOC2       = mO2_at_SOC2;
result.N_rpm             = params.N_rpm;

% Geometry + kinematics
result.V_m3                  = V_vec;
result.dVdphi_m3_per_deg     = dV_dphi;
result.x_mm                  = x_mm;
result.dx_dphi_m_per_deg     = dx_dphi;
result.ddx_dphi2_m_per_deg2  = ddx_dphi2;

% Forces/energetics
result.F_gas_N           = F_gas;
result.W_cum_J           = W_cum;
result.dQ_dphi_J_per_deg = dQ_dphi;
result.Qdot_wall_W       = Qdot_wall_W;
result.xburn1            = xburn1;
result.xburn2            = xburn2;
result.t_sec             = t_sec;
result.power_W           = power_W;
result.imep_inst_Pa      = imep_inst;

% Valve/port areas
result.A_int_m2          = A_int_m2;
result.A_exh_m2          = A_exh_m2;
result.A_scav_m2         = A_scav_m2;
result.L_int_m           = Lift_int_m;
result.L_exh_m           = Lift_exh_m;
result.L_scav_m          = Lift_scav_m;
result.spark             = spark;
result.flame             = flame;

end

%% ================= RHS: x' = dx/dphi ================================
function dx_dphi = engine_rhs_6stroke(phi_deg, x, params)
% State x = [p_cyl; m_cyl; mO2]
p_cyl = x(1);
m_cyl = max(x(2), 1e-9);
mO2   = max(min(x(3), m_cyl), 0.0);  % clamp 0 <= mO2 <= m

R      = params.R_mix;
cp     = params.cp;
% gamma  = params.gamma;
dt_dphi= 1/(6*params.N_rpm);

% Geometry & walls from geometry_6stroke
[V, Aw, dV_dph] = geom_eval(phi_deg, params);

% Gas temperature
T_cyl  = p_cyl * V / (m_cyl * R);

% Adjust cp and gamma - quick-and-safe surrogate for air (300–3000 K):
cp = cp + 0.1*(T_cyl - 300);      % J/kg/K, clamp [1005, 1400]
cp = min(max(cp,1005),1400);
gamma = cp/(cp - params.R_mix);

% ---------------- 1) Heat release (two Wiebes) ----------------
% Total chemical energy per 1080° split across burns
Q1 = params.eta_comb1 * params.m_fuel1    * params.LHV;  % [J/cycle]
% Total burn #2 energy command
Q2_cmd = params.eta_comb2 * params.m_fuel2_cmd * params.LHV;
% Window-wide O2 cap factor (updated each cycle); default to 1 if absent
if isfield(params,'lambda2')
    lambda2 = params.lambda2;
else
    lambda2 = 1.0;
end

% Wiebe derivatives (per degree)
dx1_dphi = wiebe_dxb_dphi(phi_deg, params.soc1, params.dur1, params.m1, params.a1);
dx2_shape= wiebe_dxb_dphi(phi_deg, params.soc2, params.dur2, params.m2, params.a2);

% Window-wide scaled energy rate for burn #2
dQchem_dphi = Q1 * dx1_dphi + (lambda2 * Q2_cmd) * dx2_shape;  % [J/deg]

% ---------------- 2) Wall heat loss (Woschni) ----------------
Qdot_w   = Qdot_wall(phi_deg, p_cyl, T_cyl, Aw, params);   % [W]
dQw_dphi = Qdot_w * dt_dphi;                               % [J/deg]

% ---------------- 3) Valve & port flows + algebraic scavenge -----------
[p_int,T_int,p_exh,T_exh,p_sc,T_sc] = deal(params.p_int,params.T_int, ...
                                           params.p_exh,params.T_exh, ...
                                           params.p_scav,params.T_scav);

% Areas (deg wrapped to 1080)
A_int = valve_area_1080(phi_deg, params.IVO_A, params.IVC_A, params.D_int, params.Lift_int);
% Exhaust split into blowdown in C and full exhaust in F (closed during D,E)
A_exh = valve_area_1080(phi_deg, params.EVO_C1, params.EVC_C1, params.D_exh, params.Lift_exh) ...
      + valve_area_1080(phi_deg, params.EVO_F1, params.EVC_F1, params.D_exh, params.Lift_exh);
A_scav_geom = port_area(phi_deg, params.SP_open, params.SP_close, params.SP_perim, params.SP_hmax);

% Raw nozzle flows (signed, + from upstream->downstream arg order)
mdot_int  = nozzle_mdot(p_int, T_int, p_cyl, T_cyl, A_int,  params.Cd_int,  gamma, R);
mdot_exh  = nozzle_mdot(p_cyl, T_cyl, p_exh, T_exh, A_exh,  params.Cd_exh,  gamma, R);
mdot_scav = nozzle_mdot(p_sc,  T_sc,  p_cyl, T_cyl, A_scav_geom, params.Cd_scav, gamma, R);

% Algebraic scavenge model (active only when ports are open)
eta_tr  = params.eta_tr;   % trapped fraction of fresh inflow
eta_mix = params.eta_mix;  % additional residual displacement per trapped inflow

mdot_in  = 0; Hdot_in  = 0;
mdot_out = 0; Hdot_out = 0;
% O2 flow accumulators (kg/s); convert to per-degree later
dO2_in_flow  = 0;
dO2_out_flow = 0;
Y_O2_cyl  = max(min(mO2/m_cyl,1),0);

% Intake valve (head)
if mdot_int > 0
    mdot_in  = mdot_in  + mdot_int;
    Hdot_in  = Hdot_in  + mdot_int * cp * T_int;
    % O2 entering with intake air
    dO2_in_flow = dO2_in_flow + params.YO2_air * mdot_int;
else
    mdot_out = mdot_out + (-mdot_int);
    Hdot_out = Hdot_out + (-mdot_int) * cp * T_cyl;
    % O2 leaving with backflow to intake
    dO2_out_flow = dO2_out_flow + Y_O2_cyl * (-mdot_int);
end

% Exhaust valve (head)
if mdot_exh > 0
    mdot_out = mdot_out + mdot_exh;
    Hdot_out = Hdot_out + mdot_exh * cp * T_cyl;
    % O2 leaving with exhaust
    dO2_out_flow = dO2_out_flow + Y_O2_cyl * mdot_exh;
else
    mdot_in  = mdot_in  + (-mdot_exh);
    Hdot_in  = Hdot_in  + (-mdot_exh) * cp * T_exh;
    % O2 entering with exhaust backflow (often ~0 unless YO2_exh>0)
    dO2_in_flow = dO2_in_flow + params.YO2_exh * (-mdot_exh);
end

% Scavenge ports with trapping & displacement
if A_scav_geom > 0
    if mdot_scav > 0
        % fresh scavenge air arrives from p_scav
        mdot_trap = eta_tr * mdot_scav;           % trapped in-cylinder
        % short-circuited fraction (1-eta_tr)*mdot_scav bypasses cylinder
        mdot_disp = eta_mix * mdot_trap;          % additional residual displacement (outflow)

        % Energy accounting:
        % - trapped inflow adds enthalpy at T_sc
        mdot_in  = mdot_in  + mdot_trap;
        Hdot_in  = Hdot_in  + mdot_trap * cp * T_sc;
        % - short-circuit does not change cylinder state (bypass)
        % - displacement expels cylinder gas at T_cyl
        mdot_out = mdot_out + mdot_disp;
        Hdot_out = Hdot_out + mdot_disp * cp * T_cyl;

        % Oxygen bookkeeping (fresh in, displaced out)
        dO2_in_flow  = dO2_in_flow  + params.YO2_air * mdot_trap;
        dO2_out_flow = dO2_out_flow + Y_O2_cyl * mdot_disp;
    else
        % backflow from cylinder to scavenge plenum (treat like outflow)
        mdot_out = mdot_out + (-mdot_scav);
        Hdot_out = Hdot_out + (-mdot_scav) * cp * T_cyl;

        dO2_out_flow = dO2_out_flow + Y_O2_cyl * (-mdot_scav);
    end
end

% Net mass & enthalpy per degree
dm_dphi   = (mdot_in - mdot_out) * dt_dphi;   % [kg/deg]
dHin_dphi = Hdot_in * dt_dphi;                % [J/deg]
dHout_dphi= Hdot_out* dt_dphi;                % [J/deg]

% Oxygen change per degree from flows
dmO2_flow_dphi = (dO2_in_flow - dO2_out_flow) * dt_dphi;  % [kg/deg]

% ---------------- 4) Oxygen consumption by burning ----------------
% Burn #1 O2 consumption follows Wiebe #1 shape (stoich):
dmO2_burn1_dphi = params.nu_O2 * params.m_fuel1 ...
    * dx1_dphi / max( trapz_lin(params.soc1, params.dur1, params.m1, params.a1), 1e-12);

% Burn #2 O2 consumption with the same window-wide scaling
dmO2_burn2_dphi = lambda2 * params.nu_O2 * params.m_fuel2_cmd ...
    * dx2_shape / max( trapz_lin(params.soc2, params.dur2, params.m2, params.a2), 1e-12);

% Total O2 consumption (nonzero only within each burn window)
if phi_deg < params.soc1 || phi_deg > params.soc1 + params.dur1
    dmO2_burn1_dphi = 0;
end
if phi_deg < params.soc2 || phi_deg > params.soc2 + params.dur2
    dmO2_burn2_dphi = 0;
end
dmO2_burn_dphi = dmO2_burn1_dphi + dmO2_burn2_dphi;

% ---------------- 5) Pressure ODE in phi-domain ------------------
dp_dphi = -gamma * p_cyl / V * dV_dph ...
          + (gamma - 1)/V * ( dQchem_dphi - dQw_dphi + dHin_dphi - dHout_dphi );

% States
dmO2_dphi = dmO2_flow_dphi - dmO2_burn_dphi;
dx_dphi = [dp_dphi; dm_dphi; dmO2_dphi];
end

%% =================== Combustion (Wiebe helper) ==========================
function dxb = wiebe_dxb_dphi(phi, soc, dur, m, a)
dxb = zeros(size(phi));
mask = (phi >= soc) & (phi <= soc+dur);
tau  = (phi(mask) - soc)/dur;
dxb(mask) = (m+1)*a/dur .* (tau.^m) .* exp(-a * tau.^(m+1));
end

function area = trapz_lin(soc, dur, m, a)
% normalize ∫ dxb dphi over the window (analytical would be dur*(1-exp(-a)))
% we keep numeric for robustness
phi = linspace(soc, soc+dur, max(20,ceil(dur)));
dxb = wiebe_dxb_dphi(phi, soc, dur, m, a);
area = trapz(phi, dxb);
if area <= 0, area = dur; end
end

%% =================== Heat transfer =====================================
function Qdot_w = Qdot_wall(phi_deg, p_cyl, T_cyl, Aw, params)
B      = params.B / 1e3; % convert mm to m
C1     = params.C1_wosch;
C2     = params.C2_wosch;
T_wall = params.T_wall;
p_cyl  = max(real(p_cyl), params.p_floor);
T_cyl  = max(real(T_cyl), params.T_floor);

% Segment-specific mean piston speed based on stroke used in the segment
ph = mod(phi_deg,1080);
dVdphi = params.geotab.dVdphi_s(ph);
Ap = pi*(B^2)/4;
dx_dphi = dVdphi / Ap;                      % [m/deg]
dphidt = 6*params.N_rpm;                    % deg/s
xdot = abs(dx_dphi .* dphidt);              % [m/s]
w  = C2 * max(xdot, 1e-3);                  % pseudo-velocity (avoid zero)

p_bar = p_cyl / 1e5;
h = C1 * p_bar.^0.8 .* T_cyl.^(-0.53) .* w.^0.8 .* B.^(-0.2);
Qdot_w = h .* Aw .* (T_cyl - T_wall);
end

%% =================== Flow areas ========================================
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

function A = port_area(phi_deg, open_deg, close_deg, perim, hmax)
% Scavenge port effective area ~ perimeter * effective open height (half-cosine)
ph = mod(phi_deg,1080);
o  = mod(open_deg,1080);
c  = mod(close_deg,1080);

if o < c
    active = (ph >= o) & (ph <= c);
    dph = c - o;
    tau = (ph - o)/dph;
else
    active = (ph >= o) | (ph <= c);
    dph = mod(c - o,1080);
    tau = zeros(size(ph));
    m1 = active & (ph >= o);
    m2 = active & (ph <= c);
    tau(m1) = (ph(m1) - o)/dph;
    tau(m2) = (ph(m2) + (1080 - o))/dph;
end

h = zeros(size(ph));
h(active) = 0.5*hmax*(1 - cos(2*pi*tau(active)));
A = perim .* h;
end

%% =================== Nozzle model (with upstream T) =====================
function mdot = nozzle_mdot(p1, T1, p2, T2, A, Cd, gamma, R)
if A <= 0 || p1<=0 || p2<=0 || T1<=0 || T2<=0
    mdot = 0; return;
end
if abs(p1-p2)/max(p1,p2) < 1e-4
    mdot = 0; return;
end

if p1 > p2
    pu = p1; Tu = T1; pd = p2; sign_flow = 1;
else
    pu = p2; Tu = T2; pd = p1; sign_flow = -1;
end

pr = pd/pu;
pcrit = (2/(gamma+1))^(gamma/(gamma-1));

if pr <= pcrit
    mdot_mag = Cd * A * pu * sqrt(gamma/(R*Tu)) * (2/(gamma+1))^((gamma+1)/(2*(gamma-1)));
else
    term = pr^(2/gamma) - pr^((gamma+1)/gamma);
    if term < 0, term = 0; end
    mdot_mag = Cd * A * pu * sqrt( 2*gamma/(R*Tu*(gamma-1)) * term );
end
mdot = sign_flow * mdot_mag;
end

%% =================== Torque ============================================
function T_gas = gas_torque(phi_deg, p_cyl, params)
% Torque from gas force using exact kinematics: tau = F * (dx/dtheta in m/rad)
% where dx/dtheta from geometry_6stroke is stored in [m/deg] already

B_m = params.B * 1e-3;              % convert bore from mm -> m
Ap  = pi * (B_m^2) / 4;             % piston area [m^2]
p_crank = params.p_crank;

F_gas = (p_cyl - p_crank) * Ap;     % [N] downward force on piston

% Get dx/dphi from geometry tables (already in [m/deg])
ph = mod(phi_deg, 1080);
dxdphi_m_per_deg = params.geotab.dxdphi_s(ph);  % [m/deg]

% Convert to [m/rad] and compute torque: tau = F * (dx/dtheta in m/rad)
dxdphi_m_per_rad = dxdphi_m_per_deg * (180/pi);  % [m/rad]
T_gas = F_gas .* dxdphi_m_per_rad;  % [Nm]
end

function T_gas = gas_torque_pdv(phi_deg, p_cyl, params)
% Torque from p*dV: T = (p - p_crank) * dV/dtheta, with dV/dtheta in [m^3/rad]
ph = mod(phi_deg,1080);
dVdphi_m3_per_deg = params.geotab.dVdphi_s(ph);          % [m^3/deg]
dVdtheta_m3_per_rad = dVdphi_m3_per_deg * (180/pi);      % [m^3/rad]
T_gas = (p_cyl - params.p_crank) .* dVdtheta_m3_per_rad; % [Nm]
end

%% =================== Geometry lookup via pre-tabulation ==================
function [V_m3, A_m2, dVdphi_m3_per_deg] = geom_eval(phi_deg, params)
% Evaluate V, A, dV/dphi from precomputed interpolants (consistent V/dV).
if ~isfield(params,'geotab') || ~isfield(params.geotab,'V_s')
    error('geom_eval:MissingGeoTab', 'Geometry tables not initialized in params.geotab.*');
end
ph = mod(phi_deg, 1080);
gt = params.geotab;
V_m3              = gt.V_s(ph);
A_m2              = gt.A_s(ph);
dVdphi_m3_per_deg = gt.dVdphi_s(ph);
end

%% =================== State limiter ========================================
function X = enforce_state_limits(X_in, params)
X = real(X_in);
p_floor = params.p_floor;
m_floor = params.m_floor;

if isvector(X)
    col = X(:);
    col(1) = max(col(1), p_floor);
    col(2) = max(col(2), m_floor);
    col(3) = min(max(col(3), 0), col(2));
    if isrow(X_in)
        X = col.';
    else
        X = col;
    end
    return;
end

X(:,1) = max(X(:,1), p_floor);
X(:,2) = max(X(:,2), m_floor);
X(:,3) = min(max(X(:,3), 0), X(:,2));
end

%% =================== Utility: parameter validation =======================
function ensure_fields(s, required, contextLabel)
required = cellstr(required);
missing = required(~cellfun(@(name) isfield(s, name) && ~isempty(s.(name)), required));
if ~isempty(missing)
    error('engine_cycle_steady_6stroke:MissingParam', ...
        'Missing %s field(s): %s', contextLabel, strjoin(missing, ', '));
end
end