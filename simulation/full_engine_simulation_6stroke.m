function result = full_engine_simulation_6stroke(result, params, plot_diag)
%FULL_ENGINE_SIMULATION_6STROKE  8-DOF flat-6 engine dynamics with TE.
%
%   result = full_engine_simulation_6stroke(result, params, plot_diag)
%
%   Multi-cylinder model:
%     DOFs (angles):
%       q1..q6 : theta_p1..theta_p6  [rad]  planet gear angles (one per cyl)
%       q7     : theta_c             [rad]  carrier / crank angle
%       q8     : theta_o             [rad]  output shaft angle
%
%     States:
%       z = [theta_p1..theta_p6; theta_c; theta_o; ...
%            omega_p1..omega_p6; omega_c; omega_o]   (16 x 1)
%
%   The engine is driven by 6 cylinders, each with pressure p_i(phi) obtained
%   from the thermodynamic SINGLE-cylinder simulation, phase shifted by 180°
%   in CRANK angle phi (phi = theta_c in deg, wrapped to [0,1080)).
%
%   Pressure -> piston force -> generalized torques using piston kinematics
%   (via generalized_gas_force.m).
%   Mesh stiffness acts between EACH planet and its local carrier arm
%   (via generalized_mesh_force.m).
%   A linear DMF acts between global carrier and output
%   (via generalized_DMF_force.m).
%   Output shaft is loaded by a viscous torque sized to match steady torque.
%
%   Inputs:
%     result : struct from thermodynamic single-cylinder sim, expected fields:
%                .phi_deg   [deg]  0..1080 crank angle grid (single cyl)
%                .p_cyl     [Pa]   cylinder pressure vs phi (single cyl)
%                .T_ind     [Nm]   indicated torque vs phi (single cyl)
%                (optionally .N_rpm)
%
%     params : struct; used fields include:
%       --- engine speed & sim options ---
%         .N_rpm                 nominal engine rpm (if not in result)
%         .sim.n_cycles_dmf      # cycles for averaging (default 3)
%         .sim.ode_rel_tol       (default 1e-7)
%         .sim.ode_abs_tol       (default 1e-9)
%         .sim.ode_max_step      (default 0.005)
%
%       --- geometry & inertia (single cylinder, reused for all 6) ---
%         .L, .delta, .a         [m]   piston-planet-carrier geometry
%         .m_p                   [kg]  piston mass (incl. conrod)
%         .m_pg                  [kg]  planet gear mass (translation)
%         .J_p                   [kg m^2] planet gear inertia
%         .J_c                   [kg m^2] carrier + inner DMF inertia
%         .J_o                   [kg m^2] output + outer DMF inertia
%         .Rp, .Rr               [m]   planet & ring pitch radii
%
%       --- piston area ---
%         .A_piston [m^2] or .B [m] bore diameter (used by generalized_gas_force)
%
%       --- gear mesh (constant for now) ---
%         .k_mesh   [N/m]      mesh stiffness along mesh deflection y
%         .c_mesh   [N·s/m]    mesh damping along y
%
%       --- DMF ---
%         .k_DMF    [N·m/rad]
%         .c_DMF    [N·m·s/rad]
%         .phi_DMF0 [rad]  (optional, default 0)
%
%       --- losses ---
%         .loss.T_fric_const
%         .loss.T_fric_visc
%         .loss.T_acc_const
%         .loss.T_acc_visc
%
%   Output (added to result):
%     .T_gbx_in  [Nm]  cycle-averaged DMF output torque vs phi_deg grid
%     .P_gbx_in  [W]   cycle-averaged DMF output power vs phi_deg grid
%     .dyn.*           time histories for diagnostics (angles, speeds, torques)
%
%   Dependencies (separate m-files):
%     piston_kinematics.m
%     generalized_gas_force.m
%     generalized_mesh_force.m
%     generalized_DMF_force.m
% -------------------------------------------------------------------------
% Bart Blockmans, 2024 - bart@blockmans.net
% -------------------------------------------------------------------------


if nargin < 3; plot_diag = 0; end

% Convert parameters to SI units
params = convert_to_SI(params);

% ---------------------------------------------------------------------
% 0. Basic crank-angle grid and single-cylinder data
% ---------------------------------------------------------------------
phi_single = result.phi_deg(:);      % [deg], 0..1080
T1_ind     = result.T_ind(:);        % [Nm], single-cylinder indicated torque
p1_single  = result.p_cyl(:);        % [Pa], single-cylinder pressure

% Engine speed
Nrpm   = getfield_or(result, 'N_rpm', getfield_or(params, 'N_rpm', 2000));
omega0 = Nrpm * 2*pi/60;             % [rad/s], reference speed

% 6-cylinder phasing in crank angle (allow per-cylinder offsets in degrees)
cyl_offset_deg = getfield_or(params, 'cyl_offset', zeros(1,6));
if isscalar(cyl_offset_deg)
    cyl_offset_deg = repmat(cyl_offset_deg, 1, 6);
end
% Ensure exactly 6 elements (pad or truncate) and row shape
cyl_offset_deg = reshape(cyl_offset_deg(1:min(numel(cyl_offset_deg),6)), 1, []);
if numel(cyl_offset_deg) < 6
    cyl_offset_deg = [cyl_offset_deg, zeros(1, 6-numel(cyl_offset_deg))];
end
phase_deg_vec = (0:5)*180 + cyl_offset_deg;    % [deg]
phase_rad_vec = phase_deg_vec(:) * pi/180;     % [rad]; column vector used for local carrier angles

% ---------------------------------------------------------------------
% 1. Build flat-6 SUM torque (for sizing viscous load, diagnostics)
% ---------------------------------------------------------------------
Tsum = zeros(size(T1_ind));
for iCyl = 1:6
    % Note the minus sign: each cylinder is delayed by its phase (deg)
    ph_shift = phase_deg_vec(iCyl);
    Tsum = Tsum + interp1_periodic(phi_single, T1_ind, wrap_phi(phi_single - ph_shift));
end
result.T_engine_sum = Tsum;

% Periodic pchip for summed torque vs crank angle
phi_ext     = [phi_single; phi_single(2:end)+1080];
Tsum_ext    = [Tsum;       Tsum(2:end)];
Tsum_pchip  = @(phq) interp1(phi_ext, Tsum_ext, wrap_phi(phq), 'pchip');

% ---------------------------------------------------------------------
% 2. Single-cylinder pressure interpolant (periodic over 1080°)
% ---------------------------------------------------------------------
phi_p_ext = [phi_single; phi_single(2:end)+1080];
p1_ext    = [p1_single;  p1_single(2:end)];
p1_pchip  = @(phq) interp1(phi_p_ext, p1_ext, wrap_phi(phq), 'pchip');

% Helper: flat-6 pressures as function of GLOBAL crank angle (deg)
    function p6 = flat6_pressures(phi_deg)
        % phi_deg: scalar or vector (global)
        phases = phase_deg_vec;  % cylinder phase shifts [deg], includes params.cyl_offset
        p6 = zeros(6, numel(phi_deg));
        for iCyl = 1:6
            % Local angle φ_loc,i = φ_global − phase_i
            p6(iCyl,:) = p1_pchip(phi_deg - phases(iCyl));
        end
    end

% ---------------------------------------------------------------------
% 3. Loss models (friction + accessories)
% ---------------------------------------------------------------------
L = getfield_or(params,'loss', struct());
T_fric = @(w) getfield_or(L,'T_fric_const',0) + getfield_or(L,'T_fric_visc',0)*w;
T_acc  = @(w) getfield_or(L,'T_acc_const',0)  + getfield_or(L,'T_acc_visc',0)*w;

% ---------------------------------------------------------------------
% 4. Viscous output load torque sizing
% ---------------------------------------------------------------------
% Mean engine torque after losses at nominal speed, using Tsum (flat-6)
Tmean_losses = T_fric(omega0) + T_acc(omega0);
Tmean_engine = mean(Tsum) - Tmean_losses;

% Choose viscous coefficient so that D_load*omega0 ≈ Tmean_engine
D_load = Tmean_engine / max(omega0, 1e-9);
T_load = @(wg) D_load * wg;   % [N·m]

% ---------------------------------------------------------------------
% 5. Simulation controls
% ---------------------------------------------------------------------
S        = getfield_or(params,'sim', struct());
n_cycles = getfield_or(S,'n_cycles_dmf', 3);
RelTol   = getfield_or(S,'ode_rel_tol', 1e-7);
AbsTol   = getfield_or(S,'ode_abs_tol', 1e-9);
ode_max_step = getfield_or(S,'ode_max_step', 0.005);

deg_per_sec = 6*Nrpm;                  % [deg/s], 3 rev per 1080°
tf = (n_cycles + 1) * (1080/deg_per_sec);  % +1 cycle for transients

% Progress output function for ODE integration
progress_thresholds = 0.1:0.1:1.0;  % 10%, 20%, ..., 100%
function status = ode_progress_output(t, ~, flag)
    persistent printed
    thresholds = progress_thresholds;  % Access from parent scope
    t_final = tf;  % Access tf from parent scope
    
    if strcmp(flag, 'init')
        % Initialize: reset printed flags
        printed = false(size(thresholds));
        status = 0;
    elseif isempty(flag) || strcmp(flag, '')
        % During integration: check progress
        if ~isempty(t)
            current_progress = t(end) / t_final;
            for idx = 1:numel(thresholds)
                if current_progress >= thresholds(idx) && ~printed(idx)
                    fprintf('      %d%% of the simulation completed.\n', round(thresholds(idx)*100));
                    printed(idx) = true;
                end
            end
        end
        status = 0;
    else
        % 'done' flag or other
        status = 0;
    end
end

opts = odeset('RelTol',RelTol,'AbsTol',AbsTol, 'MaxStep', ode_max_step, ...
    'OutputFcn', @ode_progress_output);

% ---------------------------------------------------------------------
% 6. Initial conditions (consistent with geometry & ideal mesh at t=0)
% ---------------------------------------------------------------------
theta_c0 = 0;       % global carrier angle at t=0
theta_o0 = 0;       % choose same reference for output (zero DMF twist)

R_p = params.Rp;
R_r = params.Rr;
k_gr = (R_r - R_p)/R_p;   % gear ratio factor used for ω_p and θ_p relation

theta_p0 = zeros(6,1);

for i = 1:6
    % Local carrier angle of cylinder i at t=0:
    %   θ_c,loc,i = θ_c(global) − phase_rad_i
    theta_ci0 = theta_c0 - phase_rad_vec(i);

    % Planet angle consistent with 6-stroke geometry relation used
    % in geometry_6stroke / piston_kinematics:
    %   theta_p = -k_gr * (theta_c_loc - 3π/2) + π/2
    theta_p0(i) = -k_gr * (theta_ci0 - 3*pi/2) + pi/2;
end

% Initial speeds:
omega_c0 = omega0;          % carrier defines nominal engine speed
omega_o0 = omega0;          % same initial speed on output (zero DMF twist)
omega_p0 = -k_gr * omega_c0 * ones(6,1);  % ideal internal gear kinematics

% State vector: [theta_p1..theta_p6; theta_c; theta_o; omega_p1..omega_p6; omega_c; omega_o]
z0 = [theta_p0; theta_c0; theta_o0; omega_p0; omega_c0; omega_o0];

% ---------------------------------------------------------------------
% 7. Integrate 8-DOF dynamics
% ---------------------------------------------------------------------
odefun = @(t,z) rhs_8dof(t, z);  % nested function below

% Display
disp('   Integrating the engine equations of motion...')
disp('   ');
disp('      Simulation started...');

% Solve EOMs
[t, z] = ode15s(odefun, [0 tf], z0, opts);

% Display
disp(' ')
disp('   Engine equations of motion integrated successfully.')
disp('   ');
disp('   Post-processing the simulation results...')

% Unpack solution
theta_p = z(:,1:6);    % [n x 6]
theta_c = z(:,7);      % [n x 1]
theta_o = z(:,8);      % [n x 1]
omega_p = z(:,9:14);   % [n x 6]
omega_c = z(:,15);
omega_o = z(:,16);

% ---------------------------------------------------------------------
% 8. Reconstruct torques for diagnostics & averaging
% ---------------------------------------------------------------------
% Crank angle in degrees (unwrapped)
tot_phi_deg = rad2deg(theta_c);     % unwrapped crank angle

% Engine torque vs crank angle (for comparison, based on thermo Tsum)
phi_c_mod = rad2deg(mod(theta_c, 6*pi));     % [0,1080)
T_eng_sum = Tsum_pchip(phi_c_mod) - T_fric(omega_c) - T_acc(omega_c);

% DMF torque and mesh torques via helper / external functions
nstep = numel(t);
T_DMF = zeros(nstep,1);
T_mesh_total = zeros(nstep,1);   % sum over all 6 meshes

% Dynamic gas torque
T_gas_dyn = zeros(nstep,1);

% Transmission error (TE) for each planet gear [m]
TE = zeros(nstep, 6);

% Reference angles for TE computation (from generalized_mesh_force.m)
% R_p and R_r are already defined earlier in the function
theta_c0 = 3*pi/2;
theta_p0 = pi/2;

for k = 1:nstep
    th_c = theta_c(k);
    w_c  = omega_c(k);
    th_o = theta_o(k);
    w_o  = omega_o(k);

    % DMF torque between carrier and output
    Q_DMF_local = generalized_DMF_force(params, th_c, w_c, th_o, w_o);
    % Q_DMF_local = [0; Q_theta_c; Q_theta_o], T_DMF acts equal & opposite
    T_DMF(k) = -Q_DMF_local(2);   % torque transmitted to gearbox/input

    % Mesh torques for each planet (with local carrier angles)
    T_mesh_k = 0;
    for iCyl = 1:6
        th_pi = theta_p(k,iCyl);
        w_pi  = omega_p(k,iCyl);

        % local carrier angle for this mesh
        theta_ci = th_c - phase_rad_vec(iCyl);

        Q_mesh_local = generalized_mesh_force(params, th_pi, w_pi, theta_ci, w_c);

        % Q_mesh_local(1) is torque on planet, Q_mesh_local(2) on carrier
        T_mesh_k = T_mesh_k + Q_mesh_local(2);   % sign consistent with planet DOF
        
        % Transmission error along line of action [m]
        % Formula from generalized_mesh_force.m line 67
        TE(k, iCyl) = (R_p - R_r) .* (theta_ci - theta_c0) - R_p .* (th_pi - theta_p0);
    end
    T_mesh_total(k) = T_mesh_k;

    % Global crank angle in deg
    phi_c_deg = rad2deg(mod(th_c, 6*pi));
    p6_now    = flat6_pressures(phi_c_deg);

    Q_gas_k = zeros(8,1);
    for iCyl = 1:6
        th_pi   = theta_p(k,iCyl);
        p_i     = p6_now(iCyl);
        theta_ci= th_c - phase_rad_vec(iCyl);

        Q_gas_local = generalized_gas_force(params, th_pi, theta_ci, p_i);
        Q_gas_k(iCyl) = Q_gas_k(iCyl) + Q_gas_local(1);
        Q_gas_k(7)    = Q_gas_k(7)    + Q_gas_local(2);
    end

    % Gas torque on carrier DOF
    T_gas_dyn(k) = Q_gas_k(7);

end

% Gearbox input power
P_gbx = T_DMF .* omega_o;

% ---------------------------------------------------------------------
% 9. Average steady cycles onto uniform φ grid
% ---------------------------------------------------------------------
phi_q = phi_single(:);   % reuse original uniform grid [0..1080]

[T_gbx_in, P_gbx_in] = average_steady_cycles( ...
    tot_phi_deg, phi_q, T_DMF, P_gbx);

T_gas_mean = average_steady_cycles( ...
    rad2deg(theta_c), phi_q, T_gas_dyn);

% Average transmission error for each planet gear over steady cycles
TE_mm = zeros(numel(phi_q), 6);
for iCyl = 1:6
    TE_mm(:, iCyl) = average_steady_cycles( ...
        rad2deg(theta_c), phi_q, TE(:, iCyl));
end
% Convert from meters to millimeters
TE_mm = TE_mm * 1000;

% Sum of all 6 transmission errors
TE_sum = sum(TE_mm, 2);

% Peak-to-peak variation (de-trended) for each planet gear
% Use windowed averaging to remove slow-varying trend (due to torque changes)
% while preserving fast-varying component (due to mesh stiffness variation)
TE_p2p = zeros(size(TE_mm));
window_deg = 30;  % Window size in degrees for averaging (adjust as needed)
dphi = phi_q(2) - phi_q(1);  % Step size in degrees
window_samples = max(1, round(window_deg / dphi));  % Window size in samples

for iCyl = 1:6
    TE_signal = TE_mm(:, iCyl);
    n = numel(TE_signal);
    
    % Compute windowed average with circular padding
    % (handles periodic boundary at 0/1080 degrees)
    TE_windowed_avg = zeros(n, 1);
    window_half = floor(window_samples / 2);
    
    for i = 1:n
        % Indices for window (with circular wrapping)
        idx_range = (i - window_half):(i + window_half);
        idx_vec = mod(idx_range - 1, n) + 1;  % Wrap to [1, n] range
        TE_windowed_avg(i) = mean(TE_signal(idx_vec));
    end
    
    % De-trended signal: subtract windowed average
    TE_p2p(:, iCyl) = TE_signal - TE_windowed_avg;
end

result.T_gbx_in = T_gbx_in(:);
result.P_gbx_in = P_gbx_in(:);
result.T_gas_mean = T_gas_mean(:);
result.TE_mm = TE_mm;
result.TE_sum = TE_sum(:);
result.TE_p2p = TE_p2p;

% ---------------------------------------------------------------------
% 10. Store diagnostics
% ---------------------------------------------------------------------
result.dyn.t        = t;
result.dyn.theta_p  = theta_p;
result.dyn.theta_c  = theta_c;
result.dyn.theta_o  = theta_o;
result.dyn.omega_p  = omega_p;
result.dyn.omega_c  = omega_c;
result.dyn.omega_o  = omega_o;
result.dyn.T_DMF    = T_DMF;
result.dyn.T_mesh   = T_mesh_total;
result.dyn.T_eng    = T_eng_sum;

% Display
disp('   Simulation results post-processed successfully.')
disp('   ');

% ---------------------------------------------------------------------
% 11. Optional diagnostic plots
% ---------------------------------------------------------------------
if plot_diag
    figure('Name','Engine dynamics diagnostics','Color','w');
    tiledlayout(3,2,'TileSpacing','compact','Padding','compact');

    % Speeds
    nexttile; hold on;
    plot(t, omega_c,'LineWidth',1.2,'DisplayName','\omega_c');
    plot(t, omega_o,'LineWidth',1.2,'DisplayName','\omega_o');
    plot(t, omega_p(:,1),'--','LineWidth',1.0,'DisplayName','\omega_{p1}');
    hold off; grid on; box on;
    xlabel('t [s]'); ylabel('\omega [rad/s]');
    title('Angular speeds'); legend('Location','best');

    % DMF torque
    nexttile; hold on;
    plot(t, T_DMF,'LineWidth',1.2);
    grid on; box on;
    xlabel('t [s]'); ylabel('T_{DMF} [N·m]');
    title('DMF torque');
    hold off;

    % Mesh torque total
    nexttile; hold on;
    plot(t, T_mesh_total,'LineWidth',1.2);
    grid on; box on;
    xlabel('t [s]'); ylabel('T_{mesh,total} [N·m]');
    title('Sum of mesh torques (all planets)');
    hold off;

    % Engine torque vs DMF torque
    nexttile; hold on;
    plot(t, T_eng_sum,'LineWidth',1.2,'DisplayName','T_{eng,sum}');
    plot(t, T_DMF,   'LineWidth',1.2,'DisplayName','T_{DMF}');
    grid on; box on;
    xlabel('t [s]'); ylabel('Torque [N·m]');
    title('Engine vs DMF torque'); legend('Location','best');
    hold off;

    % Output power
    nexttile; hold on;
    plot(t, P_gbx,'LineWidth',1.2);
    grid on; box on;
    xlabel('t [s]'); ylabel('P_{gbx} [W]');
    title('Gearbox input power');
    hold off;

    % Crank speed vs target
    nexttile; hold on;
    plot(t, omega_c,'LineWidth',1.2,'DisplayName','\omega_c');
    yline(omega0,'k--','LineWidth',1.2,'DisplayName','\omega_0 target');
    grid on; box on;
    xlabel('t [s]'); ylabel('\omega_c [rad/s]');
    title('Crank speed around target'); legend('Location','best');
    hold off;

    % Gas torque vs. thermo torque
    Tsum_flat6 = result.T_engine_sum(:);
    figure; hold on; box on; grid on;
    plot(phi_q, Tsum_flat6, 'k-', 'DisplayName','Tsum thermo');
    plot(phi_q, T_gas_mean, 'r--', 'DisplayName','Tgas dyn (carrier)');
    xlabel('\phi_c [deg]');
    ylabel('Torque [Nm]');
    legend('Location','best');
    title('Flat-6 gas torque: thermo vs dynamic');
    
end

% =====================================================================
% Nested: RHS of 8-DOF ODE system
% =====================================================================
function dzdt = rhs_8dof(~, z)
    % Unpack state (GLOBAL DOFs)
    theta_p = z(1:6);    % planet angles
    theta_c = z(7);      % carrier / crank
    theta_o = z(8);      % output shaft
    omega_p = z(9:14);   % planet speeds
    omega_c = z(15);     % carrier speed
    omega_o = z(16);     % output speed

    % Mass matrix and Coriolis for 8 DOFs (they handle local carrier angles)
    M8 = mass_matrix_8dof(params, theta_p, theta_c);
    C8 = coriolis_8dof(params, theta_p, omega_p, theta_c, omega_c);

    % Flat-6 cylinder pressures at current GLOBAL crank angle
    phi_c_deg = rad2deg(mod(theta_c, 6*pi));   % [0,1080)
    p6_now    = flat6_pressures(phi_c_deg);    % [6 x 1]

    % -------------------- Gas generalized forces ------------------
    Q_gas = zeros(8,1);
    for iCyl = 1:6
        th_pi = theta_p(iCyl);
        p_i   = p6_now(iCyl);

        % Local carrier angle for this cylinder:
        %   θ_c,loc,i = θ_c(global) − phase_rad_i
        theta_ci = theta_c - phase_rad_vec(iCyl);

        % External generalized_gas_force uses piston_kinematics internally
        % with (theta_p, theta_c_loc, p_i)
        Q_gas_local = generalized_gas_force(params, th_pi, theta_ci, p_i);

        % Map [planet; carrier; output] -> global 8D DOFs
        Q_gas(iCyl) = Q_gas(iCyl) + Q_gas_local(1);
        Q_gas(7)    = Q_gas(7)    + Q_gas_local(2);
        % Q_gas_local(3) is zero (no direct gas torque on output shaft)
    end

    % -------------------- Mesh generalized forces -----------------
    Q_mesh = zeros(8,1);
    for iCyl = 1:6
        th_pi = theta_p(iCyl);
        w_pi  = omega_p(iCyl);

        % Local carrier arm angle for this planet
        theta_ci = theta_c - phase_rad_vec(iCyl);

        Q_mesh_local = generalized_mesh_force(params, th_pi, w_pi, theta_ci, omega_c);

        Q_mesh(iCyl) = Q_mesh(iCyl) + Q_mesh_local(1);  % planet DOF
        Q_mesh(7)    = Q_mesh(7)    + Q_mesh_local(2);  % carrier DOF
    end

    % -------------------- DMF, losses, load -----------------------
    % DMF forces between global carrier and output
    Q_DMF_local = generalized_DMF_force(params, theta_c, omega_c, theta_o, omega_o);
    Q_DMF = zeros(8,1);
    Q_DMF(7) = Q_DMF_local(2);
    Q_DMF(8) = Q_DMF_local(3);

    % Loss torques on carrier
    T_loss = T_fric(omega_c) + T_acc(omega_c);
    Q_loss = zeros(8,1);
    Q_loss(7) = -T_loss;

    % Viscous load on output
    T_load_now = T_load(omega_o);
    Q_load = zeros(8,1);
    Q_load(8) = -T_load_now;

    % Total generalized forces
    Q_total = Q_gas + Q_mesh + Q_DMF + Q_loss + Q_load;

    % Accelerations
    qdd = M8 \ (Q_total - C8);   % [8 x 1]

    dzdt = zeros(16,1);
    dzdt(1:8)   = [omega_p; omega_c; omega_o];  % angle derivatives
    dzdt(9:16)  = qdd;                          % angular accelerations
end % rhs_8dof

end % full_engine_simulation_6stroke

% ========================================================================
% 8x8 Mass matrix for 6 planets + carrier + output
% ========================================================================
function M8 = mass_matrix_8dof(params, theta_p_vec, theta_c)
% theta_p_vec : [6x1]  global planet angles
% theta_c     : scalar global carrier angle

m_p  = params.mp;
m_pg = params.mpg;
J_p  = params.Jp;
J_c  = params.Jc;
J_o  = params.Jo;
a    = params.a;

% Per-cylinder phase offsets [rad] from perfect 180° spacing
cyl_offset_deg = getfield_or(params, 'cyl_offset', zeros(1,6));
if isscalar(cyl_offset_deg)
    cyl_offset_deg = repmat(cyl_offset_deg, 1, 6);
end
cyl_offset_deg = reshape(cyl_offset_deg(1:min(numel(cyl_offset_deg),6)), 1, []);
if numel(cyl_offset_deg) < 6
    cyl_offset_deg = [cyl_offset_deg, zeros(1, 6-numel(cyl_offset_deg))];
end
phase_rad_vec_local = ((0:5)*180 + cyl_offset_deg) * pi/180;

M8 = zeros(8,8);

for i = 1:6
    theta_pi = theta_p_vec(i);
    % Local carrier angle for this cylinder (consistent with RHS):
    %   θ_c,loc,i = θ_c(global) − phase_i
    theta_ci = theta_c - phase_rad_vec_local(i);

    % Use piston_kinematics to get Ac, Ap (vel/acc arguments = 0)
    [~, ~, ~, Ac, Ap] = piston_kinematics(params, ...
        theta_ci, 0, 0, ...
        theta_pi, 0, 0);

    % Piston + planet rotational inertia contribution
    M8(i,i) = M8(i,i) + m_p * Ap.^2 + J_p;
    M8(7,7) = M8(7,7) + m_p * Ac.^2;
    M8(i,7) = M8(i,7) + m_p * Ap.*Ac;
    M8(7,i) = M8(7,i) + m_p * Ap.*Ac;  % symmetry

    % Planet translational inertia from orbiting around carrier at radius a
    M8(7,7) = M8(7,7) + m_pg * a.^2;
end

% Add carrier + output inertias
M8(7,7) = M8(7,7) + J_c;
M8(8,8) = M8(8,8) + J_o;
end

% ========================================================================
% 8x1 Coriolis / centrifugal vector for 8-DOF system
% ========================================================================
function C8 = coriolis_8dof(params, theta_p_vec, omega_p_vec, theta_c, omega_c)
% Only piston mass generates configuration-dependent M; pure rotations
% (J_p, J_c, J_o, translational m_pg a^2) have constant inertia.

m_p = params.mp;
C8  = zeros(8,1);

% Per-cylinder phase offsets [rad] from perfect 180° spacing
cyl_offset_deg = getfield_or(params, 'cyl_offset', zeros(1,6));
if isscalar(cyl_offset_deg)
    cyl_offset_deg = repmat(cyl_offset_deg, 1, 6);
end
cyl_offset_deg = reshape(cyl_offset_deg(1:min(numel(cyl_offset_deg),6)), 1, []);
if numel(cyl_offset_deg) < 6
    cyl_offset_deg = [cyl_offset_deg, zeros(1, 6-numel(cyl_offset_deg))];
end
phase_rad_vec_local = ((0:5)*180 + cyl_offset_deg) * pi/180;

for i = 1:6
    theta_pi = theta_p_vec(i);
    omega_pi = omega_p_vec(i);

    % Local carrier angle for this cylinder
    theta_ci = theta_c - phase_rad_vec_local(i);

    % Get Ac, Ap and second derivatives from piston_kinematics
    [~, ~, ~, Ac, Ap, Bcc, Bpp, Bcp] = piston_kinematics(params, ...
        theta_ci, omega_c, 0, ...
        theta_pi, omega_pi, 0);

    % S = Bcc*omega_c^2 + 2*Bcp*omega_c*omega_pi + Bpp*omega_pi^2
    S_i = Bcc.*omega_c.^2 + 2.*Bcp.*omega_c.*omega_pi + Bpp.*omega_pi.^2;

    C8(i) = C8(i) + m_p .* Ap .* S_i;
    C8(7) = C8(7) + m_p .* Ac .* S_i;
end
end

% ========================================================================
% Helper functions (in-file utilities)
% ========================================================================
function y = getfield_or(s, name, defaultVal)
if isfield(s, name) && ~isempty(s.(name))
    y = s.(name);
else
    y = defaultVal;
end
end

function phw = wrap_phi(ph)
phw = mod(ph,1080);
phw(phw<0) = phw(phw<0) + 1080;
end

function y = interp1_periodic(x, y0, xi)
% periodic over 1080 deg with a non-duplicated seam and unique sample points
xx = [x; x(2:end)+1080];
yy = [y0; y0(2:end)];
[xxu, ia] = unique(xx, 'stable');
yyu = yy(ia);
y  = interp1(xxu, yyu, wrap_phi(xi), 'pchip');
end

function varargout = average_steady_cycles(tot_phi_deg, phi_q, varargin)
%AVERAGE_STEADY_CYCLES  Average one or more signals over complete 1080° cycles.
% Usage:
%   [y1_avg, y2_avg, ...] = average_steady_cycles(tot_phi_deg, phi_q, y1, y2, ...)
% Discards the first complete cycle; averages remaining cycles.

phi_q = phi_q(:);
if ~isempty(phi_q) && abs(phi_q(end) - 1080) < 1e-9
    phi_q_int = phi_q(1:end-1);
else
    phi_q_int = phi_q;
end

tot_phi = tot_phi_deg(:);
phi_min = tot_phi(1);
phi_max = tot_phi(end);

j_first = ceil((phi_min) / 1080);
j_last  = floor((phi_max - 1e-9) / 1080) - 1;

m = numel(varargin);
accum = cell(1,m);
for i = 1:m
    accum{i} = zeros(size(phi_q_int));
end
n_avg = 0;

for j = j_first+1 : j_last  % discard first complete cycle
    a = j*1080; b = (j+1)*1080;
    mask = (tot_phi >= a) & (tot_phi < b);
    if nnz(mask) < 2, continue; end
    phi_cyc = tot_phi(mask) - a;                 % [0,1080)
    [phi_u, ia] = unique(phi_cyc, 'stable');     % enforce unique φ
    for i = 1:m
        yi = varargin{i};
        yi = yi(:);
        yi_cyc = yi(mask);
        yi_u = yi_cyc(ia);
        accum{i} = accum{i} + interp1(phi_u, yi_u, phi_q_int, 'pchip','extrap');
    end
    n_avg = n_avg + 1;
end

for i = 1:m
    if n_avg > 0
        y_avg = accum{i} / n_avg;
    else
        j = j_last;
        a = j*1080; b = (j+1)*1080;
        mask = (tot_phi >= a) & (tot_phi < b);
        phi_cyc = tot_phi(mask) - a;
        [phi_u, ia] = unique(phi_cyc, 'stable');
        yi = varargin{i}(:);
        y_avg = interp1(phi_u, yi(mask(ia)), phi_q_int, 'pchip','extrap');
    end
    if numel(phi_q_int) ~= numel(phi_q)
        y_avg = [y_avg(:); y_avg(1)];
    end
    varargout{i} = y_avg;
end
end