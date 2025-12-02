function [] = plot_results_6stroke(result, params)
%PLOT_RESULTS_6STROKE  Compact visualization of a converged 6-stroke cycle.
% Usage:
%   plot_results_6stroke(result)                % minimal (no params summary)
%   plot_results_6stroke(result, params)        % prints summary using params
%
% Figures:
%   1) Engine geometry (V, x, dV/dphi)
%   2) Cylinder thermodynamics (p, T, F_gas, torque)
%   3) Torque output (single, cumulative work, flat-6 sum)
%   4) Driveline (gearbox input, DMF output)
%   5) Further thermodynamics (m_cyl, mO2, dQ/dphi, Qdot_wall)
%   6) Combustion & power (burn fractions, heat release, wall losses)
%   7) Valve and port timing (areas and lifts)
%   8) Transmission Error dynamics (optional, DMF-embedded method)
%   9) Instantaneous IMEP & Power vs Time
%   10) Piston kinematics (x, dx/dphi, d²x/dphi²)
% -------------------------------------------------------------------------
% Bart Blockmans, 2024 - bart@blockmans.net
% -------------------------------------------------------------------------


phi   = result.phi_deg(:);
pcyl  = result.p_cyl(:);
Tgas  = result.T_cyl(:);
Tq    = result.T_ind(:);

% Optional fields (compute light fallbacks if missing)
V_m3      = getfield_or(result, 'V_m3',       []);
dVdphi    = getfield_or(result, 'dVdphi_m3_per_deg', []);
x_mm      = getfield_or(result, 'x_mm',       []);
dx_dphi_m = getfield_or(result, 'dx_dphi_m_per_deg', []);
ddx_dphi2 = getfield_or(result, 'ddx_dphi2_m_per_deg2', []);
F_gas     = getfield_or(result, 'F_gas_N',    []);
W_cum     = getfield_or(result, 'W_cum_J',    []);
dQ_dphi   = getfield_or(result, 'dQ_dphi_J_per_deg', []);
Qdot_wall = getfield_or(result, 'Qdot_wall_W', []);
A_int     = getfield_or(result, 'A_int_m2',   []);
A_exh     = getfield_or(result, 'A_exh_m2',   []);
A_scav    = getfield_or(result, 'A_scav_m2',  []);
xburn1    = getfield_or(result, 'xburn1',     []);
xburn2    = getfield_or(result, 'xburn2',     []);
t_sec     = getfield_or(result, 't_sec',      []);
power_W   = getfield_or(result, 'power_W',    []);
imep_inst = getfield_or(result, 'imep_inst_Pa', []);
L_int     = getfield_or(result, 'L_int_m',    []);
L_exh     = getfield_or(result, 'L_exh_m',    []);
L_scav    = getfield_or(result, 'L_scav_m',   []);

% -------- 1) Geometry ------------------------------------------------------
fig = figure('Name','Engine Geometry','Color','w');
tiledlayout(fig,3,1,'TileSpacing','compact','Padding','compact');

% Volume
nexttile;
if ~isempty(V_m3), plot(phi, V_m3*1e3, 'LineWidth',1.5); end  % [L]
decorate_axis('Volume', '[L]');

% Displacement
nexttile;
if ~isempty(x_mm), plot(phi, x_mm, 'LineWidth',1.5); end      % [mm]
decorate_axis('Piston displacement from highest TDC', '[mm]');

% dV/dphi
nexttile;
if ~isempty(dVdphi), plot(phi, dVdphi*1e6, 'LineWidth',1.5); end % [cm^3/deg]
decorate_axis('dV/d\phi', '[cm^3/deg]');

% -------- 2) Cylinder thermodynamics --------------------------------------
fig = figure('Name','Cylinder Thermodynamics','Color','w');
tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');

% Pressure
nexttile;
plot(phi, pcyl*1e-5, 'LineWidth',1.5);         % [bar]
decorate_axis('Cylinder pressure', '[bar]');

% Temperature
nexttile;
plot(phi, Tgas, 'LineWidth',1.5);              % [K]
decorate_axis('Gas temperature', '[K]');

% Force
nexttile;
if ~isempty(F_gas), plot(phi, F_gas, 'LineWidth',1.5); end    % [N]
decorate_axis('Gas force on piston', '[N]');

% Torque
nexttile;
plot(phi, Tq, 'LineWidth',1.5);                % [Nm]
decorate_axis('Single-cylinder gas torque', '[Nm]');

% -------- 3) Torque output -------------------------------------------------
fig = figure('Name','Torque Output','Color','w');
tiledlayout(fig,3,1,'TileSpacing','compact','Padding','compact');

% Single-cylinder torque
nexttile;
plot(phi, Tq, 'LineWidth',1.5);
decorate_axis('Single-cylinder gas torque', '[Nm]');

% Cumulative work
nexttile;
W_cum_local = W_cum;
if isempty(W_cum_local)
    % fallback if not provided
    if isempty(dVdphi)
        W_cum_local = NaN(size(phi));
    else
        W_cum_local = cumtrapz(phi, pcyl .* dVdphi);
    end
end
plot(phi, W_cum_local, 'LineWidth',1.5);
decorate_axis('Cumulative work \int p dV', '[J]');

% Flat-six (6 cylinders, 180° phase offsets)
nexttile;
Tsingle = Tq(:);
if ~isempty(Tsingle)
    step = mean(diff(phi));
    phase_deg = 180;
    shift_idx = max(1, round(phase_deg/step));
    colors = lines(6);
    hold on;
    Tsum = zeros(size(Tsingle));
    for k = 1:6
        Tk = circshift(Tsingle, (k-1)*shift_idx);
        plot(phi, Tk, 'Color',colors(k,:), 'LineWidth',1.0, 'HandleVisibility','off');
        Tsum = Tsum + Tk;
    end
    plot(phi, Tsum, 'k-', 'LineWidth',2.0);
    hold off;
    % Add Y-axis margin
    yl = ylim;
    ylim([yl(1) - 0.1*(yl(2)-yl(1)), yl(2) + 0.1*(yl(2)-yl(1))]);
end
decorate_axis('Flat-six summed torque (6x, 180° phased)', '[Nm]');

% -------- 4) Driveline (gearbox input) -------------------------------------
fig = figure('Name','Gearbox Input (DMF Output)','Color','w');
tiledlayout(fig,2,1,'TileSpacing','compact','Padding','compact');

% Gearbox torque vs crank angle
nexttile;
hold on;
if isfield(result,'T_gbx_in') && ~isempty(result.T_gbx_in)
    plot(phi, result.T_gbx_in, 'LineWidth',1.8, 'DisplayName','T_{gbx,in}');
end
if isfield(result,'T_engine_sum') && ~isempty(result.T_engine_sum)
    plot(phi, result.T_engine_sum, '--', 'LineWidth',1.0, 'DisplayName','Flat-6 sum (pre-DMF)');
end
hold off;
legend('Location','northeast','AutoUpdate','off');
decorate_axis('Gearbox input torque (DMF output)', '[Nm]');

% Gearbox power vs crank angle
nexttile;
if isfield(result,'P_gbx_in') && ~isempty(result.P_gbx_in)
    plot(phi, result.P_gbx_in/1e3, 'LineWidth',1.8);
end
decorate_axis('Gearbox input power', '[kW]');

% -------- 5) Further thermodynamics ----------------------------------------
fig = figure('Name','Further Thermodynamics','Color','w');
tiledlayout(fig,2,1,'TileSpacing','compact','Padding','compact');

% Trapped mass
nexttile;
plot(phi, result.m_cyl, 'LineWidth',1.5);      % [kg]
decorate_axis('Trapped mass m_{cyl}', '[kg]');

% O2 mass
nexttile;
if isfield(result,'mO2') && numel(result.mO2)==numel(phi)
    plot(phi, result.mO2, 'LineWidth',1.5);
    decorate_axis('O_2 mass in-cylinder', '[kg]');
else
    axis off;
end

% -------- 6) Combustion & power --------------------------------------------
fig = figure('Name','Combustion & Power','Color','w');
tiledlayout(fig,3,1,'TileSpacing','compact','Padding','compact');

% Burn fractions
nexttile;
hold on;
if ~isempty(xburn1), plot(phi, xburn1, 'LineWidth',1.5, 'DisplayName','Burn #1'); end
if ~isempty(xburn2), plot(phi, xburn2, 'LineWidth',1.5, 'DisplayName','Burn #2'); end
hold off;
legend('Location','southeast', 'AutoUpdate','off');
decorate_axis('Wiebe burn fraction', '[-]');

% Heat release dQ/dphi
nexttile;
if ~isempty(dQ_dphi), plot(phi, dQ_dphi, 'LineWidth',1.5); end
decorate_axis('Heat release dQ/d\phi', '[J/deg]');

% Wall heat loss
nexttile;
if ~isempty(Qdot_wall), plot(phi, Qdot_wall, 'LineWidth',1.5); end
decorate_axis('Wall heat loss rate Q_{wall}', '[W]');

% -------- 7) Valve and port timing -----------------------------------------
fig = figure('Name','Valve and Port Timing','Color','w');
tiledlayout(fig,2,1,'TileSpacing','compact','Padding','compact');

% Areas
nexttile;
hold on;
if ~isempty(A_int),  plot(phi, A_int,  'LineWidth',1.5, 'DisplayName','Intake area'); end
if ~isempty(A_exh),  plot(phi, A_exh,  'LineWidth',1.5, 'DisplayName','Exhaust area'); end
if ~isempty(A_scav), plot(phi, A_scav, 'LineWidth',1.5, 'DisplayName','Scavenge ports'); end
hold off;
legend('Location','northeastoutside', 'AutoUpdate','off');
decorate_axis('Valve/port effective areas', '[m^2]');

% Effective lifts
nexttile;
hold on;
if ~isempty(L_int),  plot(phi, L_int*1e3,  'LineWidth',1.5, 'DisplayName','Intake lift'); end
if ~isempty(L_exh),  plot(phi, L_exh*1e3,  'LineWidth',1.5, 'DisplayName','Exhaust lift'); end
if ~isempty(L_scav), plot(phi, L_scav*1e3, 'LineWidth',1.5, 'DisplayName','Scavenge height'); end
hold off;
legend('Location','northeastoutside', 'AutoUpdate','off');
decorate_axis('Effective lifts', '[mm]');

% -------- 8) Transmission Error dynamics (optional, DMF-embedded) ----------
if isfield(result,'te_dyn_avg') && isfield(result.te_dyn_avg,'phi_deg')
    phi_te = result.te_dyn_avg.phi_deg(:);
    TE_deg_avg   = getfield_or(result.te_dyn_avg, 'TE_deg',   []);
    kmesh_avg    = getfield_or(result.te_dyn_avg, 'kmesh',   []);
    Tmesh_avg    = getfield_or(result.te_dyn_avg, 'T_mesh',  []);
    Teng_avg     = getfield_or(result.te_dyn_avg, 'T_eng',   []);

    fig = figure('Name','Transmission Error Dynamics (Averaged)','Color','w');
    tiledlayout(fig,3,1,'TileSpacing','compact','Padding','compact');

    % Transmission error across mesh (planet gear minus crank)
    nexttile;
    if ~isempty(TE_deg_avg)
        plot(phi_te, TE_deg_avg, 'LineWidth',1.5);
    end
    decorate_axis('Transmission error (\theta_{pg}-\theta_{e})', '[deg]');

    % Mesh stiffness vs crank angle
    nexttile;
    if ~isempty(kmesh_avg)
        plot(phi_te, kmesh_avg, 'LineWidth',1.5);
    end
    decorate_axis('Gear mesh stiffness', '[N·m/rad]');

    % Mesh torque (and optionally engine torque) vs crank angle
    nexttile;
    hold on;
    if ~isempty(Tmesh_avg)
        plot(phi_te, Tmesh_avg, 'LineWidth',1.5, 'DisplayName','T_{mesh}');
    end
    if ~isempty(Teng_avg)
        plot(phi_te, Teng_avg, '--', 'LineWidth',1.2, 'DisplayName','T_{eng}');
    end
    hold off;
    if ~isempty(Tmesh_avg) || ~isempty(Teng_avg)
        legend('Location','northeast','AutoUpdate','off');
    end
    decorate_axis('Torques', '[N·m]');
end

% -------- 9) Instantaneous IMEP & Power vs Time --------------------------
fig = figure('Name','IMEP & Power vs Time','Color','w');
tiledlayout(fig,2,1,'TileSpacing','compact','Padding','compact');

% Instantaneous IMEP
nexttile;
if ~isempty(imep_inst), plot(phi, imep_inst*1e-5, 'LineWidth',1.5); end
decorate_axis('Instantaneous IMEP', '[bar]');

% Mechanical power vs time
nexttile;
if ~isempty(power_W) && ~isempty(t_sec)
    plot(t_sec*1e3, power_W/1e3, 'LineWidth',1.5);
    grid on; box on;
    xlabel('t [ms]'); ylabel('Power [kW]');
    title('Mechanical power vs time');
else
    axis off;
end

% -------- 10) Piston kinematics ---------------------------------------------
fig = figure('Name','Piston Kinematics','Color','w');
tiledlayout(fig,3,1,'TileSpacing','compact','Padding','compact');

% Displacement
nexttile;
if ~isempty(x_mm), plot(phi, x_mm, 'LineWidth',1.5); end
decorate_axis('Displacement x(\phi)', '[mm]');

% Velocity
nexttile;
if ~isempty(dx_dphi_m), plot(phi, dx_dphi_m*1e3, 'LineWidth',1.5); end % [mm/deg]
decorate_axis('Velocity dx/d\phi', '[mm/deg]');

% Acceleration
nexttile;
if ~isempty(ddx_dphi2), plot(phi, ddx_dphi2*1e3, 'LineWidth',1.5); end % [mm/deg^2]
decorate_axis('Acceleration d^2x/d\phi^2', '[mm/deg^2]');

% -------- Summary (optional) ----------------------------------------------
if nargin > 1 && ~isempty(params)
fprintf('--- Converged 6-stroke cycle summary ---\n');
    if isfield(result,'N_rpm')
        fprintf('Engine speed: %.0f rpm\n', result.N_rpm);
    else
fprintf('Engine speed: %.0f rpm\n', params.N_rpm);
    end
fprintf('Indicated IMEP (vs displacement): %.2f bar\n', result.IMEP_bar);
fprintf('Indicated mean torque (1 cyl): %.2f Nm\n', result.T_ind_mean);
fprintf('Fuel per 1080° cycle: %.4e kg\n', result.m_fuel_cycle_total);
fprintf('Trapped mass range: [%.4e, %.4e] kg\n', min(result.m_cyl), max(result.m_cyl));
fprintf('O2 mass at SOC2: %.4e kg\n', result.mO2_at_SOC2);
end

end

% ==== Helpers ==============================================================
function decorate_axis(titleStr, yLabelStr)
strokeBounds  = 0:180:1080;
grid on; box on;
xlim([0 1080]); xticks(strokeBounds);
for xb = strokeBounds, xline(xb,'--','Color',[0.85 0.85 0.85]); end
xlabel('\phi [deg CA]'); ylabel(yLabelStr);
title(titleStr);
end

function val = getfield_or(s, name, defaultVal)
if isfield(s, name), val = s.(name); else, val = defaultVal; end
end