function params = plot_engine(theta_crank, params, results, gear1, gear2, pair, plot_fem)
% plot_fem requires gear1, gear2 and pair to be defined!
%
%   params = plot_engine(theta_crank, params, results, gear1, gear2, pair, plot_fem)
%
%   Inputs:
%       theta_crank : crank angle [deg]
%       params : structure with geometry parameters (from parameters.m)
%       results : structure with simulation results (from engine_cycle_steady_6stroke.m)
%       gear1 : gear 1 (not used in this version)
%       gear2 : gear 2 (not used in this version)
%       pair : pair of gears (not used in this version)
%       plot_fem : plot gear mesh (default 0)
%
%   Outputs:
%       params : modified structure with updated parameters
% -------------------------------------------------------------------------
% Bart Blockmans, 2024 - bart@blockmans.net
% -------------------------------------------------------------------------

% Check plotting mode
if nargin < 3 || isempty(results) 
    mode = 'single';
else % Animation mode
    mode = 'anim';
end

% Plot gears?
if nargin > 3 && ~isempty(gear1); plot_gears = 1; else; plot_gears = 0; end
if nargin < 7; plot_fem = 0; end

if strcmp(mode, 'single')
    % Create figure with white background
    figure('Color', 'w');
    hold on; axis equal; box on;
    
    % Title
    title('Cross-section 6-stroke combustion engine')
end

% Set results structure if in single plot mode
if strcmp(mode, 'single')
    results.IV = 0; 
    results.IV_max = 1;
    results.EV = 0; 
    results.EV_max = 1;
    results.spark = 0;
    results.flame = 0;
end

% Carrier angle [rad]
theta_c = theta_crank * (pi/180) + pi/2; 

% Compute planet angle
theta_p = theta_c * params.u_pc;

% Compute engine geometry at theta_crank
[~, ~, x_piston, ~, ~, dxdth, coords] = geometry_6stroke(theta_crank, params);

% Dimensionless temperature
if strcmp(mode, 'anim')
    T = (results.T - results.T_min) / (results.T_max - results.T_min);
end

% Plot cylinder
% -------------------------------------------------------------------------

% Draw cylinder
DrawCylinder(params);

% Draw scavenge ports
DrawScavengePorts(params);

% Plot ring gear
% -------------------------------------------------------------------------

% Draw ring gear
if plot_fem
    [t, s, p] = InterpolationParameters(theta_c + 6*pi, 0, 2*pi/gear2.z, 2*pi/(gear2.z * size(gear2.splot,1)), gear2.z);
    PlotShapeEff(gear2, gear2.splot, [0;0], pair.alpha_A2, t, s, p);
    draw_gear(gear2, [0; 0], pair.alpha_A2, params.LW);
else
    if plot_gears
        draw_gear(gear2, [0; 0], pair.alpha_A2, params.LW, [], params.colors.gear);
    else
        DrawCircle([0,0], params.Do, 2 * params.Rr, params.colors.gear);
    end
end

% Plot crankshaft
% -------------------------------------------------------------------------

% Draw crankshaft
DrawCrank([coords.X_p, coords.Y_p], params);

% Plot planet gear
% -------------------------------------------------------------------------

% Draw planet gear
if plot_fem
    PlotShapeEff(gear1, gear1.splot, [coords.X_p; coords.Y_p], pair.alpha_A1 + theta_p, t, s, p);
    draw_gear(gear1, [coords.X_p; coords.Y_p], pair.alpha_A1 + theta_p, params.LW, params.colors.gear);
    DrawSimpleCircle([coords.X_p; coords.Y_p], 0.5 * params.Di, params.colors.cyl);
else
    if plot_gears
        draw_gear(gear1, [coords.X_p; coords.Y_p], pair.alpha_A1 + theta_p, params.LW, params.colors.cyl, params.colors.gear);
    else
        DrawCircle([coords.X_p; coords.Y_p], 2 * params.Rp, params.Di, params.colors.gear);
        DrawSimpleCircle([coords.X_p; coords.Y_p], 0.5 * params.Di, params.colors.cyl);
    end
end

% Draw circle
DrawSimpleCircle([coords.X_p, coords.Y_p], params.Rcra, params.colors.cw, params.LW);

% Clim
if plot_fem; clim([0 40]); end

% Plot connecting rod
% -------------------------------------------------------------------------

% Draw rod
DrawPulley([coords.X_r, coords.Y_r], [coords.X_s, coords.Y_s], params.Rcr1, params.Rcr2, params.colors.rod, params.LW);

% Draw circle
DrawSimpleCircle([coords.X_r, coords.Y_r], params.Rcra, params.colors.cyl, params.LW);

% Plot piston
% -------------------------------------------------------------------------

% Draw piston
DrawPiston([coords.X_s, coords.Y_s], params);

% Cylinder cavity & manifolds
% -------------------------------------------------------------------------
if strcmp(mode, 'anim')

    % Compute cylinder cavity
    [Xc, Yc] = CylinderCavity([coords.X_s, coords.Y_s], params);
    
    % Valves
    IV_state = (results.IV ~= 0);
    EV_state = (results.EV ~= 0);
    SP_state = (results.SP ~= 0);
    
    % Add manifolds
    [Xc, Yc] = AddManifolds(Xc, Yc, params, IV_state, EV_state, SP_state);
    
    % Color cylinder cavity according to temperature or pressure
    ColorCavity(Xc, Yc, T, params.B/2 + params.tw);

end

% Draw particles
% -------------------------------------------------------------------------
if params.particles.anim && strcmp(mode, 'anim')

    % Did overall valve state change?
    params.valve_change = 0;
    if IV_state ~= params.memory.IV_state; params.valve_change = 1; end
    if EV_state ~= params.memory.EV_state; params.valve_change = 1; end
    if SP_state ~= params.memory.SP_state; params.valve_change = 1; end

    % Store valve state in memory
    params.memory.IV_state = IV_state;
    params.memory.EV_state = EV_state;
    params.memory.SP_state = SP_state;

    % Flame change
    dflame = max(0, results.flame - params.memory.flame);

    % Store flame
    params.memory.flame = results.flame;

    % Get particle positions & velocities at previous animation time step
    Xp = params.memory.Xp;
    Yp = params.memory.Yp;
    Up = params.memory.Up;
    Vp = params.memory.Vp;

    % Get particle colors
    Cp = params.memory.Cp;

    % Get piston Y position
    Ypiston = params.Ys_max - x_piston + 2 * params.H / 3;

    % Get piston velocity
    Vpiston = -dxdth * (params.N_rpm * pi / 30) * 180/pi * params.speed_scaling; % [mm/s]

    % Update particle positions & velocities
    [Xp, Yp, Up, Vp, Cp, params] = particle_simulation(Xp, Yp, Up, Vp, Cp, ...
                                Xc, Yc, Ypiston, Vpiston, dflame, params);

    % Set burn ratio of particles
    Cp = burned_particles(Xp, Yp, Cp, Xc, Yc, results.flame);

    % Store in memory
    params.memory.Xp = Xp;
    params.memory.Yp = Yp;
    params.memory.Up = Up;
    params.memory.Vp = Vp;
    params.memory.Cp = Cp;

    % Draw particles
    DrawParticles(Xp, Yp, Cp, params.particles.Rp, results.T, ...
                results.T_min, results.T_max, params.B/2 + params.tw);

end

% Draw valves
% -------------------------------------------------------------------------

% Draw valves
DrawValves(results.IV/results.IV_max, results.EV/results.EV_max, params);

% Draw cams
if strcmp(mode, 'anim')
    DrawCams(theta_crank, params);
end

% Draw bogie
% -------------------------------------------------------------------------

% Draw bougie
if strcmp(mode, 'anim')
    DrawBougie(params, results.spark, results.flame, Xc, Yc);
else
    DrawBougie(params);
end

% Axis limits & labels
% -------------------------------------------------------------------------

% Set axis limits
if ~isfield(params, 'xlimits')
    axis tight;
    xlimits = xlim;
    ylimits = ylim;
    params.xlimits = [round(10*(xlimits(1)-params.B/2))/10, round(10*(xlimits(2) + params.B/2))/10];
    params.ylimits = [round(10*(ylimits(1)-params.B/2))/10, round(10*(ylimits(2) + params.B/2))/10];
end
xlim(params.xlimits);
ylim(params.ylimits);

if strcmp(mode,'single')
    
    % Labels
    xlabel('[mm]')
    ylabel('[mm]')

else

    % Turn off axes
    axis off;

end

