function [Up, Vp] = SetManifoldVelocities(Xp, Yp, Up, Vp, IV, EV, SP, Y_piston, params)
%SET_MANIFOLD_VELOCITIES  Set the velocities of the particles in the intake, exhaust, and scavenge manifolds.
%
%   [Up, Vp] = SetManifoldVelocities(Xp, Yp, Up, Vp, IV, EV, SP, Y_piston, params)
%
%   Inputs:
%       Xp, Yp : particle positions [mm]
%       Up, Vp : particle velocities [mm/s]
%       IV, EV, SP : valve states (1=open, 0=closed)
%       Y_piston : piston Y position [mm]
%       params : structure with geometry parameters (from parameters.m)
%
%   Outputs:
%       Up, Vp : updated particle velocities [mm/s]
% -------------------------------------------------------------------------
% Bart Blockmans, 2024 - bart@blockmans.net
% -------------------------------------------------------------------------

% Intake manifold
% -------------------------------------------------------------------------
if IV

% Intake manifold geometry
X_intake = [params.xmi1, params.xir, fliplr(params.xmi2), params.xmi1(1)];
Y_intake = [params.ymi1, params.yir, fliplr(params.ymi2), params.ymi1(1)];

% Find indices of particles in intake manifold
[in, ~] = inpolygon(Xp, Yp, X_intake, Y_intake);

% Top & Bottom boundary
X_intake_top = [params.xmi2(1:end-1), fliplr(params.xir(3:4))];
Y_intake_top = [params.ymi2(1:end-1), fliplr(params.yir(3:4))];
X_intake_bot = [params.xmi1(1:end-1), params.xir(1:2)];
Y_intake_bot = [params.ymi1(1:end-1), params.yir(1:2)];

% Loop over all particles in intake manifold
for i = 1:numel(Xp)
    if in(i)
        % Determine velocity unit vector
        x = Xp(i);
        % Intersection of top/bottom X-ranges
        xmin = max(min(X_intake_top), min(X_intake_bot));
        xmax = min(max(X_intake_top), max(X_intake_bot));
        if xmax <= xmin
            % Fallback: keep existing velocity if geometry is degenerate
            continue
        end
        dx = max(1e-6, 1e-3 * (xmax - xmin));
        % Clamp evaluation point to avoid extrapolation artifacts
        x0 = min(max(x, xmin + dx), xmax - dx);
        % Centerline slope from averaged top/bottom derivatives
        yt_m = interp1(X_intake_top, Y_intake_top, x0 - dx, 'linear', 'extrap');
        yt_p = interp1(X_intake_top, Y_intake_top, x0 + dx, 'linear', 'extrap');
        yb_m = interp1(X_intake_bot, Y_intake_bot, x0 - dx, 'linear', 'extrap');
        yb_p = interp1(X_intake_bot, Y_intake_bot, x0 + dx, 'linear', 'extrap');
        dyt = (yt_p - yt_m) / (2*dx);
        dyb = (yb_p - yb_m) / (2*dx);
        slope = 0.5 * (dyt + dyb);
        % Unit vector along duct (positive X direction)
        nx = 1;
        ny = slope;
        inv_norm = 1 / hypot(nx, ny);
        Up(i) = nx * inv_norm;
        Vp(i) = ny * inv_norm;
        % Set velocity magnitude with small random variation
        Up(i) = 1e3 * params.particles.speedRms * abs(1 + 0.1 * randn) * Up(i);
        Vp(i) = 1e3 * params.particles.speedRms * abs(1 + 0.1 * randn) * Vp(i);
    end
end
end

% Exhaust manifold
% -------------------------------------------------------------------------
if EV

% Exhaust manifold geometry
X_exhaust = fliplr([params.xme2, params.xer, fliplr(params.xme1), params.xme2(1)]);
Y_exhaust = fliplr([params.yme2, params.yer, fliplr(params.yme1), params.yme2(2)]);

% Find indices of particles in exhaust manifold
[in, ~] = inpolygon(Xp, Yp, X_exhaust, Y_exhaust);

% Top & Bottom boundary
X_exhaust_top = [params.xme1(1:end-1), fliplr(params.xer(3:4))];
Y_exhaust_top = [params.yme1(1:end-1), fliplr(params.yer(3:4))];
X_exhaust_bot = [params.xme2(1:end-1), params.xer(1:2)];
Y_exhaust_bot = [params.yme2(1:end-1), params.yer(1:2)];

% Loop over all particles in exhaust manifold
for i = 1:numel(Xp)
    if in(i)
        % Determine velocity unit vector
        x = Xp(i);
        % Intersection of top/bottom X-ranges
        xmin = max(min(X_exhaust_top), min(X_exhaust_bot));
        xmax = min(max(X_exhaust_top), max(X_exhaust_bot));
        if xmax <= xmin
            % Fallback: keep existing velocity if geometry is degenerate
            continue
        end
        dx = max(1e-6, 1e-3 * (xmax - xmin));
        % Clamp evaluation point to avoid extrapolation artifacts
        x0 = min(max(x, xmin + dx), xmax - dx);
        % Centerline slope from averaged top/bottom derivatives
        yt_m = interp1(X_exhaust_top, Y_exhaust_top, x0 - dx, 'linear', 'extrap');
        yt_p = interp1(X_exhaust_top, Y_exhaust_top, x0 + dx, 'linear', 'extrap');
        yb_m = interp1(X_exhaust_bot, Y_exhaust_bot, x0 - dx, 'linear', 'extrap');
        yb_p = interp1(X_exhaust_bot, Y_exhaust_bot, x0 + dx, 'linear', 'extrap');
        dyt = (yt_p - yt_m) / (2*dx);
        dyb = (yb_p - yb_m) / (2*dx);
        slope = 0.5 * (dyt + dyb);
        % Unit vector along duct (positive X direction)
        nx = 1;
        ny = slope;
        inv_norm = 1 / hypot(nx, ny);
        Up(i) = nx * inv_norm;
        Vp(i) = ny * inv_norm;
        % Set velocity magnitude with small random variation
        Up(i) = 1e3 * params.particles.speedRms * abs(1 + 0.1 * randn) * Up(i);
        Vp(i) = 1e3 * params.particles.speedRms * abs(1 + 0.1 * randn) * Vp(i);
    end
end
end

% Left scavenge port
% -------------------------------------------------------------------------
if SP

% Scavenge reservoir top & bottom
Y_res_bot = params.BDC2_y + 2*params.H/3;
Y_res_top = params.BDC1_y + 2*params.H/3;

% Top scavenge port entry
Y_SP_top = Y_res_top;

% Bottom scavenge port entry
Y_SP_bot = Y_piston;

% X position scavenge port entry
Xref = -0.5 * params.B;

% Left scavenge port geometry
X_scav_l = [params.xmsl1, params.xsrl, fliplr(params.xmsl2), params.xmsl1(1)];
Y_scav_l = [params.ymsl1, params.ysrl, fliplr(params.ymsl2), params.ymsl1(1)];

% Find indices of particles in left scavenge port
[in, ~] = inpolygon(Xp, Yp, X_scav_l, Y_scav_l);

% Loop over all particles in left scavenge port
for i = 1:numel(Xp)
    if in(i)
        % Relative Y position particle between top & bottom of reservoir
        Yr_i = (Yp(i) - Y_res_bot)/(Y_res_top - Y_res_bot);
        % Reference port point
        Yref = Y_SP_bot + Yr_i * (Y_SP_top - Y_SP_bot);
        % Unit vector pointing from particle to reference port point
        nx = Xref - Xp(i);
        ny = Yref - Yp(i);
        inv_norm = 1 / hypot(nx, ny);
        Up(i) = nx * inv_norm;
        Vp(i) = ny * inv_norm;
        % Set velocity magnitude with small random variation
        Up(i) = 1e3 * params.particles.speedRms * abs(1 + 0.1 * randn) * Up(i);
        Vp(i) = 1e3 * params.particles.speedRms * abs(1 + 0.1 * randn) * Vp(i);
    end
end
end

% Right scavenge port
% -------------------------------------------------------------------------
if SP

% Scavenge reservoir top & bottom
Y_res_bot = params.BDC2_y + 2*params.H/3;
Y_res_top = params.BDC1_y + 2*params.H/3;

% Top scavenge port entry
Y_SP_top = Y_res_top;

% Bottom scavenge port entry
Y_SP_bot = Y_piston;

% X position scavenge port entry
Xref = 0.5 * params.B;

% Right scavenge port geometry
X_scav_r = fliplr([params.xmsr1, params.xsrr, fliplr(params.xmsr2), params.xmsr1(1)]);
Y_scav_r = fliplr([params.ymsr1, params.ysrr, fliplr(params.ymsr2), params.ymsr1(1)]);

% Find indices of particles in left scavenge port
[in, ~] = inpolygon(Xp, Yp, X_scav_r, Y_scav_r);

% Loop over all particles in left scavenge port
for i = 1:numel(Xp)
    if in(i)
        % Relative Y position particle between top & bottom of reservoir
        Yr_i = (Yp(i) - Y_res_bot)/(Y_res_top - Y_res_bot);
        % Reference port point
        Yref = Y_SP_bot + Yr_i * (Y_SP_top - Y_SP_bot);
        % Unit vector pointing from particle to reference port point
        nx = Xref - Xp(i);
        ny = Yref - Yp(i);
        inv_norm = 1 / hypot(nx, ny);
        Up(i) = nx * inv_norm;
        Vp(i) = ny * inv_norm;
        % Set velocity magnitude with small random variation
        Up(i) = 1e3 * params.particles.speedRms * abs(1 + 0.1 * randn) * Up(i);
        Vp(i) = 1e3 * params.particles.speedRms * abs(1 + 0.1 * randn) * Vp(i);
    end
end
end