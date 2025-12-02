function params = extrema_6stroke(params, generate_figs)
%EXTREMA_6STROKE  Compute TDC/BDC angles and positions for 6-stroke geometry.
%
%   params = extrema_6stroke(params, generate_figs)
%
%   Inputs:
%       params      : structure with geometry parameters (from parameters.m)
%       generate_figs : optional flag (default 0) to generate figures
%
%   Outputs:
%       params      : modified structure with added fields:
%                     - TDC1_th, TDC2_th: TDC angles [deg] (arrays)
%                     - BDC1_th, BDC2_th: BDC angles [deg] (arrays)
%                     - TDC1_x, TDC2_x: TDC x positions [mm] (arrays)
%                     - BDC1_x, BDC2_x: BDC x positions [mm] (arrays)
%                     - TDC1_y, TDC2_y: TDC absolute Y positions [mm] (arrays)
%                     - BDC1_y, BDC2_y: BDC absolute Y positions [mm] (arrays)
%                     - Ys_max: maximum slider Y position [mm]
%                     - theta_TDC1: angle of highest TDC [deg]
%                     - All seven dead centers: DC_th, DC_x, DC_names
% -------------------------------------------------------------------------
% Bart Blockmans, 2024 - bart@blockmans.net
% -------------------------------------------------------------------------


if nargin < 2 || isempty(generate_figs)
    generate_figs = 0;
end

% Extract geometry parameters from params structure
Rr    = params.Rr;      % ring pitch radius [mm]
Rp    = params.Rp;      % planet pitch radius [mm]
delta = params.delta;   % offset of tracked point on planet [mm]
L     = params.L;       % rod length to Y-axis slider [mm]
psi0  = params.psi0;    % body-fixed angle of tracked point [rad]

% Compute highest TDC
[Ys_max, theta_TDC1] = HighestTDC(Rr, Rp, delta, L, psi0);
if generate_figs
    fprintf('Highest TDC at θ ≈ %.4f°: Y_s,max = %.6f mm\n', theta_TDC1, Ys_max);
end

% Store in params
params.Ys_max = Ys_max;

% Crank angle grid (0..6*pi rad == 0..1080 deg)
np        = 1000;                      % discretization points
theta_rad = linspace(0, 6*pi, np);     % [rad]
theta_deg = theta_rad * 180/pi;        % [deg] for geometry_6stroke

% Common tick positions (radians) and labels (degrees)
tick_deg  = [0 180 360 540 720 900 1080];
tick_rad  = tick_deg * pi/180;
tick_lbls = compose('%d', tick_deg);

% Evaluate geometry (vectorized)
[V, A, x_piston, dVdth, dAdth, dxdth, coords] = ...
    geometry_6stroke(theta_deg, params);

% Compute max and minimum volume and store in params
params.V_max = max(V);
params.V_min = min(V);

% Generate figures if requested
if generate_figs
    % 1) x, A, V vs crank angle
    figure('Name','x / A / V vs crank angle','Color','w');
    subplot(3,1,1);
    plot(theta_rad, x_piston, 'LineWidth',1.5); grid on;
    xlim([0 6*pi]); xticks(tick_rad); xticklabels(tick_lbls);
    xlabel('\theta (deg)'); ylabel('x(\theta) [mm]');
    title('Piston displacement from highest TDC');
    
    subplot(3,1,2);
    plot(theta_rad, A, 'LineWidth',1.5); grid on;
    xlim([0 6*pi]); xticks(tick_rad); xticklabels(tick_lbls);
    xlabel('\theta (deg)'); ylabel('A(\theta) [mm^2]');
    title('Wetted surface area');
    
    subplot(3,1,3);
    plot(theta_rad, V, 'LineWidth',1.5); grid on;
    xlim([0 6*pi]); xticks(tick_rad); xticklabels(tick_lbls);
    xlabel('\theta (deg)'); ylabel('V(\theta) [mm^3]');
    title('Cylinder volume');
    
    % 2) dx/dθ, dA/dθ, dV/dθ vs crank angle
    figure('Name','derivatives vs crank angle','Color','w');
    subplot(3,1,1);
    plot(theta_rad, dxdth, 'LineWidth',1.5); grid on;
    xlim([0 6*pi]); xticks(tick_rad); xticklabels(tick_lbls);
    xlabel('\theta (deg)'); ylabel('dx/d\theta [mm/deg]');
    title('Piston speed (per deg)');
    
    subplot(3,1,2);
    plot(theta_rad, dAdth, 'LineWidth',1.5); grid on;
    xlim([0 6*pi]); xticks(tick_rad); xticklabels(tick_lbls);
    xlabel('\theta (deg)'); ylabel('dA/d\theta [mm^2/deg]');
    title('dA/d\theta');
    
    subplot(3,1,3);
    plot(theta_rad, dVdth, 'LineWidth',1.5); grid on;
    xlim([0 6*pi]); xticks(tick_rad); xticklabels(tick_lbls);
    xlabel('\theta (deg)'); ylabel('dV/d\theta [mm^3/deg]');
    title('dV/d\theta');
    
    % 3) XY trace of planet center, connection point, and slider
    figure('Name','XY traces','Color','w'); hold on; axis equal; grid on;
    plot(coords.X_p, coords.Y_p, 'k-', 'LineWidth',1.3, 'DisplayName','Planet center');
    plot(coords.X_r, coords.Y_r, 'r-', 'LineWidth',1.3, 'DisplayName','Connection point');
    plot(coords.X_s, coords.Y_s, 'b-', 'LineWidth',1.3, 'DisplayName','Slider (Y-axis)');
    xlabel('X [mm]'); ylabel('Y [mm]');
    title('Absolute XY traces (ring center at origin)');
    legend('Location','best');
    
    % 4) X(t), Y(t) of planet center, connection point, slider
    figure('Name','Coordinates vs crank angle','Color','w');
    
    subplot(2,1,1); % X vs theta
    plot(theta_rad, coords.X_p, 'k-', 'LineWidth',1.2, 'DisplayName','X_p');
    hold on;
    plot(theta_rad, coords.X_r, 'r-', 'LineWidth',1.2, 'DisplayName','X_r');
    plot(theta_rad, coords.X_s, 'b-', 'LineWidth',1.2, 'DisplayName','X_s');
    grid on; xlabel('\theta (deg)'); ylabel('X [mm]');
    xlim([0 6*pi]); xticks(tick_rad); xticklabels(tick_lbls);
    title('X-coordinates vs crank angle'); legend('Location','best');
    
    subplot(2,1,2); % Y vs theta
    plot(theta_rad, coords.Y_p, 'k-', 'LineWidth',1.2, 'DisplayName','Y_p');
    hold on;
    plot(theta_rad, coords.Y_r, 'r-', 'LineWidth',1.2, 'DisplayName','Y_r');
    plot(theta_rad, coords.Y_s, 'b-', 'LineWidth',1.2, 'DisplayName','Y_s');
    grid on; xlabel('\theta (deg)'); ylabel('Y [mm]');
    xlim([0 6*pi]); xticks(tick_rad); xticklabels(tick_lbls);
    title('Y-coordinates vs crank angle'); legend('Location','best');
    
    % 5) Y(t) of slider
    figure('Name','Y coordinate piston vs crank angle','Color','w');
    plot(theta_rad, coords.Y_s, 'k-', 'LineWidth',1.2, 'DisplayName','Y_p');
    grid on; xlabel('\theta (deg)'); ylabel('Y [mm]');
    xlim([0 6*pi]); xticks(tick_rad); xticklabels(tick_lbls);
    title('Y-coordinate piston vs crank angle');
    
    % Console check
    fprintf('x(0 deg)   = %.3f mm,  x(1080 deg) = %.3f mm\n', x_piston(1), x_piston(end));
    fprintf('V(0 deg)   = %.1f mm^3, V(1080 deg)= %.1f mm^3\n', V(1), V(end));
end

% Dead-center detection, classification, and printout

% Smooth interpolants for refinement
dx_fun = @(th) interp1(theta_deg, dxdth,   th, 'pchip');
x_fun  = @(th) interp1(theta_deg, x_piston, th, 'pchip');

% Find sign changes (exclude NaNs)
sgn = sign(dxdth);
valid_pair = ~isnan(sgn(1:end-1)) & ~isnan(sgn(2:end));
sign_change = valid_pair & (sgn(1:end-1).*sgn(2:end) < 0);
brk_idx = find(sign_change);

% Include endpoints if derivative is near zero there
end_candidates = [];
if isfinite(dx_fun(0))    && abs(dx_fun(0))    < 1e-10, end_candidates(end+1) = 0;     end
if isfinite(dx_fun(540))  && abs(dx_fun(540))  < 1e-10, end_candidates(end+1) = 540;   end
if isfinite(dx_fun(1080)) && abs(dx_fun(1080)) < 1e-10, end_candidates(end+1) = 1080;  end

% Refine each bracket with fzero
roots_deg = [];
for k = 1:numel(brk_idx)
    a = theta_deg(brk_idx(k));
    b = theta_deg(brk_idx(k)+1);
    try
        th_root = fzero(dx_fun, [a b]);
    catch
        th_root = 0.5*(a+b);  % fallback
    end
    roots_deg(end+1) = th_root; %#ok<AGROW>
end
% Add endpoint roots and tidy up
roots_deg = sort([roots_deg, end_candidates]);
if ~isempty(roots_deg)
    keep = [true, diff(roots_deg) > 1e-6];  % merge near-duplicates
    roots_deg = roots_deg(keep);
end

% Classify each root as TDC (min) or BDC (max)
kind = strings(size(roots_deg));   % "TDC" or "BDC"
eps_deg = 1e-4;
for i = 1:numel(roots_deg)
    th0 = roots_deg(i);
    xm  = x_fun(th0 - eps_deg);
    x0  = x_fun(th0);
    xp  = x_fun(th0 + eps_deg);
    if xp > x0 && xm > x0
        kind(i) = "TDC";
    elseif xp < x0 && xm < x0
        kind(i) = "BDC";
    else
        % fallback by curvature sign at the root
        kind(i) = "TDC";
        if xp - 2*x0 + xm < 0, kind(i) = "BDC"; end
    end
end

% Map each root to the nearest of the 7 canonical angles for labeling
targets = [  0, 180, 360, 540, 720, 900, 1080 ];
names   = ["TDC2","BDC1","TDC1","BDC2","TDC1","BDC1","TDC2"];

picked_th = nan(1,numel(targets));
picked_x  = nan(1,numel(targets));
picked_nm = strings(1,numel(targets));
picked_ix = nan(1,numel(targets));

available = true(size(roots_deg));
for j = 1:numel(targets)
    [~, iroot] = min( abs(roots_deg - targets(j)) + (~available)*1e6 );
    picked_ix(j) = iroot;
    picked_th(j) = roots_deg(iroot);
    picked_x(j)  = x_fun(picked_th(j));
    picked_nm(j) = names(j);
    available(iroot) = false;
end

% Extract requested DC values
th_TDC1 = picked_th(picked_nm=="TDC1");
x_TDC1  = picked_x(picked_nm=="TDC1");

th_TDC2 = picked_th(picked_nm=="TDC2");
x_TDC2  = picked_x(picked_nm=="TDC2");

th_BDC1 = picked_th(picked_nm=="BDC1");
x_BDC1  = picked_x(picked_nm=="BDC1");

th_BDC2 = picked_th(picked_nm=="BDC2");
x_BDC2  = picked_x(picked_nm=="BDC2");

% Compute absolute Y-coordinates (ring center at origin, Y-axis upward)
% Since x_piston = Ys_max - Y_s, we have Y_s = Ys_max - x_piston
y_TDC1 = Ys_max - x_TDC1;
y_TDC2 = Ys_max - x_TDC2;
y_BDC1 = Ys_max - x_BDC1;
y_BDC2 = Ys_max - x_BDC2;

% Store all computed values in params structure
params.Ys_max = Ys_max;
params.theta_TDC1 = theta_TDC1;

params.TDC1_th = th_TDC1;
params.TDC1_x  = x_TDC1;
params.TDC1_y  = y_TDC1(1);
params.TDC2_th = th_TDC2;
params.TDC2_x  = x_TDC2;
params.TDC2_y  = y_TDC2(1);
params.BDC1_th = th_BDC1;
params.BDC1_x  = x_BDC1;
params.BDC1_y  = y_BDC1(1);
params.BDC2_th = th_BDC2;
params.BDC2_x  = x_BDC2;
params.BDC2_y  = y_BDC2(1);

% Store all seven dead centers
params.DC_th    = picked_th;
params.DC_x     = picked_x;
params.DC_names = picked_nm;

% Print summary if figures are requested
if generate_figs
    fprintf('\nDead centers from d(x)/dθ = 0:\n');
    fprintf('  TDC1 (highest top):     x = %.6f mm  at ~%.2f° and ~%.2f°\n', ...
        min(x_TDC1), th_TDC1(1), th_TDC1(2));
    fprintf('  TDC2 (lowest top):      x = %.6f mm  at  %.2f° and  %.2f°\n', ...
        mean(x_TDC2), th_TDC2(1), th_TDC2(2));
    fprintf('  BDC1 (highest bottom):  x = %.6f mm  at ~%.2f° and ~%.2f°\n', ...
        max(x_BDC1), th_BDC1(1), th_BDC1(2));
    fprintf('  BDC2 (lowest bottom):   x = %.6f mm  at  %.2f°\n', ...
        x_BDC2, th_BDC2);
    
    % Also list all seven in order
    fprintf('\nSeven DCs (0..1080°):\n');
    for j = 1:numel(targets)
        fprintf('  %-4s at %7.2f° : x = %.6f mm\n', picked_nm(j), picked_th(j), picked_x(j));
    end
    
    % Mark them on the first figure (x vs crank angle)
    figure(findobj('Name','x / A / V vs crank angle'));
    subplot(3,1,1); hold on;
    plot(picked_th*pi/180, picked_x, 'ko', 'MarkerFaceColor','y', 'MarkerSize',6);
    for j = 1:numel(targets)
        text(picked_th(j)*pi/180, picked_x(j), sprintf('  %s', picked_nm(j)), ...
            'VerticalAlignment','bottom','FontSize',8);
    end
end

end

