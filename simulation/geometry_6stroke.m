function [V, A, x_piston, dVdth, dAdth, dxdth, coords] = geometry_6stroke(theta, params)
%GEOMETRY_6STROKE  6-stroke cylinder geometry from internal-gear + rod slider.
%
% Inputs (besides theta these are all defined inside params)
%   theta : crank/carrier angle [deg], scalar or vector (typically 0..1080)
%   Rr    : ring pitch radius [mm]
%   Rp    : planet pitch radius [mm]
%   delta : offset of tracked point on planet [mm]
%   L     : rod length from the red point to Y-axis slider [mm]
%   B     : bore diameter [mm]
%   Vc    : clearance volume [mm^3]
%   psi0  : (optional) body-fixed angle of tracked point [rad]  (default pi/2)
%
% Assumptions: theta0=0, phi0=0, branch s=+1, phi = theta_rad - 3*pi/2
%
% Outputs (same length as theta)
%   V       : volume [mm^3]
%   A       : wetted area [mm^2]
%   x_piston: piston displacement from highest TDC [mm] (>=0)
%   dVdth   : dV/dtheta [mm^3/deg]
%   dAdth   : dA/dtheta [mm^2/deg]
%   dxdth   : dx/dtheta [mm/deg]
%   coords  : struct with absolute coordinates (ring center at origin)
%             .X_p, .Y_p : planet center
%             .X_r, .Y_r : red point
%             .X_s, .Y_s : slider point (X_s == 0)
% -------------------------------------------------------------------------
% Bart Blockmans, 2024 - bart@blockmans.net
% -------------------------------------------------------------------------


% Get parameters
Rr = params.Rr;
Rp = params.Rp;
delta = params.delta;
L = params.L;
B = params.B;
Vc = params.Vc;
psi0 = params.psi0;
if isfield(params, "Ys_max")
    Y_s_max = params.Ys_max;
else
    Y_s_max = [];
end

% Geometry & kinematics parameters
a = Rr - Rp;                 % [mm]
k = a / Rp;                  % = (Rr - Rp)/Rp

% Angles
theta = theta(:).';                          % row
phi   = deg2rad(theta) - 3*pi/2;             % align -270..+810 deg

% Red point (planet-attached) and planet center
h   = -k.*phi + psi0;                        % inner phase
X_p = a.*cos(phi);
Y_p = a.*sin(phi);
X_r = X_p + delta.*cos(h);
Y_r = Y_p + delta.*sin(h);

% Slider (rod to Y-axis, choose upper branch s=+1)
arg = L.^2 - X_r.^2;
if any(arg < 0)
    warning('geometry_6stroke:RodTooShort', ...
        'For some theta, L^2 - x^2 < 0; returning NaN for those angles.');
end
den_sqrt = sqrt(max(arg, 0));
Y_s = Y_r + den_sqrt;                         % branch above red point
X_s = zeros(size(Y_s));                       % on Y-axis

% Highest TDC (allow standalone calls without extrema pre-pass)
if isempty(Y_s_max)
    [Y_s_max, ~] = HighestTDC(Rr, Rp, delta, L, psi0);
    % Y_s_max  = max(Y_s);
end

% Piston displacement from highest TDC
x_piston = Y_s_max - Y_s;                     % >= 0

% Volume & area
V = Vc + pi*(B^2/4) .* x_piston;              % [mm^3]
A = pi*B .* (B/2 + x_piston);                 % [mm^2]

% Derivatives wrt theta (degrees)
% First derivatives wrt phi
dXr_dphi = -a.*sin(phi) + delta.*k.*sin(h);
dYr_dphi =  a.*cos(phi) - delta.*k.*cos(h);

% dy_s/dphi = dy/dphi - (x * dx/dphi)/sqrt(L^2 - x^2)
den = den_sqrt; den_safe = den; den_safe(den_safe < eps) = eps;
dYs_dphi = dYr_dphi - (X_r .* dXr_dphi) ./ den_safe;

% Chain rule: d/dtheta = d/dphi * (pi/180)
dxdth = -(dYs_dphi) * (pi/180);               % [mm/deg]
dVdth =  (pi*(B^2/4)) .* dxdth;               % [mm^3/deg]
dAdth =  (pi*B) .* dxdth;                     % [mm^2/deg]

% Pack coordinates (column vectors)
coords = struct();
coords.X_p = X_p(:);  coords.Y_p = Y_p(:);
coords.X_r = X_r(:);  coords.Y_r = Y_r(:);
coords.X_s = X_s(:);  coords.Y_s = Y_s(:);

% Shape outputs as column vectors
V        = V(:);
A        = A(:);
x_piston = x_piston(:);
dVdth    = dVdth(:);
dAdth    = dAdth(:);
dxdth    = dxdth(:);

end

