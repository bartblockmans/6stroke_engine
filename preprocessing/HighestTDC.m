function [Ys_max, theta_TDC1_deg] = HighestTDC(Rr, Rp, delta, L, psi0)
%HIGHESTTDC  Compute highest TDC of the slider: Ys_max and its angle.
% Finds the 3rd zero of dYs/dtheta (theta in deg) over [0,1080] and
% returns Ys at that angle.
%
% Inputs (mm, rad):
%   Rr, Rp  : ring/planet pitch radii [mm]
%   delta   : offset of tracked planet point [mm]
%   L       : rod length to Y-axis slider [mm]
%   psi0    : body-fixed angle of tracked point [rad]
%
% Outputs:
%   Ys_max            : maximum Y of slider over a full 6-stroke [mm]
%   theta_TDC1_deg    : crank angle (deg) at which Ys_max occurs
% -------------------------------------------------------------------------
% Bart Blockmans, 2024 - bart@blockmans.net
% -------------------------------------------------------------------------


if nargin < 5 || isempty(psi0), psi0 = pi/2; end

% Mechanism constants
a = Rr - Rp;          % [mm]
k = a / Rp;           % (Rr - Rp)/Rp

% Analytic Ys(theta_deg) and dYs/dtheta (per degree) function handles
Ys_fun  = @(td) Ys_at_deg(td, a, k, delta, L, psi0);
dYs_fun = @(td) dYs_dtheta_deg(td, a, k, delta, L, psi0);

% Coarse scan to bracket zeros of dYs/dtheta in [0,1080] deg
th = linspace(0,1080,4097);   % dense & cheap
dY = arrayfun(dYs_fun, th);
s  = sign(dY);
ok = ~isnan(s);
signchg = ok(1:end-1) & ok(2:end) & (s(1:end-1).*s(2:end) < 0);
brk = find(signchg);

% Refine each bracket with fzero
roots_deg = zeros(1,numel(brk));
for i = 1:numel(brk)
    a_deg = th(brk(i));
    b_deg = th(brk(i)+1);
    try
        roots_deg(i) = fzero(dYs_fun, [a_deg, b_deg]);
    catch
        roots_deg(i) = 0.5*(a_deg + b_deg); % fallback
    end
end

% Add endpoints if derivative is numerically zero there
if abs(dYs_fun(0))     < 1e-12, roots_deg(end+1) = 0;      end %#ok<AGROW>
if abs(dYs_fun(1080))  < 1e-12, roots_deg(end+1) = 1080;   end %#ok<AGROW>

% Sort & unique
roots_deg = sort(roots_deg);
if ~isempty(roots_deg)
    keep = [true, diff(roots_deg) > 1e-6];
    roots_deg = roots_deg(keep);
end

% We expect 7 roots; the 3rd one is the highest TDC
if numel(roots_deg) < 3
    % Robust fallback: take maximum of Ys over dense grid
    [Ys_max, imax] = max(arrayfun(Ys_fun, th));
    theta_TDC1_deg = th(imax);
    return
end

theta_TDC1_deg = roots_deg(3);
Ys_max         = Ys_fun(theta_TDC1_deg);

end

% --- helpers (analytic, degree inputs) ---

function Ys = Ys_at_deg(theta_deg, a, k, delta, L, psi0)
    phi = deg2rad(theta_deg) - 3*pi/2;  % carrier shift
    h   = -k.*phi + psi0;
    Xr  = a.*cos(phi) + delta.*cos(h);
    Yr  = a.*sin(phi) + delta.*sin(h);
    arg = L.^2 - Xr.^2;
    if arg < 0, Ys = NaN; return; end
    Ys = Yr + sqrt(arg);
end

function dYs = dYs_dtheta_deg(theta_deg, a, k, delta, L, psi0)
    phi = deg2rad(theta_deg) - 3*pi/2;  % rad
    h   = -k.*phi + psi0;
    % Xr, Yr and their d/dphi
    Xr  = a.*cos(phi) + delta.*cos(h);
    Yr  = a.*sin(phi) + delta.*sin(h);
    dXr_dphi = -a.*sin(phi) + delta.*k.*sin(h);
    dYr_dphi =  a.*cos(phi) - delta.*k.*cos(h);
    den = L.^2 - Xr.^2;
    if den <= 0
        dYs = NaN; return
    end
    dYs_dphi = dYr_dphi - (Xr .* dXr_dphi) ./ sqrt(den);
    dYs = dYs_dphi * (pi/180);          % per degree
end
