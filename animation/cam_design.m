function params = cam_design(params, results)

% Intake cam
% -------------------------------------------------------------------------

% Get displacements & normalize
IV_disp = (results.L_int_m / max(results.L_int_m)) * params.valve_disp;

% 2pi angle
angle = linspace(0, 2*pi, length(IV_disp));

% Columnize
IV_disp = IV_disp(:);
angle = angle(:);

% Downsample
theta = linspace(0,2*pi,100);
IV_disp = interp1(angle, IV_disp, theta);

% Compute cam radius at each angle
R_IV = (0.5 * params.D_valve) + IV_disp;

% Compute cam profile
X_IV = R_IV .* cos(theta);
Y_IV = R_IV .* sin(theta);

% Smooth cam profile
[X_IV, Y_IV] = smooth_cam_profile(X_IV, Y_IV, 'Window', 21);

% Store
params.IV_cam_X = X_IV';
params.IV_cam_Y = Y_IV';

% Exhaust cam
% -------------------------------------------------------------------------

% Get displacements & normalize
EV_disp = (results.L_exh_m / max(results.L_exh_m)) * params.valve_disp;

% 2pi angle
angle = linspace(0, 2*pi, length(EV_disp));

% Columnize
EV_disp = EV_disp(:);
angle = angle(:);

% Downsample
theta = linspace(0,2*pi,100);
EV_disp = interp1(angle, EV_disp, theta);

% Compute cam radius at each angle
R_EV = (0.5 * params.D_valve) + EV_disp;

% Compute cam profile
X_EV = R_EV .* cos(theta);
Y_EV = R_EV .* sin(theta);

% Smooth cam profile
[X_EV, Y_EV] = smooth_cam_profile(X_EV, Y_EV, 'Window', 21);

% Store
params.EV_cam_X = X_EV';
params.EV_cam_Y = Y_EV';

function [X_s, Y_s] = smooth_cam_profile(X, Y, varargin)
%SMOOTH_CAM_PROFILE  Smooth/simplify a closed cam profile in XY.
%
%   [X_s, Y_s] = smooth_cam_profile(X, Y)
%   [X_s, Y_s] = smooth_cam_profile(X, Y, 'Window', 31, 'PreserveExtrema', true)
%
%   Inputs:
%     X, Y   : cam profile coordinates (closed curve; first ~ last)
%
%   Optional name-value:
%     'Window'         : smoothing window length (odd, default: 21)
%     'PreserveExtrema': true/false, preserve min & max radius (default: true)
%
%   Outputs:
%     X_s, Y_s : smoothed cam profile coordinates (same length)

    % --- Parse inputs ----------------------------------------------------
    p = inputParser;
    addParameter(p, 'Window', 21, ...
        @(w) isnumeric(w) && isscalar(w) && w >= 3);
    addParameter(p, 'PreserveExtrema', true, @(x)islogical(x) && isscalar(x));
    parse(p, varargin{:});
    win            = p.Results.Window;
    preserveExtrem = p.Results.PreserveExtrema;

    X = X(:);
    Y = Y(:);

    if numel(X) ~= numel(Y)
        error('smooth_cam_profile:SizeMismatch', ...
              'X and Y must have the same length.');
    end

    N = numel(X);

    % Ensure window is odd and not larger than N
    win = min(win, N);
    if mod(win,2) == 0
        win = win - 1;
    end
    if win < 3
        win = 3;
    end

    % --- Periodic padding to avoid edge artefacts ------------------------
    X_pad = [X; X; X];
    Y_pad = [Y; Y; Y];

    % --- Smooth using Savitzky-Golay filter ------------------------------
    X_pad_s = smoothdata(X_pad, 'sgolay', win);
    Y_pad_s = smoothdata(Y_pad, 'sgolay', win);

    % Back to original length
    X_s = X_pad_s(N+1:2*N);
    Y_s = Y_pad_s(N+1:2*N);

    % --- Optional: preserve min & max radius -----------------------------
    if preserveExtrem
        R_raw = hypot(X,   Y);
        R_s   = hypot(X_s, Y_s);

        Rmin_raw = min(R_raw);
        Rmax_raw = max(R_raw);

        Rmin_s = min(R_s);
        Rmax_s = max(R_s);

        % Avoid division by zero in pathological cases
        if abs(Rmax_s - Rmin_s) > eps
            % Affine mapping: same min & max as original
            scale = (Rmax_raw - Rmin_raw) / (Rmax_s - Rmin_s);
            R_s_aff = Rmin_raw + (R_s - Rmin_s) * scale;

            % Rescale X_s, Y_s radially
            idx_nonzero = R_s > eps;
            X_s(idx_nonzero) = X_s(idx_nonzero) ./ R_s(idx_nonzero) .* R_s_aff(idx_nonzero);
            Y_s(idx_nonzero) = Y_s(idx_nonzero) ./ R_s(idx_nonzero) .* R_s_aff(idx_nonzero);
        end
    end

    % --- Re-enforce closure ----------------------------------------------
    X_s(1) = X_s(end);
    Y_s(1) = Y_s(end);
end

end

