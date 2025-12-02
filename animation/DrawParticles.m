function DrawParticles(Xp, Yp, Cp, Rp, T, Tmin, Tmax, Xlim, varargin)
% DrawParticles
% Draws circular particles at positions (Xp, Yp) with radius Rp. Particle colors
% depend on temperature T relative to Tmin/Tmax, and can be darkened using Cp.
%
% Usage:
%   DrawParticles(Xp, Yp, Cp, Rp, T, Tmin, Tmax)
%   DrawParticles(Xp, Yp, Cp, Rp, T, Tmin, Tmax, Xlim)
%   DrawParticles(..., 'LineWidth', 0.5, ...)
%
% Inputs:
%   Xp, Yp       - Particle positions
%   Cp           - Darkening factor vector [0,1] for each particle
%                  0 = no darkening, 1 = completely black
%   Rp           - Particle radius
%   T            - Current temperature
%   Tmin, Tmax   - Temperature range for color mapping
%   Xlim         - (optional) if provided, only plot particles within [-Xlim, Xlim]
%
% Color convention (temperature-based, matching ColorCavity):
%   T = Tmin : cold -> blue
%   T = Tmax : hot  -> red
%
% Optional name-value:
%   'EdgeColor'   - [r g b] or 'none', default 'none'
%   'LineWidth'   - scalar, default 0.5
%   'ParticleScale' - scalar, scale factor on particle color intensity (default 1.0)
% -------------------------------------------------------------------------
% Bart Blockmans, 2024 - bart@blockmans.net
% -------------------------------------------------------------------------
    p = inputParser;
    addParameter(p, 'EdgeColor', 'none');
    addParameter(p, 'LineWidth', 0.5);
    addParameter(p, 'ParticleScale', 1.0);
    parse(p, varargin{:});
    edgeColor     = p.Results.EdgeColor;
    lineWidth    = p.Results.LineWidth;
    particleScale = p.Results.ParticleScale;

    if ~isfinite(Rp) || Rp <= 0
        return;
    end

    % Check if Xlim is provided and valid
    useXlim = (nargin >= 8) && ~isempty(Xlim) && isfinite(Xlim) && Xlim > 0;

    % Ensure finite data
    valid = isfinite(Xp) & isfinite(Yp);
    
    % Ensure Cp has the same length as Xp/Yp before filtering
    if numel(Cp) ~= numel(Xp)
        error('Cp must have the same length as Xp and Yp');
    end
    
    Xp = Xp(valid);
    Yp = Yp(valid);
    Cp = Cp(valid);
    Cp = max(0, min(1, Cp)); % Clamp Cp to [0, 1]
    
    % Filter particles by Xlim if specified
    if useXlim
        valid_x = (Xp >= -Xlim) & (Xp <= Xlim);
        Xp = Xp(valid_x);
        Yp = Yp(valid_x);
        Cp = Cp(valid_x);
    end
    
    np = numel(Xp);

    % Compute dimensionless temperature parameter T_norm in [0, 1]
    % Same scheme as ColorCavity: T_norm = 0 -> cold (blue), T_norm = 1 -> hot (red)
    if ~isfinite(T) || ~isfinite(Tmin) || ~isfinite(Tmax) || Tmax <= Tmin
        T_norm = 0.0;
    else
        T_norm = (T - Tmin) / (Tmax - Tmin);
        T_norm = max(0, min(1, T_norm)); % Clamp to [0, 1]
    end

    % Use the same colormap as ColorCavity: blue (cold) -> red (hot)
    % R: 0 -> 1, G: 0, B: 1 -> 0
    baseColor = [T_norm, 0, 1 - T_norm];
    baseColor = max(0, min(1, baseColor)); % Clamp to valid RGB
    
    % Apply particle scale factor
    baseColor = baseColor * particleScale;
    baseColor = max(0, min(1, baseColor));

    % Draw each particle with individual Cp darkening
    for k = 1:np
        % Apply Cp darkening: blend baseColor with black
        % Cp = 0 -> full color, Cp = 1 -> black [0 0 0]
        particleColor = baseColor * (1 - Cp(k));
        particleColor = max(0, min(1, particleColor));
        
        rectangle('Position', [Xp(k) - Rp, Yp(k) - Rp, 2 * Rp, 2 * Rp], ...
                  'Curvature', [1 1], ...
                  'FaceColor', particleColor, ...
                  'EdgeColor', edgeColor, ...
                  'LineWidth', lineWidth);
    end
end



