function DrawParticlesAndCavity(Xp, Yp, Rp, V, Vmin, Vmax, Xc, Yc, Xlim, varargin)
% DrawParticles
% Draws circular particles at positions (Xp, Yp) with radius Rp and fills the
% cavity (Xc, Yc) with a translucent color that depends on the instantaneous
% cylinder volume V relative to its min/max over the cycle.
%
% Usage:
%   DrawParticles(Xp, Yp, Rp, V, Vmin, Vmax, Xc, Yc)
%   DrawParticles(Xp, Yp, Rp, V, Vmin, Vmax, Xc, Yc, Xlim)
%   DrawParticles(..., 'LineWidth', 0.5, ...)
%
% Color convention (heuristic, volume-based):
%   V close to Vmin : strong compression  -> hotter -> more towards red
%   V close to Vmax : strong expansion    -> cooler -> more towards blue
%
% Optional inputs:
%   Xlim          - if provided, only plot particles and cylinder within [-Xlim, Xlim]
%
% Optional name-value:
%   'EdgeColor'   - [r g b] or 'none', default [0 0 0]
%   'LineWidth'   - scalar, default 0.5
%   'AlphaFill'   - scalar in [0,1], fill transparency (default 0.20)
%   'ParticleScale' - scalar, scale factor on particle color intensity (default 1.0)
% -------------------------------------------------------------------------
% Bart Blockmans, 2024 - bart@blockmans.net
% -------------------------------------------------------------------------

    p = inputParser;
    addParameter(p, 'EdgeColor', 'none');
    addParameter(p, 'LineWidth', 0.5);
    addParameter(p, 'AlphaFill', 0.20);
    addParameter(p, 'ParticleScale', 1.0);
    parse(p, varargin{:});
    edgeColor     = p.Results.EdgeColor;
    lineWidth    = p.Results.LineWidth;
    alphaFill    = p.Results.AlphaFill;
    particleScale = p.Results.ParticleScale;

    if ~isfinite(Rp) || Rp <= 0
        return;
    end

    % Check if Xlim is provided and valid
    useXlim = (nargin >= 9) && ~isempty(Xlim) && isfinite(Xlim) && Xlim > 0;

    % Ensure finite data
    valid = isfinite(Xp) & isfinite(Yp);
    Xp = Xp(valid);
    Yp = Yp(valid);
    
    % Filter particles by Xlim if specified
    if useXlim
        valid_x = (Xp >= -Xlim) & (Xp <= Xlim);
        Xp = Xp(valid_x);
        Yp = Yp(valid_x);
    end
    
    np = numel(Xp);

    % Derive a frame-wise color from volume V using global [Vmin, Vmax]
    if ~isfinite(V) || ~isfinite(Vmin) || ~isfinite(Vmax) || Vmax <= Vmin
        xi = 0.0;
    else
        % Compression ratio proxy: 0 at max volume (cool), 1 at min volume (hot)
        xi = (Vmax - V) / (Vmax - Vmin);
        xi = max(0, min(1, xi));
    end

    % Map xi in [0,1] to a blue–red ramp (blue=expanded, red=compressed)
    % t=0 -> blue (V ~ Vmax), t=1 -> red (V ~ Vmin)
    t = xi;
    cold = [0.0 0.25 0.8];            % cool blue
    hot  = [0.9 0.15 0.0];            % hot red
    baseColor = (1 - t) * cold + t * hot;
    baseColor = max(0, min(1, baseColor)); % clamp to valid RGB

    % Cavity fill: translucent, slightly desaturated
    fillColor = 0.9 * baseColor + 0.1 * [1 1 1]; % a bit lighter
    fillColor = max(0, min(1, fillColor));

    if nargin >= 7 && ~isempty(Xc) && ~isempty(Yc)
        % Ensure column vectors
        Xcv = Xc(:);
        Ycv = Yc(:);
        
        % Clip cylinder polygon to Xlim if specified
        if useXlim
            [Xcv, Ycv] = clipPolygonToXRange(Xcv, Ycv, -Xlim, Xlim);
        end
        
        % Draw cavity fill underneath particles (only if polygon is not empty)
        if ~isempty(Xcv) && ~isempty(Ycv) && numel(Xcv) >= 3
            fill(Xcv, Ycv, fillColor, ...
                 'EdgeColor', 'none', ...
                 'FaceAlpha', alphaFill);
            hold on;
        end
    end

    % Particle color: slightly more saturated/intense than fill
    faceColor = baseColor * particleScale;
    faceColor = max(0, min(1, faceColor));

    for k = 1:np
        rectangle('Position', [Xp(k) - Rp, Yp(k) - Rp, 2 * Rp, 2 * Rp], ...
                  'Curvature', [1 1], ...
                  'FaceColor', faceColor, ...
                  'EdgeColor', edgeColor, ...
                  'LineWidth', lineWidth);
    end
end

function [Xout, Yout] = clipPolygonToXRange(X, Y, xmin, xmax)
% Clips a polygon to the x-range [xmin, xmax]
% Returns the clipped polygon coordinates
    if isempty(X) || isempty(Y) || numel(X) < 3
        Xout = [];
        Yout = [];
        return;
    end
    
    % Ensure polygon is closed
    Xv = X(:);
    Yv = Y(:);
    if Xv(1) ~= Xv(end) || Yv(1) ~= Yv(end)
        Xv = [Xv; Xv(1)];
        Yv = [Yv; Yv(1)];
    end
    
    n = numel(Xv);
    Xout = [];
    Yout = [];
    
    for i = 1:n-1
        x1 = Xv(i);
        y1 = Yv(i);
        x2 = Xv(i+1);
        y2 = Yv(i+1);
        
        % Check if both points are inside
        p1_inside = (x1 >= xmin) && (x1 <= xmax);
        p2_inside = (x2 >= xmin) && (x2 <= xmax);
        
        if p1_inside && p2_inside
            % Both points inside, add the segment
            if isempty(Xout) || Xout(end) ~= x1 || Yout(end) ~= y1
                Xout = [Xout; x1];
                Yout = [Yout; y1];
            end
            Xout = [Xout; x2];
            Yout = [Yout; y2];
        elseif p1_inside && ~p2_inside
            % p1 inside, p2 outside - add p1 and intersection with boundary
            if isempty(Xout) || Xout(end) ~= x1 || Yout(end) ~= y1
                Xout = [Xout; x1];
                Yout = [Yout; y1];
            end
            % Find intersection with xmin or xmax
            if x2 < xmin
                % Intersect with xmin
                if abs(x2 - x1) > eps
                    t = (xmin - x1) / (x2 - x1);
                    y_int = y1 + t * (y2 - y1);
                    Xout = [Xout; xmin];
                    Yout = [Yout; y_int];
                end
            elseif x2 > xmax
                % Intersect with xmax
                if abs(x2 - x1) > eps
                    t = (xmax - x1) / (x2 - x1);
                    y_int = y1 + t * (y2 - y1);
                    Xout = [Xout; xmax];
                    Yout = [Yout; y_int];
                end
            end
        elseif ~p1_inside && p2_inside
            % p1 outside, p2 inside - add intersection and p2
            % Find intersection with xmin or xmax
            if x1 < xmin
                % Intersect with xmin
                if abs(x2 - x1) > eps
                    t = (xmin - x1) / (x2 - x1);
                    y_int = y1 + t * (y2 - y1);
                    Xout = [Xout; xmin];
                    Yout = [Yout; y_int];
                end
            elseif x1 > xmax
                % Intersect with xmax
                if abs(x2 - x1) > eps
                    t = (xmax - x1) / (x2 - x1);
                    y_int = y1 + t * (y2 - y1);
                    Xout = [Xout; xmax];
                    Yout = [Yout; y_int];
                end
            end
            Xout = [Xout; x2];
            Yout = [Yout; y2];
        elseif ~p1_inside && ~p2_inside
            % Both outside - check if segment crosses the range
            x_seg_min = min(x1, x2);
            x_seg_max = max(x1, x2);
            if x_seg_min < xmax && x_seg_max > xmin
                % Segment crosses the range - find intersections
                intersections = [];
                t_vals = [];
                % Check intersection with xmin
                if x_seg_min < xmin && x_seg_max >= xmin
                    if abs(x2 - x1) > eps
                        t = (xmin - x1) / (x2 - x1);
                        if t >= 0 && t <= 1
                            y_int = y1 + t * (y2 - y1);
                            intersections = [intersections; xmin, y_int];
                            t_vals = [t_vals; t];
                        end
                    end
                end
                % Check intersection with xmax
                if x_seg_min <= xmax && x_seg_max > xmax
                    if abs(x2 - x1) > eps
                        t = (xmax - x1) / (x2 - x1);
                        if t >= 0 && t <= 1
                            y_int = y1 + t * (y2 - y1);
                            intersections = [intersections; xmax, y_int];
                            t_vals = [t_vals; t];
                        end
                    end
                end
                % Add intersections in order along the segment (by parameter t)
                if size(intersections, 1) == 2
                    [~, idx] = sort(t_vals);
                    Xout = [Xout; intersections(idx(1), 1); intersections(idx(2), 1)];
                    Yout = [Yout; intersections(idx(1), 2); intersections(idx(2), 2)];
                elseif size(intersections, 1) == 1
                    Xout = [Xout; intersections(1, 1)];
                    Yout = [Yout; intersections(1, 2)];
                end
            end
        end
    end
    
    % Remove duplicate consecutive points
    if ~isempty(Xout)
        keep = true(size(Xout));
        for i = 2:numel(Xout)
            if abs(Xout(i) - Xout(i-1)) < eps && abs(Yout(i) - Yout(i-1)) < eps
                keep(i) = false;
            end
        end
        Xout = Xout(keep);
        Yout = Yout(keep);
    end
    
    % Ensure we have at least 3 points for a valid polygon
    if numel(Xout) < 3
        Xout = [];
        Yout = [];
    end
end


