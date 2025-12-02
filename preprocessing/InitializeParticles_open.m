function [Xp, Yp, Up, Vp] = InitializeParticles_open(Xc, Yc, np, Rp, dt, options)
% InitializeParticles_open
% Randomly place np circular particles inside a closed cavity polygon (Xc, Yc),
% assign initial velocities directed to the right (positive x-direction) with
% small random variations in magnitude. Used for initializing particles in the
% intake manifold that are moving towards the cylinder cavity.
%
% Inputs:
%   Xc, Yc   - cavity boundary polygon (1D arrays, closed or open; will be treated as closed)
%   np       - number of particles
%   Rp       - (optional) particle radius; if omitted, estimated from area and np
%   dt       - (optional) time step; used to scale initial speed RMS
%   options  - (optional struct) fields:
%                .speedRms  - override RMS of initial velocity distribution
%                .rngSeed   - seed for RNG
%
% Outputs:
%   Xp, Yp - particle positions (np x 1)
%   Up, Vp - particle velocities (np x 1), all Up > 0 (directed to the right)
%            Only computed and returned if dt is provided and not empty;
%            otherwise returns empty arrays []

    if nargin < 3
        error('InitializeParticles_open requires at least Xc, Yc, np.');
    end
    if nargin < 4 || isempty(Rp)
        A = polyarea(Xc, Yc);
        phi = 0.04; % target packing fraction
        Rp = sqrt(max(phi * A / (np * pi), eps));
    end
    
    % Check if dt was provided and is not empty
    dtProvided = nargin >= 5 && ~isempty(dt);
    
    speedRms = [];
    if nargin >= 6 && ~isempty(options) && isstruct(options)
        if isfield(options, 'rngSeed') && ~isempty(options.rngSeed)
            rng(options.rngSeed);
        end
        if isfield(options, 'speedRms') && ~isempty(options.speedRms)
            speedRms = options.speedRms;
        end
    end
    
    % Only compute speedRms if dt is provided (needed for velocity computation)
    if dtProvided && isempty(speedRms)
        A = polyarea(Xc, Yc);
        Lchar = sqrt(A / pi);
        speedRms = 0.05 * Lchar / max(dt, eps);
    end

    % Bounding box for rejection sampling
    xmin = min(Xc); xmax = max(Xc);
    ymin = min(Yc); ymax = max(Yc);

    Xp = zeros(np, 1);
    Yp = zeros(np, 1);
    k = 0;
    maxAttempts = max(2000, 50 * np);
    attempts = 0;

    while k < np && attempts < maxAttempts
        attempts = attempts + 1;
        x = xmin + (xmax - xmin) * rand;
        y = ymin + (ymax - ymin) * rand;

        if ~inpolygon(x, y, Xc, Yc)
            continue;
        end

        % Keep at least Rp from the wall
        [dist, ~] = closestDistanceToPolygon(x, y, Xc, Yc);
        if dist < Rp
            continue;
        end

        % Avoid initial overlaps if possible (relaxed if hard to place)
        if k > 0
            d2 = (Xp(1:k) - x).^2 + (Yp(1:k) - y).^2;
            if any(d2 < (2 * Rp)^2)
                continue;
            end
        end

        k = k + 1;
        Xp(k) = x;
        Yp(k) = y;
    end

    % Relax overlap constraint if not enough points placed
    while k < np
        x = xmin + (xmax - xmin) * rand;
        y = ymin + (ymax - ymin) * rand;
        if inpolygon(x, y, Xc, Yc)
            [dist, ~] = closestDistanceToPolygon(x, y, Xc, Yc);
            if dist >= Rp
                k = k + 1;
                Xp(k) = x;
                Yp(k) = y;
            end
        end
    end

    % Assign initial velocities only if dt was provided
    if dtProvided
        % Assign initial velocities: all directed to the right with small random variations
        % X-velocity: positive with ~10% random variation around base speed
        Up = speedRms * abs(1 + 0.1 * randn(np, 1));
        % Y-velocity: small random component (5% of base speed) for slight perpendicular motion
        Vp = 0.05 * speedRms * randn(np, 1);
    else
        % Return empty arrays if dt was not provided
        Up = [];
        Vp = [];
    end

end

function [dist, segIdx] = closestDistanceToPolygon(x, y, Xc, Yc)
% Returns the minimum distance from point (x,y) to polygon (Xc,Yc),
% and the index of the closest segment.
    n = numel(Xc);
    if n < 2
        dist = inf;
        segIdx = 0;
        return;
    end
    xc = Xc(:);
    yc = Yc(:);
    % Ensure closed
    if xc(1) ~= xc(end) || yc(1) ~= yc(end)
        xc = [xc; xc(1)];
        yc = [yc; yc(1)];
    end
    bestD = inf;
    bestIdx = 1;
    p = [x; y];
    for i = 1:numel(xc)-1
        a = [xc(i); yc(i)];
        b = [xc(i+1); yc(i+1)];
        [cp, ~] = closestPointOnSegment(a, b, p);
        d = hypot(p(1) - cp(1), p(2) - cp(2));
        if d < bestD
            bestD = d;
            bestIdx = i;
        end
    end
    dist = bestD;
    segIdx = bestIdx;
end

function [cp, t] = closestPointOnSegment(a, b, p)
% Closest point cp on segment ab to point p; also returns param t in [0,1]
    ab = b - a;
    ab2 = dot(ab, ab);
    if ab2 <= eps
        cp = a;
        t = 0;
        return;
    end
    t = dot(p - a, ab) / ab2;
    t = max(0, min(1, t));
    cp = a + t * ab;
end

