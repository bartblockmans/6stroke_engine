function [ind_in, ind_out] = SplitParticles(Xp, Yp, Xc, Yc, Rp)
% SplitParticles
% Determines which particles are inside the cylinder cavity and which are outside.
% Uses a safety margin to account for finite precision in collision detection,
% so particles very close to the boundary are still considered inside.
%
% Inputs:
%   Xp, Yp - particle coordinates (np x 1 arrays)
%   Xc, Yc - closed cylinder cavity polygon coordinates (1D arrays)
%   Rp     - (optional) particle radius; if omitted, estimated from geometry
%
% Outputs:
%   ind_in  - indices of particles inside the cylinder cavity
%   ind_out - indices of particles outside the cylinder cavity

    % Compute characteristic length and tolerance (similar to UpdateParticles.m)
    A = max(polyarea(Xc, Yc), eps);
    Lchar = sqrt(A / pi);
    tol = 1e-6 * Lchar;
    
    % Get particle radius if not provided
    if nargin < 5 || isempty(Rp)
        % Estimate from typical packing fraction
        np = numel(Xp);
        phi = 0.04;
        Rp = sqrt(max(phi * A / (np * pi), eps));
    end
    
    % Safety margin: account for particle radius + tolerance
    % In UpdateParticles.m, particles are positioned at RpEff + tol from wall,
    % so particle centers can be slightly outside but particles are still inside
    safety_margin = Rp + 2 * tol;
    
    % Ensure polygon is closed for distance calculation
    Xcv = Xc(:);
    Ycv = Yc(:);
    if Xcv(1) ~= Xcv(end) || Ycv(1) ~= Ycv(end)
        Xcv = [Xcv; Xcv(1)];
        Ycv = [Ycv; Ycv(1)];
    end
    
    % Check which particles are inside or on the boundary of the cylinder cavity
    [in, on] = inpolygon(Xp, Yp, Xc, Yc);
    
    % Consider particles on the boundary as inside
    inside = in | on;
    
    % For particles marked as outside, check if they're within safety margin
    % This accounts for finite precision in collision detection
    outside_idx = find(~inside);
    if ~isempty(outside_idx)
        % Check distance to polygon for particles marked as outside
        for i = 1:length(outside_idx)
            idx = outside_idx(i);
            dist = closestDistanceToPolygon(Xp(idx), Yp(idx), Xcv, Ycv);
            % If particle center is within safety margin of boundary, consider it inside
            % This accounts for the particle radius and numerical tolerance
            if dist < safety_margin
                inside(idx) = true;
            end
        end
    end
    
    % Get indices of particles inside and outside
    ind_in = find(inside);
    ind_out = find(~inside);

end

function dist = closestDistanceToPolygon(x, y, Xc, Yc)
% Returns the minimum distance from point (x,y) to polygon (Xc,Yc)
    n = numel(Xc);
    if n < 2
        dist = inf;
        return;
    end
    bestD = inf;
    p = [x; y];
    for i = 1:n-1
        a = [Xc(i); Yc(i)];
        b = [Xc(i+1); Yc(i+1)];
        [cp, ~] = closestPointOnSegment(a, b, p);
        d = hypot(p(1) - cp(1), p(2) - cp(2));
        if d < bestD
            bestD = d;
        end
    end
    dist = bestD;
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

