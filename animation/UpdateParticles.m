function [Xp, Yp, Up, Vp] = UpdateParticles(Xp, Yp, Up, Vp, Xc, Yc, Rp, dt, pistonEdgeIdx, v_piston, part_contact, opts)
% UpdateParticles (moving-wall & gas-like refinement)
% Advances particle positions and velocities by one time step:
% - Integrate motion with dt (sub-stepped)
% - Reflect off cavity walls; piston face treated as a moving wall (if provided)
% - Handle particle-particle near-elastic collisions with small tangential damping
% - Thermal (kinetic) velocity cap based on RMS speed
%
% Inputs:
%   Xp, Yp - positions (np x 1)
%   Up, Vp - velocities (np x 1)
%   Xc, Yc - cavity polygon (closed or open; treated as closed)
%   Rp     - draw radius (collisions use RpEff if density too high)
%   dt     - time step
%   pistonEdgeIdx - (optional) indices of segments forming piston face (in closed polygon)
%   v_piston      - (optional) piston normal speed [same units as Up,Vp], positive into gas
%   opts          - (optional) struct with fields:
%       .numRelaxIters (default 2)
%       .substepMax    (default 40)
%       .phiCap        (default 0.68)
%       .betaT         (default 0.04)
%       .eRest         (default 0.98)
%       .vCapFactor    (default 2.5)
%
% Outputs:
%   Xp, Yp - updated positions
%   Up, Vp - updated velocities
% -------------------------------------------------------------------------
% Bart Blockmans, 2024 - bart@blockmans.net
% -------------------------------------------------------------------------

    % Defaults and guards
    if nargin < 8 || isempty(dt)
        dt = 0.01;
    end
    if nargin < 7 || isempty(Rp)
        A = polyarea(Xc, Yc);
        np = numel(Xp);
        phi = 0.04;
        Rp = sqrt(max(phi * A / (np * pi), eps));
    end
    if nargin < 9 || isempty(pistonEdgeIdx)
        pistonEdgeIdx = [];
    end
    if nargin < 10 || isempty(v_piston)
        v_piston = 0;
    end
    if nargin < 11 || isempty(part_contact)
        part_contact = 1;
    end
    if nargin < 12 || isempty(opts)
        opts = struct();
    end

    np = numel(Xp);
    if any([numel(Yp) numel(Up) numel(Vp)] ~= np)
        error('Xp, Yp, Up, Vp must have the same length.');
    end

    % Tunables
    numRelaxIters = getOpt(opts, 'numRelaxIters', 2);
    substepMax    = getOpt(opts, 'substepMax', 40);
    phiCap        = getOpt(opts, 'phiCap', 0.68);
    betaT         = getOpt(opts, 'betaT', 0.04);
    eRest         = getOpt(opts, 'eRest', 0.98);
    vCapFactor    = getOpt(opts, 'vCapFactor', 2.5);

    % Ensure polygon is closed
    Xcv = Xc(:); Ycv = Yc(:);
    if Xcv(1) ~= Xcv(end) || Ycv(1) ~= Ycv(end)
        Xcv = [Xcv; Xcv(1)];
        Ycv = [Ycv; Ycv(1)];
    end
    nSeg = numel(Xcv) - 1;
    pistonMask = false(nSeg, 1);
    if ~isempty(pistonEdgeIdx)
        pistonEdgeIdx = unique(max(1, min(nSeg, round(pistonEdgeIdx(:)))));
        pistonMask(pistonEdgeIdx) = true;
    end
    A = max(polyarea(Xc, Yc), eps);
    % Effective collision radius if instantaneous packing is too high
    % This preserves solver feasibility at high compression without changing draw radius
    phiInst = (np * pi * (Rp^2)) / A;
    if phiInst > phiCap
        RpEff = Rp * sqrt(phiCap / phiInst);
    else
        RpEff = Rp;
    end

    % Characteristic length and tolerance (avoid sticking due to tiny penetrations)
    Lchar = sqrt(max(A / pi, eps));
    tol = 1e-6 * Lchar;

    % Thermal (kinetic) cap – based on current RMS speed
    vrms = sqrt(mean(Up.^2 + Vp.^2));
    if ~isfinite(vrms) || vrms == 0
        vrms = 1e-6;
    end
    vCap = vCapFactor * vrms;

    % Sub-stepping to limit displacement per step
    maxSpeed = max(hypot(Up, Vp));
    if ~isfinite(maxSpeed)
        maxSpeed = 0;
    end
    maxDisp = 0.25 * RpEff;
    nSub = max(1, ceil((maxSpeed * dt) / max(maxDisp, eps)));
    nSub = min(nSub, substepMax);
    subDt = dt / nSub;

    for sub = 1:nSub
        % Integrate motion
        Xp = Xp + Up * subDt;
        Yp = Yp + Vp * subDt;

        % Wall collisions: only reflect if moving outward or outside
        inside = inpolygon(Xp, Yp, Xc, Yc);
        for i = 1:np
            [cp, n_in, dist, eIdx] = closestPointAndNormalToPolygon(Xp(i), Yp(i), Xcv, Ycv);
            n_hat = n_in(:).';
            v = [Up(i), Vp(i)];
            % Wall velocity: piston segments move with v_piston along n_in
            if ~isempty(pistonEdgeIdx) && eIdx >= 1 && eIdx <= nSeg && pistonMask(eIdx)
                uw = v_piston * n_hat;
            else
                uw = [0, 0];
            end
            vrel = v - uw;
            vnrel = dot(vrel, n_hat);
            needReflect = (~inside(i)) || (dist < (RpEff - tol) && vnrel < 0);
            if needReflect
                % Reposition slightly inside, then reflect normal component
                Xp(i) = cp(1) + n_hat(1) * (RpEff + tol);
                Yp(i) = cp(2) + n_hat(2) * (RpEff + tol);
                % Reflect in wall frame with near-elastic restitution and tangential damping
                vt_rel = vrel - vnrel * n_hat;
                vrel_new = vt_rel - eRest * vnrel * n_hat;
                v_new = uw + vrel_new;
                % Apply small tangential damping
                t_hat = [-n_hat(2), n_hat(1)];
                vt_comp = dot(v_new, t_hat);
                v_new = v_new - betaT * vt_comp * t_hat;
                Up(i) = v_new(1); Vp(i) = v_new(2);
            end
        end

        % Particle-particle collisions (near-elastic, equal mass)
        if part_contact
            for iter = 1:numRelaxIters
                for i = 1:np-1
                    for j = i+1:np
                        dx = Xp(j) - Xp(i);
                        dy = Yp(j) - Yp(i);
                        d = hypot(dx, dy);
                        if d <= tol
                            % Deterministic separation along line of centers or pseudo-normal
                            if dx ~= 0 || dy ~= 0
                                n = [dx, dy] / max(d, eps);
                            else
                                ang = 2 * pi * frac(0.618 * (i + j));
                                n = [cos(ang), sin(ang)];
                            end
                            overlap = 2 * RpEff;
                            corr = 0.5 * overlap + 0.5 * tol;
                            Xp(i) = Xp(i) - corr * n(1);
                            Yp(i) = Yp(i) - corr * n(2);
                            Xp(j) = Xp(j) + corr * n(1);
                            Yp(j) = Yp(j) + corr * n(2);
                            % Zero approaching normal relative speed
                            rel = [Up(j) - Up(i), Vp(j) - Vp(i)];
                            un  = dot(rel, n);
                            if un < 0
                                Up(i) = Up(i) + 0.5 * un * n(1);
                                Vp(i) = Vp(i) + 0.5 * un * n(2);
                                Up(j) = Up(j) - 0.5 * un * n(1);
                                Vp(j) = Vp(j) - 0.5 * un * n(2);
                            end
                        elseif d < 2 * RpEff - tol
                            n = [dx, dy] / d;
                            overlap = (2 * RpEff - d);
                            % Separate positions
                            corr = 0.5 * overlap + 0.5 * tol;
                            Xp(i) = Xp(i) - corr * n(1);
                            Yp(i) = Yp(i) - corr * n(2);
                            Xp(j) = Xp(j) + corr * n(1);
                            Yp(j) = Yp(j) + corr * n(2);
                            % Equal-mass near-elastic collision if approaching
                            v1 = [Up(i), Vp(i)];
                            v2 = [Up(j), Vp(j)];
                            rel = v2 - v1;
                            un  = dot(rel, n);
                            if un < 0
                                Jn = -(1 + eRest) * un / 2; % equal masses
                                dv = Jn * n;
                                v1p = v1 - dv;
                                v2p = v2 + dv;
                                % Tangential damping
                                t_hat = [-n(2), n(1)];
                                vt_i = dot(v1p, t_hat);
                                vt_j = dot(v2p, t_hat);
                                v1p = v1p - betaT * vt_i * t_hat;
                                v2p = v2p - betaT * vt_j * t_hat;
                                Up(i) = v1p(1); Vp(i) = v1p(2);
                                Up(j) = v2p(1); Vp(j) = v2p(2);
                            end
                        end
                    end
                end

                % Re-enforce walls after pair resolution within the same sub-step
                inside = inpolygon(Xp, Yp, Xc, Yc);
                for i = 1:np
                    [cp, n_in, dist, eIdx] = closestPointAndNormalToPolygon(Xp(i), Yp(i), Xcv, Ycv);
                    n_hat = n_in(:).';
                    v = [Up(i), Vp(i)];
                    if ~isempty(pistonEdgeIdx) && eIdx >= 1 && eIdx <= nSeg && pistonMask(eIdx)
                        uw = v_piston * n_hat;
                    else
                        uw = [0, 0];
                    end
                    vrel = v - uw;
                    vnrel = dot(vrel, n_hat);
                    needReflect = (~inside(i)) || (dist < (RpEff - tol) && vnrel < 0);
                    if needReflect
                        Xp(i) = cp(1) + n_hat(1) * (RpEff + tol);
                        Yp(i) = cp(2) + n_hat(2) * (RpEff + tol);
                        vt_rel = vrel - vnrel * n_hat;
                        vrel_new = vt_rel - eRest * vnrel * n_hat;
                        v_new = uw + vrel_new;
                        % Tangential damping
                        t_hat = [-n_hat(2), n_hat(1)];
                        vt_comp = dot(v_new, t_hat);
                        v_new = v_new - betaT * vt_comp * t_hat;
                        Up(i) = v_new(1); Vp(i) = v_new(2);
                    end
                end
            end
        end

        % Thermal velocity cap after this substep
        sp = hypot(Up, Vp);
        over = sp > vCap & isfinite(sp);
        if any(over)
            scaleOver = vCap ./ max(sp(over), eps);
            Up(over) = Up(over) .* scaleOver;
            Vp(over) = Vp(over) .* scaleOver;
        end
    end

    % Final sanitization: replace non-finite states safely inside the cavity
    bad = ~isfinite(Xp) | ~isfinite(Yp) | ~isfinite(Up) | ~isfinite(Vp);
    if any(bad)
        [cx, cy] = polygonCentroidSafe(Xcv, Ycv);
        for i = find(bad)'
            % Place near centroid with small random offset inside
            ang = 2 * pi * rand;
            rad = 0.2 * Lchar;
            Xp(i) = cx + rad * cos(ang);
            Yp(i) = cy + rad * sin(ang);
            % Project to interior if needed
            [cp, n_in, dist] = closestPointAndNormalToPolygon(Xp(i), Yp(i), Xcv, Ycv);
            if ~inpolygon(Xp(i), Yp(i), Xcv, Ycv) || dist < RpEff
                Xp(i) = cp(1) + n_in(1) * (RpEff + tol);
                Yp(i) = cp(2) + n_in(2) * (RpEff + tol);
            end
            % Assign bounded random velocity
            vmag = 0.5 * vCap;
            angv = 2 * pi * rand;
            Up(i) = vmag * cos(angv);
            Vp(i) = vmag * sin(angv);
        end
    end
end

function val = getOpt(s, name, def)
    if isfield(s, name) && ~isempty(s.(name))
        val = s.(name);
    else
        val = def;
    end
end

function r = frac(x)
    r = x - floor(x);
end

function [cp, n_in, dist, bestIdx] = closestPointAndNormalToPolygon(x, y, Xcv, Ycv)
% Returns:
%   cp   - closest point on polygon edges to (x,y)
%   n_in - interior unit normal (points toward interior of polygon)
%   dist - distance from (x,y) to cp
%   bestIdx - index of closest segment
    p = [x; y];
    bestD = inf;
    bestCP = [NaN; NaN];
    bestIdx = 1;
    for i = 1:numel(Xcv)-1
        a = [Xcv(i);   Ycv(i)];
        b = [Xcv(i+1); Ycv(i+1)];
        [cpi, ~] = closestPointOnSegment(a, b, p);
        d = hypot(p(1) - cpi(1), p(2) - cpi(2));
        if d < bestD
            bestD = d;
            bestCP = cpi;
            bestIdx = i;
        end
    end
    cp = bestCP;
    dist = bestD;
    % Interior normal from the closest segment
    a = [Xcv(bestIdx);   Ycv(bestIdx)];
    b = [Xcv(bestIdx+1); Ycv(bestIdx+1)];
    t = b - a;
    nt = hypot(t(1), t(2));
    if nt <= eps
        n1 = [1; 0];
    else
        t = t / nt;
        n1 = [-t(2); t(1)]; % one of the two normals
    end
    % Determine which normal points inwards by testing a small offset
    eps_n = 1e-6;
    c1 = cp + eps_n * n1;
    c2 = cp - eps_n * n1;
    in1 = inpolygon(c1(1), c1(2), Xcv, Ycv);
    in2 = inpolygon(c2(1), c2(2), Xcv, Ycv);
    if in1 && ~in2
        n_in = n1;
    elseif ~in1 && in2
        n_in = -n1;
    else
        % Fallback: direct away from edge towards the point
        v = p - cp;
        nv = hypot(v(1), v(2));
        if nv <= eps
            n_in = n1;
        else
            n_in = v / nv;
        end
    end
end

function [cp, t] = closestPointOnSegment(a, b, p)
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

function [cx, cy] = polygonCentroidSafe(X, Y)
    x = X(:);
    y = Y(:);
    if x(1) ~= x(end) || y(1) ~= y(end)
        x = [x; x(1)];
        y = [y; y(1)];
    end
    A2 = 0;
    Cx = 0;
    Cy = 0;
    for i = 1:numel(x)-1
        cross = x(i) * y(i+1) - x(i+1) * y(i);
        A2 = A2 + cross;
        Cx = Cx + (x(i) + x(i+1)) * cross;
        Cy = Cy + (y(i) + y(i+1)) * cross;
    end
    if abs(A2) < eps
        cx = mean(x);
        cy = mean(y);
        return;
    end
    A = A2 / 2;
    cx = Cx / (6 * A);
    cy = Cy / (6 * A);
end