function h = ColorCavity(Xc, Yc, T, X_lim)
%COLORCYLINDER  Color the interior of a cylinder cross-section by temperature.
%
%   h = ColorCavity(Xc, Yc, T)
%
%   Inputs
%   ------
%   Xc, Yc : Vectors describing the circumference of the cylinder cavity (2D curve).
%   T      : Dimensionless temperature parameter in [0, 1]
%            T = 0 → minimal temperature (cold, blue)
%            T = 1 → maximal temperature (hot, red)
%
%   Output
%   ------
%   h      : Handle to the patch object used to fill the cylinder.
%
%   Example (inside animation loop)
%   --------------------------------
%   % Xc, Yc: precomputed cylinder boundary
%   % T: temperature parameter updated each frame
%   ColorCylinder(Xc, Yc, T);
%   axis equal; drawnow;

% Basic input check
if nargin < 3
    error('ColorCylinder requires Xc, Yc, and T as inputs.');
end

% Ensure column vectors
Xc = Xc(:);
Yc = Yc(:);

if numel(Xc) ~= numel(Yc)
    error('Xc and Yc must have the same length.');
end

% Clip cylinder cavity to X_lim [-X_lim, X_lim]
% while ensuring the polygon is closed
[Xc, Yc] = clipPolygonToXRange(Xc, Yc, -X_lim, X_lim);

% Clamp T to [0, 1] to avoid out-of-range issues
T = max(0, min(1, T));

% ---------------------------------------------------------------------
% Define a private colormap: blue (cold) → red (hot)
% ---------------------------------------------------------------------
nColors = 256;
cmap = zeros(nColors, 3);     % [R G B]
cmap(:,1) = linspace(0, 1, nColors); % R: 0 → 1
cmap(:,2) = 0;                      % G: 0
cmap(:,3) = linspace(1, 0, nColors); % B: 1 → 0

% Convert T ∈ [0,1] to a color index
idx = 1 + round(T * (nColors - 1));
idx = max(1, min(nColors, idx));
colorT = cmap(idx, :);

% ---------------------------------------------------------------------
% Draw the filled cylinder cross-section
% ---------------------------------------------------------------------
h = fill(Xc, Yc, colorT, ...
         'EdgeColor', 'none', ...   
         'FaceColor', colorT, ...
         'FaceAlpha',0.30);

% For a realistic view, you typically want:
% axis equal;  % but leave it to the caller so we don't override their settings

% Nested function: Clip polygon to X range
function [Xo, Yo] = clipPolygonToXRange(Xi, Yi, x_min, x_max)
    n = numel(Xi);
    if n == 0
        Xo = Xi; Yo = Yi;
        return;
    end
    
    % Quick check: if all points inside, return as-is
    if all(Xi >= x_min & Xi <= x_max)
        Xo = Xi; Yo = Yi;
        return;
    end
    
    Xo = [];
    Yo = [];
    
    % Process each edge
    for i = 1:n
        j = mod(i, n) + 1;
        x1 = Xi(i); y1 = Yi(i);
        x2 = Xi(j); y2 = Yi(j);
        
        % Check if point i is inside
        if x1 >= x_min && x1 <= x_max
            Xo = [Xo; x1];
            Yo = [Yo; y1];
        end
        
        % Check edge intersections with boundaries
        if x1 ~= x2
            % Left boundary
            if (x1 < x_min && x2 >= x_min) || (x1 >= x_min && x2 < x_min)
                t = (x_min - x1) / (x2 - x1);
                if t > 0 && t < 1
                    Xo = [Xo; x_min];
                    Yo = [Yo; y1 + t * (y2 - y1)];
                end
            end
            % Right boundary
            if (x1 > x_max && x2 <= x_max) || (x1 <= x_max && x2 > x_max)
                t = (x_max - x1) / (x2 - x1);
                if t > 0 && t < 1
                    Xo = [Xo; x_max];
                    Yo = [Yo; y1 + t * (y2 - y1)];
                end
            end
        end
    end
    
    % Ensure closed polygon
    if ~isempty(Xo) && (Xo(1) ~= Xo(end) || Yo(1) ~= Yo(end))
        Xo = [Xo; Xo(1)];
        Yo = [Yo; Yo(1)];
    end
end

end
