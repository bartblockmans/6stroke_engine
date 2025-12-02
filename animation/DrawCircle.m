function [] = DrawCircle(p, Do, Di, face_color, np, line_color)
% 
% Inputs:
%   p          - [x, y] coordinates of the center
%   R          - Outer radius of the circle
%   t          - Thickness of the ring (must be < R)
%   line_color - Color of the circle border lines (e.g., 'k' or [0 0 0])
%   face_color - Fill color of the thick ring (e.g., 'r' or [1 0 0])

if nargin < 4; face_color = 'w'; end
if nargin < 5; np = 100; end
if nargin < 6; line_color = 'k'; end

% Number of points to create smooth circle
theta = linspace(0, 2*pi, np);

% Outer and inner circle coordinates
x_outer = p(1) + 0.5 * Do * cos(theta);
y_outer = p(2) + 0.5 * Do * sin(theta);

x_inner = p(1) + 0.5 * Di * cos(flip(theta));  % Flip for proper fill
y_inner = p(2) + 0.5 * Di * sin(flip(theta));

% Combine outer and inner ring path for fill
x_ring = [x_outer, x_inner];
y_ring = [y_outer, y_inner];

% Fill the thick ring
fill(x_ring, y_ring, face_color, 'EdgeColor', 'none');

% Plot outer and inner border
plot(x_outer, y_outer, 'Color', line_color, 'LineWidth',1.5);
plot(x_inner, y_inner, 'Color', line_color, 'LineWidth',1.5);

end
