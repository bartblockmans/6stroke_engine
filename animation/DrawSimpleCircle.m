function [] = DrawSimpleCircle(p, R, face_color, LW, np)

% Check input arguments
if nargin < 3; face_color = [1 1 1]; end
if nargin < 4; LW = 1.5; end
if nargin < 5; np = 100; end

ang = linspace(0,2*pi,np);
x = p(1) + R * cos(ang);
y = p(2) + R * sin(ang);
fill(x,y,face_color, 'EdgeColor', [0 0 0], 'LineWidth', LW);