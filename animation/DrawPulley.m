function [] = DrawPulley(p1, p2, R1, R2, face_color, LW, np)
% x1, y1 = coordinates starting point
% x2, y2 = coordinates end points
% t = thickness of the line
% np = number of points per unit

if nargin < 7; np = 100; end
if nargin < 5; face_color = 'w'; end
if nargin < 6; LW = 1; end

% Get x and y
x1 = p1(1); y1 = p1(2);
x2 = p2(1); y2 = p2(2);

% Length
L = sqrt((x2-x1)^2 + (y2-y1)^2);

% Angle theta
theta = 2 * acos((R1 - R2)/L);

% POINTS
% -------------------------------------------------------------------------

c1x = R1 * cos(theta/2);
c1y = R1 * sin(theta/2);

c2x = L + R2 * cos(theta/2);
c2y =     R2 * sin(theta/2);

c3x = c2x;
c3y =-c2y;

c4x = c1x;
c4y =-c1y;

% DRAW LINE HORIZONTALLY
% -------------------------------------------------------------------------

% Upper horizontal line from 0 to L
x_hu = linspace(c1x,c2x,np);
y_hu = linspace(c1y,c2y,np);

% Right circular line
ang = fliplr(linspace(-theta/2,theta/2,np));
x_cr = L + R2 * cos(ang);
y_cr = R2 * sin(ang);

% Lower horizontal line from L to 0
x_hl = linspace(c3x,c4x,np);
y_hl = linspace(c3y,c4y,np);

% Left circular line
ang = linspace(-theta/2,theta/2,np);
x_cl = - R1 * cos(ang);
y_cl = R1 * sin(ang);

% Assemble
x_h = [x_hu x_cr x_hl x_cl x_hu(1)];
y_h = [y_hu y_cr y_hl y_cl y_hu(1)];
xy_h = [x_h;y_h];

% ROTATE AND TRANSLATE LINE
% -------------------------------------------------------------------------

% Angle
theta = atan2((y2-y1),(x2-x1));

% Rotation matrix
R = [cos(theta) -sin(theta); sin(theta) cos(theta)];

% Rotate
xy_r = R * xy_h;
x_r = xy_r(1,:);
y_r = xy_r(2,:);

% Translate
x = x1 + x_r;
y = y1 + y_r;

% PLOT
% -------------------------------------------------------------------------

% Plot line
fill(x,y,face_color, 'LineStyle','None');
plot(x,y,'Color', [0 0 0],'LineWidth',LW);



