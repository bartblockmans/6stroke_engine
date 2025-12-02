function [] = DrawCrank(p, params)

% Get parameters
Ds = params.Ds;
Rw = params.Rcw;
colors = params.colors;
LW = params.LW;

% Angle
theta_c = atan2(p(2), p(1)) - pi/2;

% Rotation matrix
R = [cos(theta_c) -sin(theta_c); sin(theta_c), cos(theta_c)];

% Length
L = sqrt(p(1)^2 + p(2)^2);

% Top circle
ang = linspace(0, pi, 100);
x_top = (Ds/2) * cos(ang);
y_top = L + (Ds/2) * sin(ang);
x_top = -x_top;

% Counterweight
ang = linspace(pi/6, pi - pi/6, 100);
x_cw = Rw * cos(ang);
y_cw = Rw * sin(ang);
x_cw = -x_cw;
y_cw = -y_cw;

% Assemble whole crankshaft
x = [x_top, x_top(end), fliplr(x_cw), x_top(1), x_top(1)];
y = [y_top, 0         , fliplr(y_cw), 0       , y_top(1)];

% Rotate
xy = R * [x; y];
x = xy(1,:);
y = xy(2,:);

% Plot
fill(x,y,colors.cw,'EdgeColor',[0 0 0],'LineWidth',LW);
% plot(x,y,'k');
