function [] = DrawValves(IV, EV, params)

% Get parameters
Xip = params.xip;
Yip = params.yip;
Xep = params.xep;
Yep = params.yep;
colors = params.colors;
LW = params.LW;

% Valve diameter
D = sqrt((Xip(2)-Xip(1))^2 + (Yip(2)-Yip(1))^2);

% Opening displacement
d = D/2;

% Valve length
L_valve = 2.75 * D;

% Valve head
x_head = [-D/2, -D/2 + D/6, D/2 - D/6, D/2, -D/2];
y_head = [0, D/6, D/6, 0, 0];

% Little circular corner
ang = linspace(-pi/2,0,20);
x_corner = (D/5) * cos(ang);
y_corner = (D/5) * sin(ang);
x_corner = x_corner - D/5;
y_corner = y_corner + D/5;

% Coordinates
x_rod = [x_corner - D/10, -D/10, D/10, -fliplr(x_corner) + D/10, x_corner(1) - D/10];
y_rod = [y_corner + D/6, L_valve, L_valve, fliplr(y_corner) + D/6, y_corner(1) + D/6];

% Rotate and translate for intake
% -------------------------------------------------------------------------

% Intake angle
theta_in = atan2(Yip(1)-Yip(2), Xip(1)-Xip(2));

% Rotation matrix
R = [cos(theta_in), -sin(theta_in); sin(theta_in), cos(theta_in)];

% Location
x_intake = mean(Xip);
y_intake = mean(Yip);

% Rotate & translate
xy_intake_head = [x_intake, y_intake] + [x_head', y_head'] * R';
x_intake_head = xy_intake_head(:,1);
y_intake_head = xy_intake_head(:,2);
xy_intake_rod = [x_intake, y_intake] + [x_rod', y_rod'] * R';
x_intake_rod = xy_intake_rod(:,1);
y_intake_rod = xy_intake_rod(:,2);

% Exhaust valve
% -------------------------------------------------------------------------

x_exhaust_rod = -x_intake_rod;
y_exhaust_rod =  y_intake_rod;
x_exhaust_head = -x_intake_head;
y_exhaust_head =  y_intake_head;

% Open valves if needed
% -------------------------------------------------------------------------

% Open input valve if needed:
if IV > 0
    dx_in =  IV * d * sin(theta_in);
    dy_in = -IV * d * cos(theta_in);

else; dx_in = 0; dy_in = 0; 
end

if EV > 0
    % Exhasut angle
    theta_ex = atan2(Yep(2)-Yep(1), Xep(2)-Xep(1));

    dx_ex =  EV * d * sin(theta_ex);
    dy_ex = -EV * d * cos(theta_ex);

else; dx_ex = 0; dy_ex = 0; 
end

% Draw  valve
fill(dx_in + x_intake_head, dy_in + y_intake_head, colors.valve,'EdgeColor',[0 0 0], 'LineWidth',LW);
fill(dx_in + x_intake_rod, dy_in + y_intake_rod, colors.valve,'EdgeColor',[0 0 0], 'LineWidth',LW);
fill(dx_ex + x_exhaust_head, dy_ex + y_exhaust_head, colors.valve,'EdgeColor',[0 0 0], 'LineWidth',LW);
fill(dx_ex + x_exhaust_rod, dy_ex + y_exhaust_rod, colors.valve,'EdgeColor',[0 0 0], 'LineWidth',LW);