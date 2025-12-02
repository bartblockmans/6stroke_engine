function [x, y, xi, yi, xe, ye] = CylinderCavity(pp, params)
%CYLINDER_CAVITY  Compute the cylinder cavity geometry.
%
%   [x, y, xi, yi, xe, ye] = CylinderCavity(pp, params)
%
%   Inputs:
%       pp : piston position [mm]
%       params : structure with geometry parameters (from parameters.m)
%
%   Outputs:
% -------------------------------------------------------------------------
% Bart Blockmans, 2024 - bart@blockmans.net
% -------------------------------------------------------------------------

% Get parameters
H = params.H;
B = params.B;
Y_max = params.Ys_max + (2/3) * params.H + params.Hc;
Rc = params.Rc;

% Cylinder cavity
% -------------------------------------------------------------------------

% Height of the upper piston surface
yp = pp(2) + (2/3) * H;

% Top coordinates
xp = linspace(-B/2, B/2, 10);
yp = ones(1,10) * yp;

% Upper surface
ang_u = asin((B/2)/ Rc);
ang = linspace(pi/2 - ang_u, pi/2 + ang_u, 100);
xu = Rc * cos(ang);
yu = Y_max - Rc * cos(ang_u) + Rc * sin(ang);

% Left surface
xl = ones(1,10) * (-B/2);
yl = linspace(yp(1),yu(1),10);

% Right surface
xr = ones(1,10) * (B/2);
yr = linspace(yp(end),yu(end),10);

% Combine
x = [xl, fliplr(xu), fliplr(xr), fliplr(xp)];
y = [yl, fliplr(yu), fliplr(yr), fliplr(yp)];

% Plot
% plot(x,y,'r','LineWidth',2);

% Intake & exhaust
% -------------------------------------------------------------------------

% Angular span intake & exhaust
dang = abs(ang_u/3);

% Angular location intake & exhaust
ang_ie = (2/3)*abs(ang_u);

% Intake coordinates
ang_i1 = pi/2 + ang_ie + dang/2;
ang_i2 = pi/2 + ang_ie - dang/2;
xi = Rc * cos([ang_i1, ang_i2]);
yi = Y_max - Rc * cos(ang_u) + Rc * sin([ang_i1, ang_i2]);
% scatter(xi,yi,'filled');

% Exhaust coordinates
ang_e1 = pi/2 - (ang_ie + dang/2);
ang_e2 = pi/2 - (ang_ie - dang/2);
xe = Rc * cos([ang_e1, ang_e2]);
ye = Y_max - Rc * cos(ang_u) + Rc * sin([ang_e1, ang_e2]);
% scatter(xe,ye,'filled');





