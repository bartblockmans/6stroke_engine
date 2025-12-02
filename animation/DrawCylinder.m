function [xmi1, ymi1, xmi2, ymi2] = DrawCylinder(params)
%DRAW_CYLINDER  Draw the cylinder for the 6-stroke engine.
%
%   [xmi1, ymi1, xmi2, ymi2] = DrawCylinder(params)
%
%   Inputs:
%       params : structure with geometry parameters (from parameters.m)
%
%   Outputs:
% -------------------------------------------------------------------------
% Bart Blockmans, 2024 - bart@blockmans.net
% -------------------------------------------------------------------------

B = params.B;
Do = params.Do;
to = params.to;
tw = params.tw;
Y_max = params.Ys_max + 2*params.H/3 + params.Hc;
Rc = params.Rc;
colors = params.colors;
LW = params.LW;

% Crank part of the cylinder
% -------------------------------------------------------------------------

% Effective outside diameter
D = Do + 2 * to;

% Limit angles outer diameter
ang_i = asin((B/2) / (D/2));
ang_o = asin((B/2 + tw) / (D/2));

% Outer circles
ang = linspace(pi/2 - ang_i, pi/2 + ang_i, 100);
x_ci = (D/2) * cos(ang);
y_ci = (D/2) * sin(ang);
ang = linspace(pi/2 + ang_o, pi/2 + 2 * pi - ang_o, 500);
x_co = (D/2) * cos(ang);
y_co = (D/2) * sin(ang);

% Inner circle
x_ic = (Do/2) * cos(linspace(0,2*pi,200));
y_ic = (Do/2) * sin(linspace(0,2*pi,200));

% Top of the cylinder
% -------------------------------------------------------------------------

% Angle
ang_ti = asin((B/2)/ Rc);

% Outer circles
ang = linspace(pi/2 - ang_ti, pi/2 + ang_ti, 100);
x_i = Rc * cos(ang);
y_i = Y_max - Rc * cos(ang_ti) + Rc * sin(ang);

% Intake & exhaust
% -------------------------------------------------------------------------

% Angular span intake & exhaust
dang = abs(ang_ti/2.5);

% Angular location intake & exhaust
ang_ie = (2/3)*abs(ang_ti);

% Intake coordinates
ang_i1 = pi/2 + ang_ie + dang/2;
ang_i2 = pi/2 + ang_ie - dang/2;
xii = Rc * cos([ang_i1, ang_i2]);
yii = Y_max - Rc * cos(ang_ti) + Rc * sin([ang_i1, ang_i2]);

% Radius of the intake/exhaust manifold
Rm = (xii(1) - x_co(1)) / cos(pi/6);

% Center of rotation
xc = xii(1) - Rm * cos(pi/6);
yc = yii(1) - Rm * sin(pi/6);

% Intake
ang = linspace(pi/6, pi/2, 100);
xmi1 = xc + Rm * cos(ang);
ymi1 = yc + Rm * sin(ang);
dist = sqrt((xii(2)-xc)^2 + (yii(2)-yc)^2);
xmi2 = xc + dist * cos(ang);
ymi2 = yc + dist * sin(ang);
xmi2 = [xii(2), xmi2];
ymi2 = [yii(2), ymi2];
x_in = [xmi1, fliplr(xmi2), xmi1(1)];
y_in = [ymi1, fliplr(ymi2), ymi1(1)];

% Outer lines
% -------------------------------------------------------------------------

% Margin above manifold
dman = 10;

% Bougie height
Yb = 0.65 * 20 + max(y_i);

% Angle
angb = asin((Yb - yc) / (dist + dman));

% ang
ang = linspace(angb, pi/2, 50);
x_top = xc + (dist + dman) * cos(ang);
y_top = yc + (dist + dman) * sin(ang);
x_top = [fliplr(x_top), -x_top];
y_top = [fliplr(y_top),  y_top];

% Fill cylinder
% -------------------------------------------------------------------------

% Construct x & y
x_outer = [fliplr(x_co), x_top, x_co(end)];
y_outer = [fliplr(y_co), y_top, y_co(end)];

% Inner curve
x_inner = [x_ci, fliplr(x_i), x_ci(1)];
y_inner = [y_ci, fliplr(y_i), y_ci(1)];

% Combine
x_fill = [x_outer x_inner];
y_fill = [y_outer y_inner];

% Fill
fill(x_fill, y_fill, colors.cyl, 'EdgeColor','none');

% Plot circles around crank shaft
plot(x_ci, y_ci, 'k','LineWidth',LW);
plot(x_co, y_co, 'k','LineWidth',LW);

% Plot top of cylinder
plot(x_i, y_i, 'k','LineWidth',LW);

% Plot vertical lines
plot(x_ci(1)*ones(2,1), [y_ci(1) y_i(1)],'k','LineWidth',LW);
plot(x_ci(end)*ones(2,1), [y_ci(1) y_i(1)],'k','LineWidth',LW);
plot(x_co(1)*ones(2,1), [y_co(1) ymi2(end) + dman],'k', 'LineWidth',LW);
plot(-x_co(1)*ones(2,1), [y_co(1) ymi2(end) + dman],'k', 'LineWidth',LW);

% Plot top of manifold
plot(x_top, y_top, 'k', 'LineWidth',LW);

% Plot inner circle
fill(x_ic, y_ic, [1 1 1],'EdgeColor',[0 0 0],'LineWidth',LW);

% Fill chamber
fill(x_inner, y_inner, colors.cham,'EdgeColor',[0 0 0],'LineWidth',LW);

% Fill manifolds
% -------------------------------------------------------------------------

% Fill manifolds
fill(x_in, y_in, colors.man, 'EdgeColor',[0 0 0],'LineWidth',LW)

% Exhaust
fill(-x_in, y_in, colors.man, 'EdgeColor',[0 0 0],'LineWidth',LW)






