function params = manifold_geometry(params)
%MANIFOLD_GEOMETRY  Compute manifold geometry for 6-stroke engine.
%
%   params = manifold_geometry(params)
%
%   Inputs:
%       params      : structure with geometry parameters (from parameters.m)
%
%   Outputs:
%       params      : modified structure with added fields:
% -------------------------------------------------------------------------
% Bart Blockmans, 2024 - bart@blockmans.net
% -------------------------------------------------------------------------

B = params.B;
Do = params.Do;
to = params.to;
tw = params.tw;
Y_max = params.Ys_max + 2*params.H/3 + params.Hc;
Rc = params.Rc;
Y_bdc2 = params.BDC2_y + 2*params.H/3;
Y_bdc1 = params.BDC1_y + 2*params.H/3;
np = params.particles.np;
np_scav = params.particles.np_scav;

% Crank part of the cylinder
% -------------------------------------------------------------------------

% Effective outside diameter
D = Do + 2 * to;

% Limit angles outer diameter
ang_o = asin((B/2 + tw) / (D/2));

% Outer circles
x_co = (D/2) * cos(pi/2 + ang_o);

% Top of the cylinder
% -------------------------------------------------------------------------

% Angle
ang_ti = asin((B/2)/ Rc);

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

% Store in params
params.xmi1 = xmi1;
params.ymi1 = ymi1;
params.xmi2 = xmi2;
params.ymi2 = ymi2;
params.xme1 = -xmi2;
params.yme1 =  ymi2;
params.xme2 = -xmi1;
params.yme2 =  ymi1;

% Get port coordinates
params.xip = [xmi2(1); xmi1(1)];
params.yip = [ymi2(1); ymi1(1)];
params.xep = -[xmi2(1); xmi1(1)];
params.yep =  [ymi2(1); ymi1(1)];

% Scavenge ports
% -------------------------------------------------------------------------

% Angle between TDC1 and BDC2
theta_deg = linspace(params.TDC1_th(1), params.BDC2_th,1000);

% Compute engine kinematics
[~, ~, x_piston] = geometry_6stroke(theta_deg, params);

% Find angle where 
[~, ind] = min(abs(x_piston - params.BDC1_x(1)));

% Start angle opening scavenge ports
params.SP_open = theta_deg(ind);

% Angle between BDC2 and TDC1(2)
theta_deg = linspace(params.BDC2_th,params.TDC1_th(2),1000);

% Compute engine kinematics
[~, ~, x_piston] = geometry_6stroke(theta_deg, params);

% Find angle where 
[~, ind] = min(abs(x_piston - params.BDC1_x(1)));

% Start angle opening scavenge ports
params.SP_close = theta_deg(ind);

% Max effective port height [m]
% params.SP_hmax = (Y_bdc1 - Y_bdc2) * 1e-3;

% Get coordinates
xms1 = [-0.5*B, -0.5*B-tw];
xms2 = [-0.5*B, -0.5*B-tw];
yms1 = [Y_bdc2, Y_bdc2];
yms2 = [Y_bdc1, Y_bdc1];

% Store
params.xmsl1 = xms1;
params.xmsl2 = xms2;
params.ymsl1 = yms1;
params.ymsl2 = yms2;
params.xmsr1 = -xms1;
params.xmsr2 = -xms2;
params.ymsr1 = yms1;
params.ymsr2 = yms2;

% Reservoirs
% -------------------------------------------------------------------------

% Cavity boundary & area at theta_crank = 0
[~, ~, ~, ~, ~, ~, coords] = geometry_6stroke(0, params);
[Xc, Yc] = CylinderCavity([coords.X_s, coords.Y_s], params);
A0  = polyarea(Xc, Yc);

% Imaginary intake reservoir
Hr = ymi2(end) - ymi1(end);
Br = 0.5 * A0 / Hr;
xir = [xmi1(end), xmi1(end)-Br, xmi2(end)-Br, xmi2(end)];
yir = [ymi1(end), ymi1(end), ymi2(end), ymi2(end)];

% Store
params.xir = xir;
params.yir = yir;

% Imaginary exhaust reservoir
params.xer = -xir;
params.yer = yir;

% Imaginary left scavenge port reservoir
Hrs = Y_bdc1 - Y_bdc2;
Brs = 2 * (np_scav/np) * A0 / Hrs;
xsrl = [-0.5*B-tw, -0.5*B-tw-Brs, -0.5*B-tw-Brs, -0.5*B-tw];
ysrl = [Y_bdc2, Y_bdc2, Y_bdc1, Y_bdc1];

% Store
params.xsrl = xsrl;
params.ysrl = ysrl;

% Imaginary right scavenge port reservoir
params.xsrr = -xsrl;
params.ysrr = ysrl;

% Valve geometry
% -------------------------------------------------------------------------

% Get parameters
Xip = params.xip;
Yip = params.yip;
Xep = params.xep;
Yep = params.yep;

% Valve diameter
D = sqrt((Xip(2)-Xip(1))^2 + (Yip(2)-Yip(1))^2);

% Opening displacement
d = D/2;

% Valve angle
theta_in = atan2(Yip(1)-Yip(2), Xip(1)-Xip(2));
theta_ex = atan2(Yep(2)-Yep(1), Xep(2)-Xep(1));

% Valve length
L = 2.75 * D;

% Cam center location
x_cam_IV = mean(Xip) - (L + D/2) * sin(theta_in);
y_cam_IV = mean(Yip) + (L + D/2) * cos(theta_in);
x_cam_EV = mean(Xep) - (L + D/2) * sin(theta_ex);
y_cam_EV = mean(Yep) + (L + D/2) * cos(theta_ex);

% Store
params.D_valve = D;
params.valve_disp = d;
params.theta_IV = theta_in;
params.theta_EV = theta_ex;
params.L_valve = L;
params.x_cam_IV = x_cam_IV;
params.y_cam_IV = y_cam_IV;
params.x_cam_EV = x_cam_EV;
params.y_cam_EV = y_cam_EV;

% Valve lift for animation
params = ValveLiftForAnimation(params);