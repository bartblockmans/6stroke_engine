function [x, dxdt, ddxdt, Ac, Ap, Bcc, Bpp, Bcp] = piston_kinematics( ...
    params, theta_c, omega_c, alpha_c, theta_p, omega_p, alpha_p)
%PISTON_KINEMATICS  Piston kinematics for 6-stroke internal-gear slider.
%
%   [x, dxdt, ddxdt, Ac, Ap, Bcc, Bpp, Bcp] = piston_kinematics( ...
%       params, theta_c, omega_c, alpha_c, theta_p, omega_p, alpha_p)
%
%   Geometry is made consistent with geometry_6stroke.m:
%     - theta_c: carrier/crank angle [rad], same as deg2rad(theta)
%     - Internally we use the same orbit of the planet center as:
%         phi = theta_c - 3*pi/2
%       which is equivalent to rotating the carrier frame by +90°.
%     - theta_p: angle of the line from planet center to rod attachment
%       (red point) w.r.t. the global X-axis [rad]; this is the same
%       quantity called "h" in geometry_6stroke when psi0 and the gear
%       relation are enforced.
%
%   Inputs:
%     params : struct with fields
%                .L     [m]   rod length
%                .delta [m]   eccentricity of red point from planet center
%                .a     [m]   orbit radius of planet center (= Rr - Rp)
%     theta_c : carrier/crank angle [rad]
%     omega_c : carrier angular velocity [rad/s]
%     alpha_c : carrier angular acceleration [rad/s^2]
%     theta_p : planet eccentric-point angle [rad]
%     omega_p : planet eccentric-point angular velocity [rad/s]
%     alpha_p : planet eccentric-point angular acceleration [rad/s^2]
%
%   Outputs:
%     x      : piston vertical position (Y_s) [m] w.r.t. crank axis
%     dxdt   : piston vertical velocity [m/s]
%     ddxdt  : piston vertical acceleration [m/s^2]
%     Ac     : ∂x/∂theta_c
%     Ap     : ∂x/∂theta_p
%     Bcc    : ∂²x/∂theta_c²
%     Bpp    : ∂²x/∂theta_p²
%     Bcp    : ∂²x/(∂theta_c∂theta_p)
%
%   Note:
%     geometry_6stroke.m uses x_piston = Y_s_max - Y_s (displacement from
%     highest TDC). Here we use x = Y_s itself. This only differs by a
%     constant offset; all derivatives (and thus forces/torques) remain
%     consistent. If you need x_piston, you can subtract params.Ys_max
%     externally if you have it.
% -------------------------------------------------------------------------
% Bart Blockmans, 2024 - bart@blockmans.net
% -------------------------------------------------------------------------

% Unpack parameters
L     = params.L;
delta = params.delta;
a     = params.a;

% Trig shortcuts
sc = sin(theta_c);
cc = cos(theta_c);
sp = sin(theta_p);
cp = cos(theta_p);

% Eccentric point coordinates (consistent with geometry_6stroke)
% Planet center: X_p = -a*sin(theta_c), Y_p = a*cos(theta_c)
Xr = -a.*sc + delta.*cp;
Yr =  a.*cc + delta.*sp;

% Slider position on Y-axis: (0, Y_s), rod length L
arg = L.^2 - Xr.^2;
if any(arg(:) < 0)
    % Numerical safety; you can choose to warn or clip
    warning('piston_kinematics:RodTooShort', ...
        'Some L^2 - Xr^2 < 0; clipping to zero for sqrt.');
    arg = max(arg, 0);
end
D = sqrt(arg);              % = sqrt(L^2 - Xr^2)
D_safe = D;
eps_val = eps;
D_safe(D_safe < eps_val) = eps_val;

% Piston vertical position (same Y_s as in geometry_6stroke)
x = Yr + D;

% --- First derivatives wrt angles (Jacobian) ------------------------------

% ∂x/∂theta_c
Ac = -a.*sc + (a.*Xr.*cc)./D_safe;

% ∂x/∂theta_p
Ap =  delta.*cp + (delta.*Xr.*sp)./D_safe;

% --- Velocity --------------------------------------------------------------
dxdt = Ac.*omega_c + Ap.*omega_p;

% --- Second derivatives wrt angles ---------------------------------------- 
% ddx/dt^2 = Bcc*omega_c^2 + Bpp*omega_p^2 + 2*Bcp*omega_c*omega_p ...
%            + Ac*alpha_c + Ap*alpha_p

% Bcc = ∂²x/∂theta_c²
Bcc = -a.*cc ...
    + a.*(-a.*cc.^2 - Xr.*sc)./D_safe ...
    - (a.^2 .* Xr.^2 .* cc.^2)./(D_safe.^3);

% Bpp = ∂²x/∂theta_p²
Bpp = -delta.*sp ...
    + delta.*(-delta.*sp.^2 + Xr.*cp)./D_safe ...
    - (delta.^2 .* Xr.^2 .* sp.^2)./(D_safe.^3);

% Bcp = ∂²x/(∂theta_c∂theta_p)
Bcp = a .* ( (-delta.*sp.*cc)./D_safe ...
           - (delta.*Xr.^2.*cc.*sp)./(D_safe.^3) );

% --- Acceleration ---------------------------------------------------------
ddxdt = Bcc.*omega_c.^2 + Bpp.*omega_p.^2 + 2.*Bcp.*omega_c.*omega_p ...
      + Ac.*alpha_c + Ap.*alpha_p;

end
