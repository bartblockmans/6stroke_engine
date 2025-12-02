function Q_mesh = generalized_mesh_force(params, theta_p, omega_p, theta_c, omega_c)
%GENERALIZED_MESH_FORCE  Generalized torques from planet–ring gear mesh.
%
%   Q_mesh = generalized_mesh_force(params, theta_p, omega_p, theta_c, omega_c)
%
%   Generalized coordinates:
%     q1 = theta_p  [rad]  planet gear angle
%     q2 = theta_c  [rad]  carrier / crank angle
%     q3 = theta_o  [rad]  output shaft angle (not affected by mesh)
%
%   Inputs:
%     params   : struct with fields
%                  .Rp          [m]     planet pitch radius
%                  .Rr          [m]     ring pitch radius
%                  .k_mesh_avg  [N/m]   average mesh stiffness
%                  .k_mesh_amp  [N/m]   amplitude of 1st harmonic variation
%                  .zr          [-]     # teeth of ring gear
%                  .c_mesh      [N·s/m] mesh damping along mesh deflection y
%
%     theta_p  : [rad]      planet angle
%     omega_p  : [rad/s]    planet angular velocity
%     theta_c  : [rad]      carrier angle
%     omega_c  : [rad/s]    carrier angular velocity
%
%   Output:
%     Q_mesh   : [3x1] generalized torque vector (q1,q2,q3):
%                  Q_mesh(1) = Q_theta_p^mesh
%                  Q_mesh(2) = Q_theta_c^mesh
%                  Q_mesh(3) = 0 (mesh does not act directly on theta_o)
%
%   Theory (deflection & damping as in your existing implementation):
%     Mesh deflection coordinate:
%       y      = (R_p - R_r)*(theta_c - theta_c0) - R_p*(theta_p - theta_p0)
%       y_dot  = (R_p - R_r)*omega_c - R_p*omega_p
%
%     Mesh stiffness (1st harmonic vs carrier angle):
%       k_mesh(theta_c) = k_mesh_avg + k_mesh_amp * sin( z_r * (theta_c - theta_c0) )
%
%     Mesh force along y:
%       F_m = -k_mesh(theta_c) * y - c_mesh * y_dot
%
%     Generalized torques:
%       Q_theta_p =  F_m * (∂y/∂theta_p) =  F_m * (-R_p)
%       Q_theta_c =  F_m * (∂y/∂theta_c) =  F_m * (R_p - R_r)
%
%     or explicitly (same as your original, but with k_mesh(theta_c)):
%       Q_theta_p =  R_p        * (k_mesh*y + c_mesh*y_dot)
%       Q_theta_c = -(R_p - R_r)* (k_mesh*y + c_mesh*y_dot)

    % Unpack parameters
    R_p = params.Rp;
    R_r = params.Rr;
    k_avg = params.k_mesh_avg;           % average stiffness
    if isfield(params, 'k_mesh_amp') && ~isempty(params.k_mesh_amp)
        k_amp = params.k_mesh_amp;       % amplitude of variation
    else
        k_amp = 0.0;                     % default: no variation
    end
    c_mesh = params.c_mesh;
    zr     = params.zr;                  % # teeth ring gear

    % Reference angles used to define zero-deflection configuration
    theta_c0 = 3*pi/2;
    theta_p0 = pi/2;

    % Mesh deflection and rate (as in your working implementation)
    y     = (R_p - R_r) .* (theta_c - theta_c0) - R_p .* (theta_p - theta_p0);
    y_dot = (R_p - R_r) .* omega_c           - R_p .* omega_p;

    % Instantaneous mesh stiffness: 1st harmonic in carrier angle.
    % Period in theta_c = 2*pi / zr  => use sin(zr * (theta_c - theta_c0))
    k_mesh = k_avg + k_amp * sin( zr * (theta_c - theta_c0) );

    % Convenience term: k(theta_c)*y + c*y_dot
    KyCy = k_mesh .* y + c_mesh .* y_dot;

    % Generalized torques
    Q_theta_p =  R_p        .* KyCy;    % torque on planet DOF
    Q_theta_c = -(R_p - R_r).* KyCy;    % torque on carrier DOF
    Q_theta_o = 0;                      % no direct mesh torque on output

    % Assemble vector
    Q_mesh = [Q_theta_p; Q_theta_c; Q_theta_o];
end
