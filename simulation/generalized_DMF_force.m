function Q_dmf = generalized_DMF_force(params, theta_c, omega_c, theta_o, omega_o)
%GENERALIZED_dmf_FORCE  Generalized torques from DMF spring-damper.
%
%   Q_dmf = generalized_dmf_force(params, theta_c, omega_c, theta_o, omega_o)
%
%   Generalized coordinates:
%     q1 = theta_p  [rad]  planet gear angle
%     q2 = theta_c  [rad]  carrier / crank angle  (inner DMF side)
%     q3 = theta_o  [rad]  output shaft angle     (outer DMF side)
%
%   Inputs:
%     params   : struct with fields
%                  .k_dmf     [N·m/rad]     DMF torsional stiffness
%                  .c_dmf     [N·m·s/rad]   DMF torsional damping
%                 (optional) .phi_dmf0 [rad] equilibrium twist (preload)
%
%     theta_c  : [rad]      carrier angle
%     omega_c  : [rad/s]    carrier angular velocity
%     theta_o  : [rad]      output shaft angle
%     omega_o  : [rad/s]    output shaft angular velocity
%
%   Output:
%     Q_dmf   : [3x1] generalized DMF torques:
%                Q_dmf(1) = 0
%                Q_dmf(2) = Q_theta_c^DMF
%                Q_dmf(3) = Q_theta_o^DMF
%
%   Theory:
%     Relative twist and rate:
%       phi     = theta_c - theta_o
%       phi_dot = omega_c - omega_o
%
%     DMF torque:
%       T_dmf = k_dmf * (phi - phi0) + c_dmf * phi_dot
%
%     Generalized torques:
%       Q_theta_c = -T_dmf
%       Q_theta_o = +T_dmf
%       Q_theta_p =  0

    % Unpack parameters
    k_dmf = params.k_dmf;
    c_dmf = params.c_dmf;

    if isfield(params, 'phi_dmf0')
        phi0 = params.phi_dmf0;
    else
        phi0 = 0.0;
    end

    % Relative twist and rate
    phi     = theta_c - theta_o;
    phi_dot = omega_c - omega_o;

    % DMF torque (positive from output on carrier if phi > phi0)
    T_dmf = k_dmf .* (phi - phi0) + c_dmf .* phi_dot;

    % Generalized forces
    Q_theta_p = 0;           % no direct action on planet DOF
    Q_theta_c = -T_dmf;      % torque on carrier
    Q_theta_o =  T_dmf;      % opposite torque on output

    % Assemble vector
    Q_dmf = [Q_theta_p; Q_theta_c; Q_theta_o];
end
