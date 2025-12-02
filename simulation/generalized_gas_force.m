function Q_gas = generalized_gas_force(params, theta_p, theta_c, p_cyl)
%GENERALIZED_GAS_FORCE  Generalized torques from cylinder gas pressure.
%
%   Q_gas = generalized_gas_force(params, theta_p, theta_c, p_cyl)
%
%   Generalized coordinates:
%     q1 = theta_p  [rad]  planet gear angle
%     q2 = theta_c  [rad]  carrier / crank angle
%     q3 = theta_o  [rad]  output shaft angle
%
%   Inputs:
%     params   : struct with at least:
%                  .L, .delta, .a   (for piston_kinematics)
%                and either:
%                  .A_piston        [m^2] piston area
%                or:
%                  .B               [m]   cylinder bore diameter
%
%     theta_p  : [rad]      planet angle
%     theta_c  : [rad]      carrier angle
%     p_cyl   : [Pa]       cylinder pressure (can be scalar or vector)
%
%   Output:
%     Q_gas   : [3x1] generalized gas-force torques:
%                Q_gas(1) = Q_theta_p^gas
%                Q_gas(2) = Q_theta_c^gas
%                Q_gas(3) = 0
%
%   Sign convention:
%     - x (piston coordinate) is positive upward.
%     - Gas pressure acts downward, so
%         F_gas = -p_cyl * A_piston.
%     - Generalized forces follow from
%         Q_theta_p = F_gas * (∂x/∂theta_p) = F_gas * Ap
%         Q_theta_c = F_gas * (∂x/∂theta_c) = F_gas * Ac.

    % ---------------------------------------------------------------------
    % Piston area
    % ---------------------------------------------------------------------
    if isfield(params, 'A_piston')
        A_piston = params.A_piston;
    elseif isfield(params, 'B')
        % Assume 'B' is bore diameter
        A_piston = 0.25 * pi * params.B.^2;
    else
        error(['generalized_gas_force: piston area not specified. ', ...
               'Provide params.A_piston or params.B (bore diameter).']);
    end

    % ---------------------------------------------------------------------
    % Get Ac, Ap from piston_kinematics
    % ---------------------------------------------------------------------
    % Velocities and accelerations are irrelevant here, so set to zero.
    [~, ~, ~, Ac, Ap] = piston_kinematics(params, ...
        theta_c, 0, 0, ...
        theta_p, 0, 0);

    % ---------------------------------------------------------------------
    % Gas force along piston axis (x positive upward)
    % ---------------------------------------------------------------------
    % Downward force => negative w.r.t. +x:
    F_gas = -p_cyl .* A_piston;  % [N]

    % ---------------------------------------------------------------------
    % Generalized torques due to gas force
    % ---------------------------------------------------------------------
    Q_theta_p = F_gas .* Ap;
    Q_theta_c = F_gas .* Ac;
    Q_theta_o = 0;

    % Assemble vector
    Q_gas = [Q_theta_p; Q_theta_c; Q_theta_o];
end
