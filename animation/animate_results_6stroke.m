function [] = animate_results_6stroke(params, results, gear1, gear2, pair, nf)
%ANIMATE_RESULTS_6STROKE  Animate the results of a 6-stroke engine simulation.
%
%   animate_results_6stroke(params, results, gear1, gear2, pair, nf)
%
%   Inputs:
%       params : structure with geometry parameters (from parameters.m)
%       results : structure with simulation results (from engine_cycle_steady_6stroke.m)
%       gear1 : gear 1 (not used in this version)
%       gear2 : gear 2 (not used in this version)
%       pair : pair of gears (not used in this version)
%       nf : number of frames (default 100)
%
%   Outputs:
%       None
% -------------------------------------------------------------------------
% Bart Blockmans, 2024 - bart@blockmans.net
% -------------------------------------------------------------------------

% Check input arguments
if nargin < 3 || isempty(gear1); gear1 = []; end
if nargin < 4 || isempty(gear2); gear2 = []; end
if nargin < 5 || isempty(pair); pair = []; end
if nargin < 6; nf = 100; end

% Crank angular velocity
crank_rpm = params.N_rpm;
crank_omega = crank_rpm * pi / 30;

% Time required for 3 full revolutions of the crank shaft
t = 6 * pi / crank_omega; params.t_anim = t;

% Time step
dt = t / nf; params.dt_anim = dt; 

% Crank angle step
theta_crank_step = crank_omega * dt * 180/pi; % [deg]

% Crank angle
theta_crank = 0 : theta_crank_step : 1080; % [deg]

% Reverse crank angle for clockwise rotation
theta_crank = -theta_crank;

% Get max temperature and max pressure and store in results_i
results_i.T_max = max(results.T_cyl);
results_i.T_min = min(results.T_cyl);
results_i.P_max = max(results.p_cyl);
results_i.P_min = min(results.p_cyl);

% Max valve & scavenge ports displacements
results_i.IV_max = max(results.L_int_m);
results_i.EV_max = max(results.L_exh_m);
results_i.SP_max = max(results.L_scav_m);

% Initialize animation memory
params.memory.IV_state = 1; 
params.memory.EV_state = 0;
params.memory.SP_state = 0;
params.memory.Xp = params.particles.Xp;
params.memory.Yp = params.particles.Yp;
params.memory.Up = 0*params.particles.Xp;
params.memory.Vp = 0*params.particles.Yp;
params.memory.Cp = 0*params.particles.Xp;
params.memory.scavenge = 0; 
params.memory.flame = 0;

% Set speed scaling factor to 1
params.speed_scaling = 1;

% Design cam
params = cam_design(params, results);

% Create figure with white background
figure('Color', 'w');
hold on; axis equal;

% Loop over all crank angles
for i = 1 : length(theta_crank)

    % Clear figure
    cla; hold on; axis equal;

    % Current crank angle
    theta_crank_i = theta_crank(i);

    % Pressure
    P_i = interp1(results.phi_deg, results.p_cyl, -theta_crank_i, "pchip","extrap");

    % Temperature
    T_i = interp1(results.phi_deg, results.T_cyl, -theta_crank_i, "pchip","extrap");

    % Valve & scavenge ports displacement
    IV_i = interp1(results.phi_deg, results.L_int_m, -theta_crank_i, "pchip","extrap");
    EV_i = interp1(results.phi_deg, results.L_exh_m, -theta_crank_i, "pchip","extrap");
    SP_i = interp1(results.phi_deg, results.L_scav_m, -theta_crank_i, "pchip","extrap");

    % Spark & flame
    spark_i = interp1(results.phi_deg, results.spark, -theta_crank_i, "pchip","extrap");
    flame_i = interp1(results.phi_deg, results.flame, -theta_crank_i, "pchip","extrap");

    % Store in results_i structure
    results_i.P = P_i;
    results_i.T = T_i;
    results_i.IV = IV_i;
    results_i.EV = EV_i;
    results_i.SP = SP_i;
    results_i.spark = spark_i;
    results_i.flame = flame_i;

    % Draw engine
    params = plot_engine(theta_crank_i, params, results_i, gear1, gear2, pair);

    % Draw now
    drawnow();

end