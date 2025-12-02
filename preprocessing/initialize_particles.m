function params = initialize_particles(params)
%INITIALIZE_PARTICLES  Initialize particles for 6-stroke engine.
%
%   params = initialize_particles(params)
%
%   Inputs:
%       params      : structure with geometry parameters (from parameters.m)
%
%   Outputs:
%       params      : modified structure with added fields:
% -------------------------------------------------------------------------
% Bart Blockmans, 2024 - bart@blockmans.net
% -------------------------------------------------------------------------

% Get parameters
np = params.particles.np;
np_scav = params.particles.np_scav;
phi_A = params.particles.phi_A;

% % Velocity & time
% % -------------------------------------------------------------------------
% 
% % Crank angular velocity
% crank_rpm = params.N_rpm;
% crank_omega = crank_rpm * pi / 30;
% 
% % Time required for 3 full revolutions of the crank shaft
% t = 6 * pi / crank_omega; params.t_anim = t;
% 
% % Time step
% dt = t / params.nf_anim; params.dt_anim = dt; 

% Compute particle radius
% -------------------------------------------------------------------------

% Cavity boundary & area at theta_crank = 0
[~, ~, ~, ~, ~, ~, coords] = geometry_6stroke(0, params);
[Xc, Yc] = CylinderCavity([coords.X_s, coords.Y_s], params);
A0  = polyarea(Xc, Yc);

% Particle radius
params.particles.Rp = sqrt(max(phi_A * A0 / (np * pi), eps));
Rp = params.particles.Rp;

% Initialize particles in intake manifold
% -------------------------------------------------------------------------

% Compute initial positions and velocities
[Xp, Yp] = InitializeParticles_open(params.xir, params.yir, np, Rp);

% Store
params.particles.Xp = Xp;
params.particles.Yp = Yp;

% Initialize particles in scafold reservoir
% -------------------------------------------------------------------------

% Initialize scavenging particles left port
[Xpsl, Ypsl] = InitializeParticles_open(params.xsrl, params.ysrl, np_scav, Rp);

% Initialize scavenging particles right port
[Xpsr, Ypsr] = InitializeParticles_open(params.xsrr, params.ysrr, np_scav, Rp);

% Combine
Xps = [Xpsl; Xpsr];
Yps = [Ypsl; Ypsr];

% Store
params.particles.Xps = Xps;
params.particles.Yps = Yps;