function [Xp, Yp, Up, Vp, Cp, params] = particle_simulation(Xp, Yp, Up, Vp, Cp, Xc, Yc, Ypiston, Vpiston, dflame, params)
%PARTICLE_SIMULATION  Simulate the particle motion in the 6-stroke engine.
%
%   [Xp, Yp, Up, Vp, Cp, params] = particle_simulation(Xp, Yp, Up, Vp, Cp, Xc, Yc, Ypiston, Vpiston, dflame, params)
%
%   Inputs:
%       Xp, Yp : particle positions [mm]
%       Up, Vp : particle velocities [mm/s]
%       Cp : particle colors [0,1]
%       params : structure with geometry parameters (from parameters.m)
%
%   Outputs:
%       Xp, Yp : updated particle positions [mm]
%       Up, Vp : updated particle velocities [mm/s]
%       Cp : updated particle colors [0,1]
%       params : modified structure with updated parameters
% -------------------------------------------------------------------------
% Bart Blockmans, 2024 - bart@blockmans.net
% -------------------------------------------------------------------------

% Get valve states
IV = params.memory.IV_state;
EV = params.memory.EV_state;
SP = params.memory.SP_state;

% Get parameters
Rp = params.particles.Rp;
dt = params.dt_anim;

% Add scavenge particle flow when necessary
% -------------------------------------------------------------------------
if SP && (params.memory.scavenge == 0)

    % Set scavenge
    params.memory.scavenge = 1;

    % Add scavenge particles
    Xp = [Xp; params.particles.Xps];
    Yp = [Yp; params.particles.Yps];
    Up = [Up; 0*params.particles.Xps];
    Vp = [Vp; 0*params.particles.Yps];
    Cp = [Cp; 0*params.particles.Xps];

end

% Split particles when valve state has changed
% -------------------------------------------------------------------------
if params.valve_change

    % Split particles between "in" and "out" particles
    [~, ind_out] = SplitParticles(Xp, Yp, Xc, Yc, params.particles.Rp);

    % Eliminate particles outside cylinder cavity
    Xp(ind_out) = [];
    Yp(ind_out) = [];
    Up(ind_out) = [];
    Vp(ind_out) = [];
    Cp(ind_out) = [];

end

% Add exhaust suction velocity when exhaust is open
% -------------------------------------------------------------------------
if EV

    % Add exhaust suction
    [Up, Vp] = ExhaustSuction(Xp, Yp, Up, Vp, Xc, Yc, params.xep, params.yep, 6e5); 

    % Eliminate particles that are outside the cylinder
    ind_out = (Xp > (0.5*params.B + params.tw)) & (Yp > params.TDC1_y);
    Xp(ind_out) = [];
    Yp(ind_out) = [];
    Up(ind_out) = [];
    Vp(ind_out) = [];
    Cp(ind_out) = [];

end

% Set manifold velocities
% -------------------------------------------------------------------------

% Set manifold velocities
[Up, Vp] = SetManifoldVelocities(Xp, Yp, Up, Vp, IV, EV, SP, Ypiston, params);

% Increase particle kinetic energy during combustion
% -------------------------------------------------------------------------
if dflame ~= 0

    Up = Up * (1 + dflame);
    Vp = Vp * (1 + dflame);

end

% Physics-based update of particle positions & velocities
% -------------------------------------------------------------------------

% Identify piston edges: last 10 segments of the closed polygon
Npoly = numel(Xc);
pistonEdgeIdx = (Npoly - 10) : (Npoly - 1);

% Update particle positions & velocities (moving-wall piston)
if ~isempty(Xp)
    [Xp, Yp, Up, Vp] = UpdateParticles(Xp, Yp, Up, Vp, Xc, Yc, Rp, dt, ...
        pistonEdgeIdx, Vpiston, params.particles.contactP);
end

