% Main_engine.m
% -------------------------------------------------------------------------
% 
% This script performs a simplified steady-state simulation of a 6-stroke 
% combustion engine (as patented by Porsche in US2024/0301817 A1)
% 
% -------------------------------------------------------------------------
% Bart Blockmans, 2024 - bart@blockmans.net
% -------------------------------------------------------------------------

%% ========================================================================
%                         ADD PATHS TO FUNCTIONS
% =========================================================================

% Get the directory where Main_engine.m is located
main_dir = fileparts(mfilename('fullpath'));

% Add preprocessing, simulation, and animation folders to path
addpath(fullfile(main_dir, 'preprocessing'));
addpath(fullfile(main_dir, 'simulation'));
addpath(fullfile(main_dir, 'animation'));

%% ========================================================================
%                            LOAD PARAMETERS
% =========================================================================

% Display progress
display_progress(1);

% Load engine parameters
parameters;

% Compute derived parameters
params = derived_parameters(params);

% Load gear parameters - not used in this version

%% ========================================================================
%                       GEOMETRIC ANALYSIS ENGINE
% =========================================================================

% Display progress
display_progress(2);

% Compute engine extrema, second argument is to generate figures
params = extrema_6stroke(params);

% Manifold geometry
params = manifold_geometry(params);

% Plot engine geometry
params = plot_engine(0, params, []); pause(0.1);

% Initialize particles (for animation only)
params = initialize_particles(params);

%% ========================================================================
%                     SOLVE FOR PERIODIC STEADY STATE
% ========================================================================= 

% Display progress
display_progress(3);

% Integrate EOMs & find periodic steady state
results = engine_cycle_steady_6stroke(params);

% Display progress
display_progress(4);

% Simulate full 6-cylinder engine
results = full_engine_simulation_6stroke(results, params);

% Display progress
display_progress(5);

% Plot results
plot_results_6stroke(results);

%% ========================================================================
%                           ANIMATE RESULTS
% ========================================================================= 

% Display progress
display_progress(6);

% Animate results
animate_results_6stroke(params, results);

% -------------------------------------------------------------------------