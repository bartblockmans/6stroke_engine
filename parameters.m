% parameters for the Main_engine.m script
% 
% Angle convention for 6-stroke:
%   phi = 0 deg  : TDC1 (start of Stroke A)
%   Strokes (each 180°):
%     A:   0–180   intake (head intake valve)           [long stroke]
%     B: 180–360   compression #1                       [long stroke]
%     C: 360–540   combustion/expansion #1 (+ late EVO) [long stroke]
%     D: 540–720   scavenge+mix + compression #2        [short stroke]
%     E: 720–900   combustion/expansion #2              [short stroke]
%     F: 900–1080  exhaust (head exhaust valve)         [long stroke]
% -------------------------------------------------------------------------
% Bart Blockmans, 2024 - bart@blockmans.net
% -------------------------------------------------------------------------

% =========================================================================
%                   OPERATING CONDITIONS & ENGINE LAYOUT
% =========================================================================

params.N_rpm = 2000;    % Engine speed [rpm]
params.N_cyl = 1;       % Number of cylinders [-]

% =========================================================================
%                       GEOMETRICAL PARAMETERS
% =========================================================================

% Parameters used in simulation
% -------------------------------------------------------------------------
params.B     = 100;         % cylinder bore [mm]
params.Rr    = 100;         % ring gear pitch radius [mm]
params.Rp    = 60;          % planet gear pitch radius [mm]
params.mn    = 4;           % Normal modulus gears [mm]
params.delta = 20;          % offset of connecting rod [mm]
params.L     = 200;         % lenth of connecting rod [mm]
params.Vc    = 8e4;         % clearance volume (at highest TDC) [mm³]
params.psi0  = pi/2;        % body-fixed angle of connecting rod [rad]

% Parameters used for animation
% -------------------------------------------------------------------------
params.H    = 60;           % total piston height [mm]
params.Hc   = 5;            % cylinder clearance height [mm]
params.Rcw  = 80;           % counter weight crank shaft radius [mm]
params.Rc   = 100;          % top of the cylinder radius [mm]
params.Ds   = 30;           % crank shaft diameter [mm]
params.Do   = 224;          % ring gear outside diameter [mm]
params.Di   = 96;           % planet gear inside diameter [mm]
params.to   = 10;           % rim thickness ring gear [mm]
params.tw   = 15;           % wall thickness cylinder [mm]
params.Rcra = 12.5;         % crank shaft radius [mm]
params.Rcr1 = 20;           % connecting rod large radius [mm]
params.Rcr2 = 10;           % connecting rod small radius [mm]
params.Rpis = 10;           % piston hinge radius [mm]
params.yb = 330;            % igniter mount height [mm]
params.bb = 20;             % igniter height [mm]

% =========================================================================
%                        FLAT-SIX ENGINE LAYOUT
% =========================================================================

% Cylinder offset angles (from 180° spacing)
params.cyl_offset = [0, 0, 0, 0, 0, 0]; % [deg, deg, deg, deg, deg, deg]

% =========================================================================
%                   MANIFOLDS & SCAVENGE PORT PARAMETERS
% =========================================================================

params.p_int    = 1.05e5;   % intake plenum [Pa]
params.T_int    = 315;      % intake temp [K]
params.p_exh    = 1.15e5;   % exhaust plenum [Pa]
params.T_exh    = 800;      % exhaust temp [K]
params.p_scav   = 1.20e5;   % scavenge supply [Pa] (can be boosted)
params.T_scav   = 320;      % scavenge temp [K]
params.p_crank  = 1.00e5;   % crankcase (for net piston force)

% =========================================================================
%                 FUEL, CHEMISTRY & SPLIT BETWEEN BURNS
% =========================================================================

params.LHV        = 44e6;       % gasoline-like [J/kg]
params.AFR_st     = 14.7;       % stoichiometric air/fuel (mass) [-]
params.YO2_air    = 0.232;      % O2 mass fraction in air [-]
params.YO2_exh    = 0.00;       % O2 fraction in exhaust plenum [-]
params.eta_comb1  = 0.97;       % combustion efficiency burn #1 [-]
params.eta_comb2  = 0.97;       % combustion efficiency burn #2 [-]
params.mdot_fuel  = 4.5e-4;     % [kg/s] per cylinder 
params.fuel_split = 0.70;       % fraction to burn #1 (rest to burn #2) [-]

% Compute O2 required for stoichiometric combustion
params.nu_O2    = params.AFR_st * params.YO2_air; % ~3.41 kgO2/kg fuel

% =========================================================================
%                   COMBUSTION TIMING (WIEBE PER BURN)
% =========================================================================

% Burn #1 (stroke C: 360–540°)
params.soc1       = 345;        % start of combustion #1 [deg]
params.dur1       = 45;         % duration #1 [deg]
params.m1         = 2.0;        % Wiebe m
params.a1         = 6.9;        % Wiebe a (~99% mass-fraction burned)

% Burn #2 (stroke E: 720–900°)
params.soc2       = 725;        % start of combustion #2 [deg]
params.dur2       = 40;         % duration #2 [deg]
params.m2         = 2.0;        % Wiebe m
params.a2         = 6.9;        % Wiebe a

% =========================================================================
%                   VALVE & PORT TIMING & GEOMETRY
% =========================================================================

% Head intake valve (stroke A)
params.IVO_A   = -20;           % (wraps to 1060) open just before 0
params.IVC_A   = 140;           % closes well before 180
params.D_int   = 0.032;         % inlet diameter [m]
params.Lift_int= 0.008;         % inlet lift [m]
params.Cd_int  = 0.80;          % discharge coefficient for inlet

% Head exhaust valve: late blowdown in C + full exhaust in F
% Split into two windows: C blowdown and F exhaust (kept closed during D,E)
params.EVO_C1  = 440;           % open blowdown late in C
params.EVC_C1  = 580;           % close before scavenge / compression-2
params.EVO_F1  = 900;           % open for exhaust stroke F
params.EVC_F1  = 1080;          % close at end of cycle
params.D_exh   = 0.027;         % exhaust diameter [m]
params.Lift_exh= 0.007;         % exhaust lift [m]
params.Cd_exh  = 0.80;          % discharge coefficient for exhaust

% Scavenge ports: around BDC2 into early D (540–660)
% params.SP_open = 480;         % open at BDC2 (automatically computed)
% params.SP_close= 600;         % close early in D (automatically computed)
params.SP_perim= 0.140;         % total perimeter of port windows [m] 
params.SP_hmax = 0.009;         % max effective port height [m]
params.Cd_scav = 0.85;          % discharge coefficient for ports

% =========================================================================
%                     GAS PROPERTIES & HEAT TRANSFER
% =========================================================================

% Gas properties
params.R_mix   = 287;           % air-like gas constant [J/(kg K)]
params.cp      = 1005;          % constant specific heat capacity [J/(kg K)]

% Woschni heat-transfer coefficients
params.C1_wosch= 3.26;          % correlation constant (SI units, p in bar)
params.C2_wosch= 6.18;          % pseudo-velocity constant
params.T_wall  = 450;           % average wall temperature [K]

% Adiabatic index
params.gamma   = params.cp/(params.cp - params.R_mix); % gamma = cp/cv

% =========================================================================
%                   SCAVENGE ALGEBRAIC MODEL PARAMETERS
% =========================================================================

% Trapping efficiency for fresh scavenge air (fraction of inflow that stays)
params.eta_tr   = 0.90;        % 0.8–0.95 typical with good uniflow

% Residual displacement efficiency: additional mass expelled per unit trapped inflow
params.eta_mix  = 0.35;        % 0.2–0.5 (tunes residual removal without CFD)

% Under-relaxation for burn-2 O2 cap between cycles
params.lambda2_relax = 0.5;    % relaxation factor [-]

% =========================================================================
%                          NUMERICAL PARAMETERS
% =========================================================================

params.maxCycles = 100;         % maximum number of cycles to iterate
params.cycleTol  = 1e-3;        % convergence tolerance on state norm [Pa + kg]
params.phiSpan   = [0 1080];    % crank-angle span [deg]
params.phiGrid   = linspace(0,1080,2161); % desired output grid [deg] (0.5 deg)

% =========================================================================
%                          FULL ENGINE PARAMETERS
% =========================================================================

% Mass and inertia
params.mp           = 0.65;     % Piston mass, including conrod [kg]
params.mpg          = 1.2;      % Planetary gear mass [kg]
params.Jp           = 1.5e-3;   % Planetary gear inertia [kg·m^2]
params.Jc           = 0.15;     % Carrier inertia [kg·m^2]
params.Jo           = 0.045;    % Output side inertia [kg·m^2]

% Gear mesh
params.k_mesh_avg   = 2.5e7;    % Average mesh stiffness along LOA [N/m]
params.k_mesh_amp   = 0.25e7;   % 1st-order amplitude mesh stiffness [N/m]
params.c_mesh       = 2.0e2;    % Mesh damping along LOA [Ns/m]

% Dual-mass flywheel (DMF)
params.k_dmf        = 3.0e4;    % DMF torsional stiffness [Nm/rad]
params.c_dmf        = 10;       % DMF torsional damping [Nms/rad]
params.phi_dmf0     = 0;        % DMF preload [rad]

% Pre-DMF losses (engine friction + accessories)
params.loss.T_fric_const = 10;   % constant friction [N·m]
params.loss.T_fric_visc  = 0.02; % viscous friction [N·m·s]
params.loss.T_acc_const  = 2;    % accessory drag [N·m]
params.loss.T_acc_visc   = 0.01; % accessory viscous drag [N·m·s]

% DMF simulation controls
params.sim.n_cycles_dmf  = 4;     % number of 1080° windows to settle DMF
params.sim.ode_rel_tol   = 1e-5;  % relative tolerance for ODE solver
params.sim.ode_abs_tol   = 1e-5;  % absolute tolerance for ODE solver
params.sim.ode_max_step  = 1e-3;  % maximum step size for ODE solver [s]
params.sim.ode_show_progress = true;

% =========================================================================
%                          ANIMATION PARAMETERS
% =========================================================================

% Number of frames
params.nf_anim = 100; 

% Line width
params.LW = 1.5; 

% Colors
params.colors.cyl = [186 197 204] / 255;    % Cylinder color
params.colors.rod = [204 203 172] / 255;    % Connecting rod color
params.colors.cw = [159 163 130] / 255;     % Counterweight color
params.colors.gear = [112 146 190] / 255;   % Planet gear body color
params.colors.pis = [119 114 113] / 255;    % Piston color
params.colors.cham = [202 208 234] / 255;   % Chamber underneath piston color
params.colors.valve = [215 245 246] / 255;  % Valve color
params.colors.boog = [149 142 141] / 255;   % Boogie color
params.colors.cam = [186 197 204] / 255;    % Cam color
params.colors.man = [202 208 234] / 255;    % Empty manifold color

% Particle simulation
params.particles.anim = 0;          % Animate particles? 1 = yes
params.particles.np = 200;          % # particles regular intake
params.particles.np_scav = 120;     % # particles per scavenge port
params.particles.phi_A = 0.2;       % target area fraction of particles in cavity
params.particles.speedRms = 10;     % RMS of initial velocity distribution [m/s]
params.particles.contactP = 1;      % Contact between particles? 1 = yes