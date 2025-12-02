function [] = DrawPiston(p, params)
%DRAW_PISTON  Draw the piston for the 6-stroke engine.
%
%   DrawPiston(p, params)
%
%   Inputs:
%       p : piston position [mm]
%       params : structure with geometry parameters (from parameters.m)
%
%   Outputs:
% -------------------------------------------------------------------------
% Bart Blockmans, 2024 - bart@blockmans.net
% -------------------------------------------------------------------------

% Get parameters
W = params.B;
H = params.H;
D = 2 * params.Rpis;
C = 0.025 * W;
colors = params.colors;
LW = params.LW;

% Width piston
Wp = W - C; 

% Piston center
pp = [p(1); p(2) + H/6];

% Draw main rectangle
draw_rotated_rectangle(pp, Wp, H, 0, [0 0 0], colors.pis, LW);

% Coloring
draw_rotated_rectangle([p(1) + W/4; pp(2)], Wp/5, H, 0, 'none',1.5*colors.pis, 1);
draw_rotated_rectangle(pp, Wp, H, 0, [0 0 0], 'none', LW);

% Draw circle
DrawSimpleCircle(p, D/2, colors.rod, LW);

% Height rings
Hr = (H/3) / 9;

% Draw 3 piston rings
for i = 1 : 3

    % Center
    pr = p;
    pr(2) = pr(2) + H/3 + Hr + (i-1) * (3*Hr);

    % Left & Right coordinates
    prl = pr;
    prl(1) = prl(1) - Wp/2;
    prr = pr;
    prr(1) = prr(1) + Wp/2;

    % Plot piston ring
    DrawPulley(prl, prr, C/2, C/2, [0.2 0.2 0.2], 1, 10);

end