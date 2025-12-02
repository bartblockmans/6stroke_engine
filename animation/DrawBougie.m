function [] = DrawBougie(params, spark, flame, Xc, Yc)
%DRAW_BOUJIE  Draw the bougie (spark plug) for the 6-stroke engine.
%
%   DrawBougie(params, spark, flame, Xc, Yc)
%
%   Inputs:
%       params : structure with geometry parameters (from parameters.m)
%       spark : spark strength (0-1)
%       flame : flame level (0-1)
%       Xc, Yc : cylinder cavity coordinates [mm]
% -------------------------------------------------------------------------
% Bart Blockmans, 2024 - bart@blockmans.net
% -------------------------------------------------------------------------

% Check input arguments / defaults
if nargin < 2 || isempty(spark); spark = 0; end
if nargin < 3 || isempty(flame); flame = 0; end
if nargin < 4; Xc = []; end
if nargin < 5; Yc = []; end

% Get parameters
p = [0, params.yb];
B = params.bb;
colors = params.colors;
LW = params.LW;

% Temporary flame magnification factor
flame_mag = 7.5;

% Draw flame
% -------------------------------------------------------------------------
spark_strength = max(min(spark, 1), 0);
flame_level    = max(min(flame, 1), 0);
if spark_strength > 0

    % Flame position (anchor at plug tip)
    x_flame = p(1);
    y_flame = p(2) - 20;

    % Base footprint (L0) and scaled growth
    base_half_w = 100;
    base_half_h = 150;
    size_scale  = 1 + (flame_mag - 1) * flame_level;
    half_w      = base_half_w * size_scale;
    half_h      = base_half_h * size_scale;

    % Generate flame texture
    Nx = 500; Ny = 1000;
    [img_flame, alpha_flame] = generate_blurred_flame(Nx, Ny, 0.15, 0.30, 0.3);

    % Scale color intensity and transparency by spark strength
    color_scale = min(0.3 + 0.7 * spark_strength, 1);
    img_scaled   = img_flame * color_scale;
    alpha_scaled = alpha_flame * spark_strength;

    % Clip flame to cylinder cavity boundary (if provided)
    if ~isempty(Xc) && ~isempty(Yc)
        x_vec = linspace(x_flame - half_w, x_flame + half_w, Nx);
        y_vec = linspace(y_flame - half_h, y_flame + half_h, Ny);
        [Xgrid, Ygrid] = meshgrid(x_vec, y_vec);
        inside = inpolygon(Xgrid, Ygrid, Xc(:), Yc(:));
        
        % Crop image to only include visible region (prevents axis limit issues)
        [row_idx, col_idx] = find(inside);
        if ~isempty(row_idx) && ~isempty(col_idx)
            row_min = max(1, min(row_idx) - 1);
            row_max = min(Ny, max(row_idx) + 1);
            col_min = max(1, min(col_idx) - 1);
            col_max = min(Nx, max(col_idx) + 1);
            
            % Crop arrays
            img_cropped = img_scaled(row_min:row_max, col_min:col_max, :);
            alpha_cropped = alpha_scaled(row_min:row_max, col_min:col_max);
            inside_cropped = inside(row_min:row_max, col_min:col_max);
            
            % Set outside pixels to NaN so they're not plotted
            alpha_cropped(~inside_cropped) = NaN;
            
            % Adjust spatial coordinates to match cropped region
            x_crop_min = x_vec(col_min);
            x_crop_max = x_vec(col_max);
            y_crop_min = y_vec(row_min);
            y_crop_max = y_vec(row_max);
            
            % Use cropped data for plotting
            img_plot = img_cropped;
            alpha_plot = alpha_cropped;
            x_data = [x_crop_min, x_crop_max];
            y_data = [y_crop_min, y_crop_max];
        else
            % No visible pixels, skip plotting
            img_plot = [];
            alpha_plot = [];
            x_data = [];
            y_data = [];
        end
    else
        % No clipping needed, use full image
        img_plot = img_scaled;
        alpha_plot = alpha_scaled;
        x_data = x_flame + half_w * [-1, 1];
        y_data = y_flame + half_h * [-1, 1];
    end

    % Plot flame with dynamic footprint (only if there's visible data)
    if ~isempty(img_plot)
        image('XData', x_data, ...
              'YData', y_data, ...
              'CData', img_plot, ...
              'AlphaData', alpha_plot);
    end

    % Draw cylinder cavity
    plot(Xc, Yc, 'k', 'LineWidth', params.LW);

end

% -------------------------------------------------------------------------

% Height main rectangle
H = 0.75 * B;

% Main rectangle
x_mr = p(1) + [-B/2 -B/2 B/2 B/2 -B/2];
y_mr = p(2) + [-H/2  H/2 H/2 -H/2 -H/2];
fill(x_mr, y_mr, colors.boog,'EdgeColor',[0 0 0], 'LineWidth',LW);

% Thread
x_sr = p(1) + [-B/3 -B/3 B/3 B/3 -B/3];
y_sr = p(2) + [-H/2  H/2 H/2 -H/2 -H/2] - (H/2 + H/2);
fill(x_sr, y_sr, colors.boog,'EdgeColor',[0 0 0], 'LineWidth',LW);

% Thread Y values
y_t = linspace(min(y_sr), max(y_sr),10);
y_t = fliplr(y_t);

% Loop over all threads
for i = 1 : 7

    y1 = y_t(i+1);
    y2 = y_t(i+2);
    yc = (y1 + y2)/2;
    xc = B/2 - B/15;
    x_ti = [-xc -B/3 B/3 xc B/3 -B/3 -xc];
    y_ti = [yc y1 y1 yc y2 y2 yc];
    fill(x_ti, y_ti, colors.boog, 'EdgeColor',[0 0 0],'LineWidth',LW-0.5);

end

% Igniter
x_ig = p(1) + [-B/20 -B/20 B/20 B/20 -B/20];
y_ig = p(2) + [-B/40  B/40 B/40 -B/40 -B/40] - (H/2 + 2*H/2 + B/40);
fill(x_ig, y_ig, colors.boog,'EdgeColor',[0 0 0], 'LineWidth',LW-0.5);
x_ig = p(1) + [-B/3, -B/3, B/20, B/20, -B/3+B/20, -B/3+B/20, -B/3];
y_ig = p(2) + [0, -B/9, -B/9, -B/9+B/40, -B/9+B/40, 0, 0] - (H/2 + 2*H/2);
fill(x_ig, y_ig, colors.boog,'EdgeColor',[0 0 0], 'LineWidth',LW-0.5);

% Bolt head
x_b = p(1) + [-1.3 * B/2, -1.3 * B/2, -B/2, B/2, 1.3 * B/2, 1.3 * B/2, -1.3 * B/2];
y_b = p(2) + [0, H/4, H/4+H/8, H/4+H/8, H/4, 0, 0] + H/2;
fill(x_b, y_b, colors.boog,'EdgeColor',[0 0 0], 'LineWidth',LW);

% Circle
ang = linspace(0,2*pi,100);
x_circle = p(1) + (B/6) * cos(ang);
y_circle = p(2) + (B/6) * sin(ang) + 2.5 * H;
fill(x_circle, y_circle, colors.boog,'EdgeColor',[0 0 0], 'LineWidth',LW);

% Top
x_top = p(1) + [-0.75*B/2, -B/4, B/4, 0.75*B/2, -0.75*B/2];
y_top = p(2) + [0, 1.5*H, 1.5*H, 0, 0] + H/2 + H/4 + H/8;
fill(x_top, y_top, colors.boog,'EdgeColor',[0 0 0], 'LineWidth',LW);


