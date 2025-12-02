function h = draw_rotated_rectangle(center, width, height, angle_deg, edgeColor, faceColor, lineWidth, ndiv_width, ndiv_height, axis_xzy)
% DRAW_ROTATED_RECTANGLE Draws a rotated rectangle with optional internal lines.
%
%   h = draw_rotated_rectangle(center, width, height, angle_deg, edgeColor, faceColor, lineWidth, ndiv_width, ndiv_height)
%
% Inputs:
%   center        - [xc, yc] rectangle center
%   width         - width of the rectangle (horizontal pre-rotation)
%   height        - height of the rectangle (vertical pre-rotation)
%   angle_deg     - rotation in degrees (CCW from x-axis)
%   edgeColor     - border color
%   faceColor     - fill color
%   lineWidth     - border line thickness
%   ndiv_width    - (optional) # of internal lines in width direction (default = 0)
%   ndiv_height   - (optional) # of internal lines in height direction (default = 0)
%
% Output:
%   h             - handle to outer patch

    arguments
        center (1,2) double
        width (1,1) double {mustBePositive}
        height (1,1) double {mustBePositive}
        angle_deg (1,1) double
        edgeColor
        faceColor
        lineWidth (1,1) double {mustBePositive}
        ndiv_width (1,1) double {mustBeNonnegative} = 0
        ndiv_height (1,1) double {mustBeNonnegative} = 0
        axis_xzy (1,1) double = 0
    end


    % Half-dimensions
    hw = width / 2;
    hh = height / 2;

    % Rectangle corners before rotation (CCW)
    corners = [-hw, -hh;
                hw, -hh;
                hw,  hh;
               -hw,  hh]';

    % Rotation matrix
    theta = deg2rad(angle_deg);
    R = [cos(theta), -sin(theta); sin(theta), cos(theta)];

    % Rotate and translate corners
    rotatedCorners = R * corners;
    rotatedCorners(1, :) = rotatedCorners(1, :) + center(1);
    rotatedCorners(2, :) = rotatedCorners(2, :) + center(2);

    % Plot rectangle patch
    h = patch('XData',rotatedCorners(1, :), 'YData', rotatedCorners(2, :), ...
          'FaceColor', faceColor, ...
          'EdgeColor', edgeColor, ...
          'LineWidth', lineWidth);

    % h = patch(rotatedCorners(1, :), rotatedCorners(2, :), 'FaceColor', faceColor, ...
    %           'EdgeColor', edgeColor, 'LineWidth', lineWidth);

    hold_state = ishold;
    hold on;

    % Draw internal width-direction lines (vertical in local frame)
    if ndiv_width > 0
        for i = 1:ndiv_width
            alpha = (i / (ndiv_width + 1)) * 2 - 1;  % from -1 to 1
            x_offset = alpha * hw;

            % Line from bottom to top at fixed x offset
            p1 = [x_offset; -hh];
            p2 = [x_offset;  hh];

            % Rotate and translate
            p = R * [p1, p2];
            p(1,:) = p(1,:) + center(1);
            p(2,:) = p(2,:) + center(2);

            if axis_xzy
                plot3(p(1,:), 0*p(1,:), p(2,:), 'Color', edgeColor, 'LineWidth', 1);
            else
                plot(p(1,:), p(2,:), 'Color', edgeColor, 'LineWidth', 1);
            end
        end
    end

    % Draw internal height-direction lines (horizontal in local frame)
    if ndiv_height > 0
        for j = 1:ndiv_height
            beta = (j / (ndiv_height + 1)) * 2 - 1;  % from -1 to 1
            y_offset = beta * hh;

            % Line from left to right at fixed y offset
            p1 = [-hw; y_offset];
            p2 = [ hw; y_offset];

            % Rotate and translate
            p = R * [p1, p2];
            p(1,:) = p(1,:) + center(1);
            p(2,:) = p(2,:) + center(2);

            if axis_xzy
                plot3(p(1,:), 0*p(1,:), p(2,:), 'Color', edgeColor, 'LineWidth', 1);
            else
                plot(p(1,:), p(2,:), 'Color', edgeColor, 'LineWidth', 1);
            end
        end
    end

    if ~hold_state
        hold off;
    end
end
