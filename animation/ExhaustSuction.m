function [Up, Vp] = ExhaustSuction(Xp, Yp, Up, Vp, Xc, Yc, Xe, Ye, strength)
% ExhaustSuction
% Applies a suction velocity to particles based on their location:
% - Particles inside the cylinder: pulled towards the exhaust port (inverse distance)
% - Particles outside the cylinder (in exhaust manifold): pushed to the right
% The computed suction velocity is averaged with the current velocity to provide
% gradual changes and prevent velocity blow-up.
%
% Inputs:
%   Xp, Yp - particle positions (np x 1 arrays)
%   Up, Vp - particle velocities (np x 1 arrays) - will be modified
%   Xc, Yc - cylinder cavity polygon coordinates (closed polygon)
%   Xe, Ye - exhaust port boundary coordinates (2-element arrays: [left, right])
%   strength - (optional) suction strength parameter (default: 1000)
%             Higher values = stronger suction. Tune by trial and error.
%
% Outputs:
%   Up, Vp - updated particle velocities (averaged with computed suction velocity)
% -------------------------------------------------------------------------
% Bart Blockmans, 2024 - bart@blockmans.net
% -------------------------------------------------------------------------

    % Default suction strength (tunable parameter)
    if nargin < 9 || isempty(strength)
        strength = 1000;
    end
    
    % Get number of particles
    np = numel(Xp);
    if np == 0
        return;
    end
    
    % Check which particles are inside the cylinder cavity
    [in, on] = inpolygon(Xp, Yp, Xc, Yc);
    inside_cylinder = in | on;
    outside_cylinder = ~inside_cylinder;
    
    % Initialize suction velocity components
    Up_suction = zeros(np, 1);
    Vp_suction = zeros(np, 1);
    
    % For particles inside the cylinder: apply suction towards exhaust port
    if any(inside_cylinder)
        % Calculate exhaust port center (midpoint between left and right boundaries)
        if numel(Xe) ~= 2 || numel(Ye) ~= 2
            error('ExhaustSuction: Xe and Ye must each have 2 elements (left and right boundary points)');
        end
        
        Xe_center = mean(Xe);
        Ye_center = mean(Ye);
        
        % Small epsilon to avoid division by zero for particles exactly at the port
        eps_dist = 1e-6;
        
        % Calculate distance and direction from each particle to exhaust port center
        dx = Xe_center - Xp(inside_cylinder);
        dy = Ye_center - Yp(inside_cylinder);
        dist = hypot(dx, dy);
        
        % Avoid division by zero
        dist = max(dist, eps_dist);
        
        % Calculate suction velocity magnitude: inversely proportional to distance
        % The strength parameter controls the overall magnitude
        v_suction_mag = strength ./ dist;
        
        % Normalize direction vectors
        dx_norm = dx ./ dist;
        dy_norm = dy ./ dist;
        
        % Calculate suction velocity components for inside particles
        Up_suction(inside_cylinder) = v_suction_mag .* dx_norm;
        Vp_suction(inside_cylinder) = v_suction_mag .* dy_norm;
    end
    
    % For particles outside the cylinder (in exhaust manifold): push to the right
    if any(outside_cylinder)
        % Apply velocity boost to the right (positive x-direction)
        % Use a fraction of the strength parameter for the exhaust manifold push
        exhaust_boost = 0.5 * strength;  % Can be tuned separately if needed
        
        Up_suction(outside_cylinder) = exhaust_boost;
        % No vertical component for exhaust manifold particles
        Vp_suction(outside_cylinder) = 0;
    end
    
    % Average current velocity with computed suction velocity
    % This provides gradual changes and prevents velocity blow-up
    Up = 0.5 * (Up + Up_suction);
    Vp = 0.5 * (Vp + Vp_suction);

end
