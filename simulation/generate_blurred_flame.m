function [img, alpha] = generate_blurred_flame(Nx, Ny, radius_x_frac, radius_y_frac, blur_sigma)
    % Nx, Ny: image size
    % radius_x_frac, radius_y_frac: relative size of the "core"
    % blur_sigma: standard deviation of Gaussian blur

    [X, Y] = meshgrid(linspace(-1,1,Nx), linspace(-1,1,Ny));
    R = sqrt((X/radius_x_frac).^2 + (Y/radius_y_frac).^2);

    % Alpha mask with elliptical Gaussian falloff
    alpha = exp(-R.^2 / (2 * blur_sigma^2));
    alpha(R > 1) = 0;

    % Flame color: white/yellow core fading to orange
    img = zeros(Ny, Nx, 3);
    img(:,:,1) = 1.0;      % Red channel
    img(:,:,2) = 0.6;      % Green
    img(:,:,3) = 0.1;      % Blue
end