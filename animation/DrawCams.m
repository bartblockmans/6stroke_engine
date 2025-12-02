function [] = DrawCams(theta_crank, params)

% Cam shaft angle
theta_cam = (theta_crank * pi/180) / 3;

% Add intake & exhaust angles
theta_cam_IV = theta_cam + (params.theta_IV - pi/2);
theta_cam_EV = theta_cam + (params.theta_EV - pi/2);

% Rotation matrices
R_IV = [cos(theta_cam_IV) -sin(theta_cam_IV); sin(theta_cam_IV) cos(theta_cam_IV)];
R_EV = [cos(theta_cam_EV) -sin(theta_cam_EV); sin(theta_cam_EV) cos(theta_cam_EV)];

% Rotate cam profiles
XY_IV = [params.IV_cam_X', params.IV_cam_Y'] * R_IV'; X_IV = XY_IV(:,1); Y_IV = XY_IV(:,2);
XY_EV = [params.EV_cam_X', params.EV_cam_Y'] * R_EV'; X_EV = XY_EV(:,1); Y_EV = XY_EV(:,2);

% Intake cam coordinates
X_IV = params.x_cam_IV + X_IV;
Y_IV = params.y_cam_IV + Y_IV;

% Exhaust cam coordinates
X_EV = params.x_cam_EV + X_EV;
Y_EV = params.y_cam_EV + Y_EV;

% Draw intake cam
fill(X_IV, Y_IV, params.colors.cam,'EdgeColor','none');
plot(X_IV, Y_IV, 'k', 'LineWidth',params.LW);
scatter(params.x_cam_IV, params.y_cam_IV, 30, [0.2 0.2 0.2], 'filled');

% Draw exhaust cam
fill(X_EV, Y_EV, params.colors.cam,'EdgeColor','none');
plot(X_EV, Y_EV, 'k', 'LineWidth',params.LW);
scatter(params.x_cam_EV, params.y_cam_EV, 30, [0.2 0.2 0.2], 'filled');