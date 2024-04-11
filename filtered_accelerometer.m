%FILTERED_ACCELEROMETER Computes IMU accelerometer signals without gravity.
%   Computes the accelerometer signals for IMU 0 and IMU 1 in the IMU 0
%   coordinate frame.
function [d2r_E0_E0,d2r_E1_E0] = filtered_accelerometer(param, x_arr,t_arr,flag_arr, F_e, k_t)
% © Nils C. A. Wilhelmsen
% 11/04/2024
%% Extract parameters
r_So = param.r_So;
r_w = param.r_w;
mu_k = param.mu_k;
m = param.m;
rho = param.rho;
V = param.V;
I_S = param.I_S;
r_Si = param.r_Si;
g = param.g;
%% Extract states
theta = x_arr(:,1);
d_theta = x_arr(:,2);
phi_S = x_arr(:,3);
d_phi_S = x_arr(:,4);
phi_ds = x_arr(:,5);
d_phi_ds = x_arr(:,6);
%% Reconstruct angular accelerations
% Relative slippage
alpha = -r_So/(r_w - r_So);
Delta_dot_theta = alpha*d_phi_S - d_theta;

% F_f under condition of slipping
F_f = mu_k*abs((m-rho*V)*g*cos(theta) + m*(r_w - r_So)*d_theta.^2).*sign(Delta_dot_theta);

% Accelerations under conditions of slipping
d2_theta_slip = (1/(m*(r_w - r_So)))*(-(m-rho*V)*g*sin(theta) + F_f + F_e);
d2_phi_S_slip = (1/I_S)*(r_Si*k_t.*(d_phi_ds - d_phi_S) + r_So*F_f);

% Accelerations under conditions of not slipping
d2_theta_noslip = (1/((m + (I_S/(r_So^2)))*(r_w - r_So)))*(-(m-rho*V)*g*sin(theta) + F_e - (r_Si/r_So)*k_t.*(d_phi_ds - d_phi_S));
d2_phi_S_noslip = (1/(I_S + m*r_So^2))*(r_Si*k_t.*(d_phi_ds - d_phi_S) + r_So*(m-rho*V)*g*sin(theta) - r_So*F_e);

% Combine vectors with flag array
d2_theta = flag_arr.*d2_theta_slip + (~flag_arr).*d2_theta_noslip;
d2_phi_S = flag_arr.*d2_phi_S_slip + (~flag_arr).*d2_phi_S_noslip;

%% Build accelerometer signals
% Acceleration of sub centre of mass in wellbore-fixed coordinates
d2r_B_A = [((r_w - r_So)*(-d2_theta.*sin(theta) - (d_theta.^2).*cos(theta)))';
              ((r_w - r_So)*(d2_theta.*cos(theta) - (d_theta.^2).*sin(theta)))';];

% Acceleration of rotation matrix between sub-fixed and wellbore-fixed
% coordinates
d2R_B_A = zeros(2,2,length(t_arr));
d2R_B_A(1,1,:) = -d2_theta.*sin(theta)-(d_theta.^2).*cos(theta);
d2R_B_A(1,2,:) = -d2_theta.*cos(theta) + (d_theta.^2).*sin(theta);
d2R_B_A(2,1,:) = d2_theta.*cos(theta) - (d_theta.^2).*sin(theta);
d2R_B_A(2,2,:) = -d2_theta.*sin(theta) - (d_theta.^2).*cos(theta);

% Location of IMU 0 wrt centre of sub
r_E0_B = [r_So*cos(phi_S - theta)'; r_So*sin(phi_S-theta)';];

% Velocity of rotation matrix between sub-fixed and wellbore-fixed
% coordinates
dR_B_A = zeros(2,2,length(t_arr));
dR_B_A(1,1,:) = -d_theta.*sin(theta);
dR_B_A(1,2,:) = -d_theta.*cos(theta);
dR_B_A(2,1,:) = d_theta.*cos(theta);
dR_B_A(2,2,:) = -d_theta.*sin(theta);

% Velocity of IMU 0 wrt centre of sub
dr_E0_B = [(-r_So*(d_phi_S - d_theta).*sin(phi_S - theta))';
            (r_So*(d_phi_S-d_theta).*cos(phi_S - theta))';];

% Rotation matrix between sub-fixed and wellbore-fixed coordinates
R_B_A = zeros(2,2,length(t_arr));
R_B_A(1,1,:) = cos(theta);
R_B_A(1,2,:) = -sin(theta);
R_B_A(2,1,:) = sin(theta);
R_B_A(2,2,:) = cos(theta);

% Acceleration of IMU wrt centre of sub
d2r_E0_B = [-r_So*((d2_phi_S - d2_theta).*sin(phi_S-theta) + ((d_phi_S - d_theta).^2).*cos(phi_S-theta))'; 
            r_So*((d2_phi_S - d2_theta).*cos(phi_S-theta) - ((d_phi_S - d_theta).^2).*sin(phi_S-theta))';];

% Rotation matrix between wellbore-fixed and sub-fixed coordinates
R_A_B = zeros(2,2,length(t_arr));
R_A_B(1,1,:) = cos(theta);
R_A_B(1,2,:) = sin(theta);
R_A_B(2,1,:) = -sin(theta);
R_A_B(2,2,:) = cos(theta);

% Rotation matrix between sub-fixed and IMU-fixed coordinates
R_B_E0 = zeros(2,2,length(t_arr));
R_B_E0(1,1,:) = cos(phi_S-theta);
R_B_E0(1,2,:) = sin(phi_S - theta);
R_B_E0(2,1,:) = -sin(phi_S - theta);
R_B_E0(2,2,:) = cos(phi_S - theta);

% Calculate accelerations of IMUs
d2r_E0_A = zeros(2,length(t_arr));
d2r_E1_A = zeros(2,length(t_arr));

d2r_E0_E0 = zeros(2,length(t_arr));
d2r_E1_E0 = zeros(2,length(t_arr));

for idx=1:length(t_arr)
    % Wellbore-fixed coordinates
    d2r_E0_A(:,idx) = d2r_B_A(:,idx) + d2R_B_A(:,:,idx)*r_E0_B(:,idx) + 2*dR_B_A(:,:,idx)*dr_E0_B(:,idx) + R_B_A(:,:,idx)*d2r_E0_B(:,idx);
    d2r_E1_A(:,idx) = d2r_B_A(:,idx) - d2R_B_A(:,:,idx)*r_E0_B(:,idx) - 2*dR_B_A(:,:,idx)*dr_E0_B(:,idx) - R_B_A(:,:,idx)*d2r_E0_B(:,idx);

    % Convert to IMU0-fixed coordinates
    d2r_E0_E0(:,idx) = R_B_E0(:,:,idx)*R_A_B(:,:,idx)*d2r_E0_A(:,idx);
    d2r_E1_E0(:,idx) = R_B_E0(:,:,idx)*R_A_B(:,:,idx)*d2r_E1_A(:,idx);
end
end