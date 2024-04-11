%% Main code file for "Modeling and Detection of Slipping in Damping Subs for Drilling"
% © Nils C. A. Wilhelmsen
% 11/04/2024
clear all; close all;
%% Model parameters
% Damping sub
m_S = 1e3;                                                                  % Mass of damping sub [kg]
l_S = 0.2;                                                                  % Length of damping sub [m]
r_So = 0.08;                                                                % Outer radius of damping sub [m]
r_Si = 0.06;                                                                % Inner radius of damping sub [m]
I_S = 0.5*m_S*(r_Si + r_So)^2;                                              % Sub moment of inertia [kg*m^2]

% Drillstring
m_ds = 1e4;                                                                 % Mass of drillstring segment [kg]
l_ds = 800;                                                                 % Length of drillstring segment [m]
r_ds = 0.055;                                                               % Outer radius of drillstring segment [m]
I_ds = m_ds*r_ds^2;                                                         % Drillstring moment of inertia [kg*m^2]
I_td = 2900;                                                                % Top-drive moment of inertia [kg*m^2]
c = 1.81e3;                                                                 % Drillstring rotational spring constant [N*m/rad]
n = 6;                                                                      % Number of blades in drag bit [-]
a = 0.159;                                                                  % Radius of drag bit [m]
xi = 1;                                                                     % Parameter characterizing spatial distribution of wearflats [-]
mu = 0.6;                                                                   % Bit friction coefficient [-]
l = 1e-3;                                                                   % Length of wearflats [m]
sigma_bar = 6e7;                                                            % Bit contact stress [Pa]
t_drill = 10;                                                               % Time at which drillbit comes in contact with rocks [s]
k_P = 2000;                                                                 % Top drive controller P gain [N*m*s/rad]
k_I = 1000;                                                                 % Top drive controller I gain [N*m/rad]
d_phi_td_setpoint = 15;                                                     % Top drive velocity setpoint [rad/s]

% Well
r_w = 0.1;                                                                  % Inner radius of well [m]
mu_s = 0.3;                                                                 % Static friction coefficient of wellbore [-]
mu_k = 0.2;                                                                 % Kinetic friction coefficient of wellbore [-]

% Miscallaneous
g = 9.8;                                                                    % Acceleration of gravity [m/s^2]
rho = 1200;                                                                 % Density of drilling mud [kg/m^3]
m = m_S + m_ds;                                                             % Combined sub and drillstring segment mass [kg]
V = pi*((l_ds - l_S)*r_ds^2 + l_S*r_So^2);                                  % Net sub and drillstring segment volume [m^3]
k_t = 200;                                                                  % Magnetic damping coefficient [N*s/rad]
alpha = -r_So/(r_w - r_So);                                                 % Proportionality constant during rolling [-]
%% Initialize simulation
% Define parameters for integration function
param.m = m;
param.V = V;
param.rho = rho;
param.g = g;
param.r_w = r_w;
param.r_So = r_So;
param.r_Si = r_Si;
param.r_ds = r_ds;
param.I_S = I_S;
param.I_ds = I_ds;
param.r_w = r_w;
param.mu_k = mu_k;
param.mu_s = mu_s;
param.alpha = alpha;
param.I_td = I_td;
param.c = c;
param.n = n;
param.a = a;
param.xi = xi;
param.mu = mu;
param.l = l;
param.sigma_bar = sigma_bar;
param.t_drill = t_drill;
param.k_P = k_P;
param.k_I = k_I;
param.d_phi_td_setpoint = 15;

% Define initial conditions
theta0 = 0;                                                                 % Initial sub location in wellbore
d_theta0 = 0;                                                               % Initial sub angular velocity in wellbore
phi_S0 = 0;                                                                 % Initial sub rotational angle
d_phi_S0 = 0;                                                               % Initial sub rotational velocity
phi_ds0 = 0;                                                                % Initial drillstring rotational angle
d_phi_ds0 = 0;                                                              % Initial drillstring rotational velocity
phi_td0 = 0;                                                                % Initial topdrive rotational angle
d_phi_td0 = 0;                                                              % Intiial topdrive rotational velocity
integrator0 = 0;                                                            % Initial value of integrator term in PI controller

x0 = [theta0;d_theta0;phi_S0;d_phi_S0;phi_ds0;d_phi_ds0;%]; 
        phi_td0; d_phi_td0;integrator0;];                                   % Vector containing initial conditions

% Simulation initial and final times
t0 = 0;                                                                     % Simulation intitial time
tf = 20;                                                                    % Simulation final time     

% Define force signal
F_e = @(t) 0*t;                                                             % External force

% Initialize integration stop time to initial time
tend = t0;

% Global storage arrays
x_arr = [];                                                                 % Storage array for states
t_arr = [];                                                                 % Storage array for time
flag_arr = [];                                                              % Storage array for flags indicating whether slipping or not (1=slipping, 0=not slipping)

% Compute initial friction and normal forces
d2_theta0 = (1/((m + I_S/(r_So^2))*(r_w - r_So))*(-(m-rho*V)*g*sin(theta0) + F_e(0) - (r_Si/r_So)*k_t*(d_phi_ds0 - d_phi_S0)));
F_f0 = (-I_S/(r_So^2))*(r_w - r_So)*d2_theta0 - (r_Si/r_So)*k_t*(d_phi_ds0 - d_phi_S0);
F_N0 = (m-rho*V)*g*cos(theta0) + m*(r_w - r_So)*d_theta0^2;

% Initialize slipping flag
if(abs(F_f0) > mu_s*F_N0)
    slipping = true;
else
    slipping = false;
end

%% Simulation via events loop
while(tend < tf)
    % Specify which event to monitor for
    if(slipping)
        options = odeset('Events', @(t,x) rolling_event(t,x,param));
    else
        options = odeset('Events', @(t,x) slipping_event(t,x, param, k_t, F_e(t)));
    end
    
    % Run simulation
    [t,x,te,xe,ie] = ode23(@(t,x) drillstring_sub_coupled_model(t,x,param,slipping, k_t, F_e(t)), [tend tf], x0, options);

    % Calculate recent vector of flags
    flags = slipping*ones(size(t));
    
    % Truncate true flags to interior points
    flags(1)=false;flags(end)=false;

    % Flip slipping flag
    if(slipping)
        % Obtain latest state vector
        x_end = x(end,:);        
        
        % Extract latest states
        theta_end = x_end(1);
        d_theta_end = x_end(2);
        d_phi_S_end = x_end(4);
        d_phi_ds_end = x_end(6);

        % Compute friction and normal forces
        d2_theta_end = (1/((m + I_S/(r_So^2))*(r_w - r_So))*(-(m-rho*V)*g*sin(theta_end) + F_e(t(end)) - (r_Si/r_So)*k_t*(d_phi_ds_end - d_phi_S_end)));
        F_f_end = (-I_S/(r_So^2))*(r_w - r_So)*d2_theta_end - (r_Si/r_So)*k_t*(d_phi_ds_end - d_phi_S_end);
        F_N_end = (m-rho*V)*g*cos(theta_end) + m*(r_w - r_So)*d_theta_end^2;

        % Check whether friction small enough for rolling without slipping
        if(abs(F_f_end) <= mu_s*abs(F_N_end))
            slipping = false;
        end
    else
        slipping = true;
    end

    % Update tend
    if(~isempty(te))
        % Set next simulation start time to last time instant
        tend = t(end);
    else
        tend = tf;
    end

    % Store simulation output in global arrays
    x_arr = [x_arr;x;];
    t_arr = [t_arr;t;];
    flag_arr = [flag_arr;flags;];

    % Set new initial conditions to final conditions
    x0 = x(end,:);
end
% Extract states
theta = x_arr(:,1);
d_theta = x_arr(:,2);
phi_S = x_arr(:,3);
d_phi_S = x_arr(:,4);
phi_ds = x_arr(:,5);
d_phi_ds = x_arr(:,6);
phi_td = x_arr(:,7);
d_phi_td = x_arr(:,8);
%% Slip detection
% Calculate filtered accelerometer signals
[d2r_E0_E0, d2r_E1_E0] = filtered_accelerometer(param,x_arr,t_arr,flag_arr,F_e(t_arr), k_t);

% Add noise to accelerometer signals
[col, row] = size(d2r_E0_E0);
d2r_E0_E0_noise = d2r_E0_E0 + wgn(col, row,1, 1e-3);
d2r_E1_E0_noise = d2r_E1_E0 + wgn(col, row,1, 1e-3);


% Filter accelerometer signals
filt_order = 10;
d2r_E0_E0_filtered = [filter(ones(filt_order, 1)/filt_order, 1, d2r_E0_E0_noise(1,:));
                       filter(ones(filt_order, 1)/filt_order, 1, d2r_E0_E0_noise(2,:))];

d2r_E1_E0_filtered = [filter(ones(filt_order, 1)/filt_order, 1, d2r_E1_E0_noise(1,:));
                       filter(ones(filt_order, 1)/filt_order, 1, d2r_E1_E0_noise(2,:))];

% Calculate sums and differences of accelerometer signals
m0 = 0.5*(d2r_E0_E0_filtered + d2r_E1_E0_filtered);
m1 = 0.5*(d2r_E0_E0_filtered - d2r_E1_E0_filtered);

% Norm containing only theta-based signals
m0_norm = sqrt(sum(m0.^2));

% Reconstruct norm containing only theta-based signals in case of no slip
m1_x = m1(1,:);
m1_y = m1(2,:);

% Compute sigma signal
sigma = (alpha*m1_x).^2 + m1_y.^2;

% Compute residual
residual = sqrt(sigma) - m0_norm;
%% Make plots
% Initialize plot counter
i=1;

% Last time of plot
t_plot = tf;

% Plot topdrive velocity
figure(i); i=i+1;
subplot(2,1,1);
hold on;
plot(t_arr, d_phi_td, 'linewidth', 2);
box on; grid on;
ylim([0 35]);
xlim([0 t_plot]);
ylabel('Angular velocity $[\frac{rad}{s}]$', 'interpreter', 'latex', 'fontsize', 16);
xlabel('Time $[s]$', 'interpreter', 'latex', 'fontsize', 16);
legend('$\dot\phi_{td}(t)$','interpreter', 'latex', 'fontsize', 16, 'location', 'northwest');
%
% Plot drillstring velocity
subplot(2,1,2);
hold on;
plot(t_arr, d_phi_ds, 'linewidth', 2);
box on; grid on;
ylim([0 35]);
xlim([0 t_plot]);
ylabel('Angular velocity $[\frac{rad}{s}]$', 'interpreter', 'latex', 'fontsize', 16);
xlabel('Time $[s]$', 'interpreter', 'latex', 'fontsize', 16);
legend('$\dot\phi_{ds}(t)$', '$\dot\phi_{ds}(t)$ with damping','interpreter', 'latex', 'fontsize', 16, 'location', 'northwest');

% Plot theta and theta_dot
figure(i);i=i+1;
%
subplot(2,1,1);
plot(t_arr, theta, 'linewidth', 2);
box on; grid on;
ylim([-pi pi]);
xlim([0 t_plot]);
ylabel('Angle $[rad]$', 'interpreter', 'latex', 'fontsize', 16);
xlabel('Time $[s]$', 'interpreter', 'latex', 'fontsize', 16);
legend('$\theta(t)$', 'interpreter', 'latex', 'fontsize', 16, 'location', 'northwest');
%
subplot(2,1,2);
plot(t_arr, d_theta, 'linewidth', 2);
box on; grid on;
ylim([-3*pi 3*pi]);
xlim([0 t_plot]);
ylabel('Angular velocity $[\frac{rad}{s}]$', 'interpreter', 'latex', 'fontsize', 16);
xlabel('Time $[s]$', 'interpreter', 'latex', 'fontsize', 16);
legend('$\dot\theta(t)$', 'interpreter', 'latex', 'fontsize', 16, 'location', 'northwest');

% Plot phi_S and phi_S_dot
figure(i);i=i+1;
%
subplot(2,1,1);
plot(t_arr, mod(phi_S+pi, 2*pi)-pi, 'linewidth', 2);
box on; grid on;
ylim([-pi pi]);
xlim([0 t_plot]);
ylabel('Angle $[rad]$', 'interpreter', 'latex', 'fontsize', 16);
xlabel('Time $[s]$', 'interpreter', 'latex', 'fontsize', 16);
legend('$\phi_S(t)$', 'interpreter', 'latex', 'fontsize', 16, 'location', 'northwest');
%
subplot(2,1,2);
plot(t_arr, d_phi_S, 'linewidth', 2);
box on; grid on;
ylim([-pi pi]);
xlim([0 t_plot]);
ylabel('Angular velocity $[\frac{rad}{s}]$', 'interpreter', 'latex', 'fontsize', 16);
xlabel('Time $[s]$', 'interpreter', 'latex', 'fontsize', 16);
legend('$\dot\phi_S(t)$', 'interpreter', 'latex', 'fontsize', 16, 'location', 'northwest');

% Plot relative slippage
figure(i);i=i+1;
plot(t_arr, -alpha*d_phi_S + d_theta, 'linewidth', 2);
box on; grid on;
xlim([0 t_plot]);
ylim([-1 6]);
ylabel('Relative slippage $[\frac{rad}{s}]$', 'interpreter', 'latex', 'fontsize', 16);
xlabel('Time $[s]$', 'interpreter', 'latex', 'fontsize', 16);
legend('$\Delta\dot\theta(t)$', 'interpreter', 'latex', 'fontsize', 16, 'location', 'northwest');

% Plot gravity-filtered accelerometer signals
figure(i);i=i+1;
%
subplot(2,1,1);
plot(t_arr, d2r_E0_E0, 'linewidth', 2);
box on; grid on;
ylim([-3 3]);
xlim([0 t_plot]);
ylabel('Acceleration $E_0~[\frac{m}{s^2}]$', 'interpreter', 'latex', 'fontsize', 16);
xlabel('Time $[s]$', 'Interpreter','latex','FontSize',16);
legend('$\ddot{r}_{E_0}^{E_0}(t)\cdot {\bf i}_x^{E_0}$', '$\ddot{\bf r}_{E_0}^{E_0}(t)\cdot {\bf i}_y^{E_0}$', 'interpreter', 'latex', 'fontsize', 16, 'location', 'northwest', 'orientation', 'horizontal');
%
subplot(2,1,2);
plot(t_arr, d2r_E1_E0, 'linewidth', 2);
box on; grid on;
ylim([-3 3]);
xlim([0 t_plot]);
ylabel('Acceleration $E_1~[\frac{m}{s^2}]$', 'interpreter', 'latex', 'fontsize', 16);
xlabel('Time $[s]$', 'Interpreter','latex','FontSize',16);
legend('$\ddot{\bf r}_{E_1}^{E_0}(t)\cdot {\bf i}_x^{E_0}$','$\ddot{\bf r}_{E_1}^{E_0}(t)\cdot {\bf i}_y^{E_0}$', 'interpreter', 'latex', 'fontsize', 16, 'location', 'northwest', 'orientation', 'horizontal');

% Find indices corresponding to when the sub is slipping
counter=1;
sub_array_cell = {};
start_ind = 1;
true_flag_ind = find(flag_arr==true);
for element=2:length(true_flag_ind)
    if(true_flag_ind(element) - true_flag_ind(element-1) ~= 1)
        sub_array = true_flag_ind(start_ind:element-1);
        sub_array_cell{counter} = sub_array;
        counter=counter+1;
        start_ind = element;
    end
end
% Extract last sub-array
sub_array_cell{counter} = true_flag_ind(start_ind:end);

% Plot residual
figure(i);i=i+1;
hold on
for cell=sub_array_cell
    cell_arr = cell{:};
    rectangle(Position=[t_arr(cell_arr(1)),-1,t_arr(cell_arr(end))-t_arr(cell_arr(1)),6], FaceColor=[0.5 0 0 0.1], EdgeColor=[0.5 0 0 0.2]);
end
plot(t_arr, residual, 'linewidth', 2);
plot(t_arr, 0.25*ones(size(t_arr)), '--', 'linewidth', 2, 'color', 'black');
plot(t_arr, -0.25*ones(size(t_arr)), '--', 'linewidth', 2, 'color', 'black');
box on;grid on;
ylabel('Residual signal $[\frac{m}{s^2}]$', 'Interpreter','latex','FontSize',16);
xlabel('Time $[s]$', 'Interpreter','latex','FontSize',16);
legend('$\varrho(t)$', 'Threshold', 'interpreter', 'latex', 'fontsize', 16, 'location', 'northwest');
xlim([0 t_plot]);
ylim([-1 2]);