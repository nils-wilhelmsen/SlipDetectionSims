%DRILLSTRING_SUB_COUPLED_MODEL Computes time-derivative for coupled model.
%   Couples low-order drillstring model to model describing rolling motion.
function dx = drillstring_sub_coupled_model(t,x,param,slipping,k_t, F_e)
% © Nils C. A. Wilhelmsen
% 11/04/2024
%% Extract model parameters
I_td = param.I_td;                                                          % Top-drive moment of inertia [kg*m^2]
c = param.c;                                                                % Drillstring rotational spring constant [...]
n = param.n;                                                                % Number of blades in drag bit [-]
a = param.a;                                                                % Radius of drag bit [m]
xi = param.xi;                                                              % Parameter characterizing spatial distribution of wearflats [-]
mu = param.mu;                                                              % Bit friction coefficient [-]
l = param.l;                                                                % Length of wearflats [m]
sigma_bar = param.sigma_bar;                                                % Bit contact stress [...]
t_drill = param.t_drill;                                                    % Time at which bit comes in contact with rocks [s]
k_P = param.k_P;                                                            % Top drive controller P gain [N*m*s/rad]
k_I = param.k_I;                                                            % Top drive controller I gain [N*m/rad]
d_phi_td_setpoint = param.d_phi_td_setpoint;                                % Top drive velocity setpoint [rad/s]
%% Extract paramters
phi_ds = x(5);
phi_td = x(7);
d_phi_td = x(8);
integrator = x(9);
%% Compute torques and forces acting on system
tau_d = k_P*(d_phi_td_setpoint-d_phi_td) + k_I*integrator;                  % Top drive torque
tau_f = 0.5*n*a^2*xi*mu*l*sigma_bar*heaviside(t-t_drill);                   % Frictional torque on bit
tau_e = c*(phi_td - phi_ds) - tau_f;                                        % Net external torque on drillstring
%% Compute derivatives
% Rolling motion model time derivatives
dx_sub = rolling_motion(t,x,param,slipping,k_t,tau_e,F_e);
% Second-order time derivative of top drive
d2_phi_td = (1/I_td)*(c*(phi_ds - phi_td) + tau_d);
% First-order time-derivative of integrator
d_integrator = d_phi_td_setpoint-d_phi_td;
%% Set up and return dx
dx = [dx_sub;d_phi_td;d2_phi_td;d_integrator];
end