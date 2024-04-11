%SLIPPING_EVENT Monitors for a slipping event.
%   Monitors for the transition from rolling without slipping to rolling
%   with slipping. For use with the events framework via 'odeset' in
%   specifying the ODE integration options.
function [value, isterminal, direction] = slipping_event(t, x, param, k_t, F_e)
% © Nils C.A. Wilhelmsen
% 11/04/2024
%% Extract parameters
m = param.m;
I_S = param.I_S;
r_So = param.r_So;
r_Si = param.r_Si;
r_w = param.r_w;
g = param.g;
mu_s = param.mu_s;
rho = param.rho;
V = param.V;
%% Extract states
theta = x(1);
d_theta = x(2);
d_phi_S = x(4);
d_phi_ds = x(6);
%% Compute friction and normal forces
d2_theta = (1/((m + I_S/(r_So^2))*(r_w - r_So))*(-(m-rho*V)*g*sin(theta) + F_e - (r_Si/r_So)*k_t*(d_phi_ds - d_phi_S)));
F_f = (-I_S/(r_So^2))*(r_w - r_So)*d2_theta - (r_Si/r_So)*k_t*(d_phi_ds - d_phi_S);
F_N = (m-rho*V)*g*cos(theta) + m*(r_w - r_So)*d_theta^2;
%% Calculate event condition
value = abs(F_f) - mu_s*abs(F_N);                                           % Condition for transitioning to slipping
isterminal = 1;                                                             % Stop the integration
direction = 1;                                                              % Slippage occurs when 'value' crosses the origin from negative to positive value
end