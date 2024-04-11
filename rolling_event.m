%ROLLING_EVENT Monitors for a rolling event.
%   Monitors for the transition from rolling with slipping to rolling
%   without slipping. For use with the events framework via 'odeset' in
%   specifying the ODE integration options.
function [value, isterminal, direction] = rolling_event(t,x,param)
% © Nils C.A. Wilhelmsen
% 11/04/2024
%% Extract parameters
r_So = param.r_So;
r_w = param.r_w;
%% Extract states
d_theta = x(2);
d_phi_S = x(4);
%% Compute alpha
alpha = -r_So/(r_w - r_So);
%% Calculate event condition
value = alpha*d_phi_S - d_theta;                                            % Detect when relative slippage reaches zero
isterminal = 1;                                                             % Stop the integration
direction = 0;                                                              % Any direction
end