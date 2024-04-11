%ROLLING_MOTION Computes time-derivative of sub motion for simulation.
%   Given net external torque 'tau_e' and force 'F_e' acting on the
%   drillstring segment along with magnetic damping coefficient 'k_t', the
%   first- and second-order time derivatives of the sub angle in wellbore
%   'theta', sub angle of rotation 'psi_S' and drillstring angle of
%   rotation 'phi_ds' are calculated for use in a MATLAB ODE solver. The
%   function can also be used as a building block in a subroutine for
%   computing the time-derivative of a larger system in which the sub is a
%   component.
function dx = rolling_motion(t,x,param,slipping,k_t,tau_e,F_e)
% © Nils C. A. Wilhelmsen
% 11/04/2024
%% Extract model parameters
m = param.m;                                                                % Combined sub and drillstring segment mass [kg]
V = param.V;                                                                % Net sub and drillstring segment volume [m^3]
rho = param.rho;                                                            % Drilling mud density [kg/m^3]
g = param.g;                                                                % Acceleration of gravity [m/s^2]
r_w = param.r_w;                                                            % Radius of wellbore [m]
r_So = param.r_So;                                                          % Outer radius of sub [m]
r_Si = param.r_Si;                                                          % Inner radius of sub [m]
r_ds = param.r_ds;                                                          % Outer radius of drillstring [m]
I_S = param.I_S;                                                            % Sub moment of inertia [kg*m^2]
I_ds = param.I_ds;                                                          % Drillstring moment of inertia [kg*m^2]
mu_k = param.mu_k;                                                          % Kinetic friction coefficient [-]
alpha = param.alpha;                                                        % Proportionality constant during rolling without slipping
%% Extract states
theta = x(1);
d_theta = x(2); 
phi_S = x(3);
d_phi_S = x(4); 
phi_ds = x(5); 
d_phi_ds = x(6);
%% Compute second-order derivatives
if(slipping)
    % Compute relative slippage
    alpha = -r_So/(r_w - r_So);
    Delta_dot_theta = alpha*d_phi_S - d_theta;

    % Compute F_f
    F_f = mu_k*abs((m-rho*V)*g*cos(theta) + m*(r_w - r_So)*d_theta^2)*sign(Delta_dot_theta);
    
    % Compute d2_theta, d2_psi_S in case of slipping
    d2_theta = (1/(m*(r_w - r_So)))*(-(m-rho*V)*g*sin(theta) + F_f + F_e);
    d2_phi_S = (1/I_S)*(r_Si*k_t*(d_phi_ds - d_phi_S) + r_So*F_f);
else
    % Compute d2_theta, d2_psi_S in case of not slipping
    d2_theta = (1/((m + (I_S/(r_So^2)))*(r_w - r_So)))*(-(m-rho*V)*g*sin(theta) + F_e - (r_Si/r_So)*k_t*(d_phi_ds - d_phi_S));
    d2_phi_S = (1/(I_S + m*r_So^2))*(r_Si*k_t*(d_phi_ds - d_phi_S) + r_So*(m-rho*V)*g*sin(theta) - r_So*F_e);
end
% Compute d2_phi_ds
d2_phi_ds = (1/I_ds)*(r_ds*k_t*(d_phi_S - d_phi_ds) + tau_e);
%% Set up and return dx
dx = [d_theta; d2_theta; d_phi_S; d2_phi_S; d_phi_ds; d2_phi_ds;];
end