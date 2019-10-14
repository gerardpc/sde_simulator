%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot_theory.m
% 
% Compare simulated trace variance with analytical expression of variance
% to check the code. This is just an example of a free Brownian particle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[] = plot_theory(tt)

% Compare with theory
% Experimental parameters
m = 1.495e-17;
gamma = 4.075e-11;
T = 295;
k_B = 1.38065e-23;
sigma_noise = sqrt(2*k_B*T*gamma);
w = 128e3*2*pi;

%% Stochastic harmonic oscillator, 2nd order
% C = w^2;
% A = gamma/m;
% B = sigma_noise/m;
% C2 = (C - A^2/4);
% 
% % f(1) = v; 
% % f(2) = -C*x - A*v;
% % sigma = [0; B]
% 
% % Then, variance of x will be
% var_short = B^2/(2*A*C)*(1 - exp(-A*tt).*(C/C2 - A^2/(4*C2)*cos(2*sqrt(C2)*tt) + A/(2*sqrt(C2))*sin(2*sqrt(C2)*tt)));
% % Plot variance expression
% nice_plot(tt, var_short, "", "", "");

%% PAUL TRAP
var_short_t = sigma_noise^2/gamma^2*tt;
nice_plot(tt, var_short_t, "", "", "");
end