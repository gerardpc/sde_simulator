%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot_theory.m
% 
% Compare simulated trace variance with analytical expression of variance
% to check the code. This is just an example of a free Brownian particle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[] = plot_theory(tt)

% Compare with theory
% Experimental parameters
m = 9.2e-18;
gamma = 3.5e-11;
T = 295;
k_B = 1.38065e-23;
eps = 6.3e-9;
sigma_noise = sqrt(2*k_B*T*gamma);

% Plot variance expression
nice_plot(tt, (sigma_noise/m)^2*tt, "", "", "");
end