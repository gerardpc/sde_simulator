%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% python script: call numerical_sde_cpp.py (simulate trace + plot result)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters:
%
%         dt             -- time step size, e.g. 1e-5
%         t_interval     -- vector with [t_ini t_end]
%         initial_values -- vector with initial values, e.g. [0; 0],
%                           can also be 'rand_ini'
%         subs_f         -- Subsampling factor. For every N points generated
%                           by the C++ routine, import only N/subs_f
%         eq_type        -- Equation type: "ode" (deterministic) or 
%                           "sde" (stochastic)
%         num_traces     -- Number of traces for average to estimate moments
%                           e.g. 100
%         n_dim          -- problem dimension, e.g. 2
%         initial_values -- vector with initial values, e.g. [0;0]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dt = 2*1e-7; % seconds
t_interval = [0, 1e-3];
num_traces = 100;
eq_type = "sde";
subs_f = 1;
Ito = 1;
n_dim = 2; % 2D: x and v
initial_values = [1e-7, 0]; % start at origin

[tt, fun_avg, fun_var] = numerical_sde_cpp(dt, t_interval, num_traces, subs_f, eq_type, Ito, n_dim, initial_values);
min(min(fun_var))
max(max(fun_var))