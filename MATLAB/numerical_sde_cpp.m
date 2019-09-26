%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% numerical_sde_cpp.m
% 
% Calculate sample path/estimate moments of SDE by using Runge_kutta + 
% averaging on a fast C++ compiled routine.
% (x500 speed increase with respect to equivalent MATLAB routine).
% Plot output with nice_plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numerically solve a langevin type SDE dY = a(t,Y)dt + b*dW_t, Y(t0) = Y0
%
% inputs: dt             -- time step size, e.g. 1e-5
%         t_interval  -- vector with [t_ini t_end]
%         num_traces     -- Number of traces for average to estimate moments
%                           e.g. 100
%         n_dim          -- problem dimension, e.g. 2
%         initial_values -- vector with initial values, e.g. [0;0]
%
% output: tt             -- Time vector
%         avg_var        -- E(y^2) estimated from num_traces runs 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GP Conangla, 26.09.2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[tt, avg_var] = numerical_sde_cpp(dt, t_interval, num_traces, n_dim, initial_values)

%% Initial checks/allocations
% Check if folder where trace is going to be exists
if ~exist("./simulated_traces", 'dir')
    mkdir("./simulated_traces")
end

% Sanity check: num_traces not multiple of 4, make it
num_traces = num_traces - mod(num_traces, 4);

% Put initial values to a string to call ./sde.out routine
ini_values_string = "";
for i = 1:n_dim
    ini_values_string = ini_values_string + num2str(initial_values(i)) + " ";
end

%% Print info on stdout
% Info on the method
fprintf("====================================================================\n");
fprintf("Estimate moments of SDE with modified Runge_kutta of strong order 1.\n");
% print number of traces on stdout
fprintf("Number of traces that will be generated: %d\n\n", num_traces);

% delete files with old traces (if it already exists)
if isfile('./simulated_traces/sde_sample_path_0.txt')
    delete './simulated_traces/sde_sample_path_0.txt';
    delete './simulated_traces/sde_sample_path_1.txt';
end

% call compiled C++ routine to generate traces
fprintf("Calling C++ routine...\n");
status = system("../C++/sde.out " + ...
    num2str(dt) + " " + num2str(t_interval(1)) + " " + num2str(t_interval(2)) + " " ...
    + num2str(num_traces) + " " + num2str(n_dim) + " " + ini_values_string);

%% import traces from files
fprintf("\nImporting traces...\n");
ini_t = toc;
avg_var_x = importdata("./simulated_traces/sde_sample_path_0.txt");
avg_var_v = importdata("./simulated_traces/sde_sample_path_1.txt");
avg_var = [avg_var_x; avg_var_v];
end_t = toc;
fprintf("Elapsed time to import traces: %f s\n", end_t - ini_t);
fprintf("\n");

%% Plot results
fprintf("Plotting results.\n\n");
% Generate time vector
tt = [t_interval(1):dt:t_interval(2)];
tt = tt(1:length(avg_var_x)); % Check that both vectors are of same length

% Plot x
figure(1);
clf;
hold on;
nice_plot(tt, avg_var_x, "time (s)", "$\mathbf{E}[x^2(t)]$", "Simulated variance");

% Plot v
figure(2);
clf;
hold on;
% Plot theory comparison
plot_theory(tt);
% Plot simulated variance
nice_plot(tt, avg_var_v, "time (s)", "$\mathbf{E}[v^2(t)]$", "Simulated variance");

end



















