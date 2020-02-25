%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% numerical_sde_cpp.py
% 
% Calculate sample path/estimate moments of a dynamical system defined by
% a SDE by using Runge_kutta + averaging on a fast C++ compiled routine.
%
% Alternatively, simulate a deterministic ODE.
%
% Plot output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numerically solve a langevin type SDE dY = a(t,Y)dt + b*dW_t, Y(t0) = Y0
%
% inputs: dt             -- time step size, e.g. 1e-5
%         t_interval     -- vector with [t_ini t_end]
%         num_traces     -- Number of traces for average to estimate moments
%                           e.g. 100
%         subs_f         -- subsampling factor (for every N generated plots
%                           in C++ routine, import only N/subs_f)
%         eq_type        -- Choose between ode and sde simulation (i.e., 
%                           only drift or drift & diffusion).
%         Ito            -- Boolean type. True if Ito equation, false if
%                           Stratonovich equation.
%         n_dim          -- problem dimension, e.g. 2
%         initial_values -- vector with initial values, e.g. [0;0]
%
% output: tt             -- Time vector
%         avg_var        -- E(y^2) estimated from num_traces runs 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GP Conangla, 26.09.2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[tt, fun_avg, fun_var] = numerical_sde_cpp(dt, t_interval, num_traces, subs_f, eq_type, Ito, n_dim, initial_values)

%% Initial checks/allocations
% Check if folder where trace is going to be exists
if ~exist("./simulated_traces", 'dir')
    mkdir("./simulated_traces")
end

% Put initial values to a string to call ./sde.out routine
ini_values_string = "";
for i = 1:n_dim
    ini_values_string = ini_values_string + num2str(initial_values(i)) + " ";
end

%% Print info on stdout
% Info on the method
fprintf("====================================================================\n");
fprintf("Estimate moments of SDE with modified Runge_kutta of strong order 1.\n\n");

% delete files with old traces (if it already exists)
if isfile('./simulated_traces/sde_sample_path_avg_0.txt')
    delete './simulated_traces/sde_sample_path_avg_0.txt';
    delete './simulated_traces/sde_sample_path_avg_1.txt';
    delete './simulated_traces/sde_sample_path_var_0.txt';
    delete './simulated_traces/sde_sample_path_var_1.txt';
end

% recompile + link the eq_parameters.cpp file
fprintf("Updating drift and diffusion equations...\n");
cd ../C++/
status = system("g++ -c eq_definitions.cpp");
status = system("g++ -o sde.out sde.o sde_library.o num_vector.o physics.o eq_definitions.o -O3 -pthread");
cd ../MATLAB/

% call compiled C++ routine to generate traces
fprintf("Calling C++ routine...\n");
status = system("../C++/sde.out " + ...
    num2str(dt) + " " + num2str(t_interval(1)) + " " + num2str(t_interval(2)) + " " ...
    + num2str(num_traces) + " " + num2str(subs_f) + " " + eq_type + " " + num2str(Ito) + ...
    " " + num2str(n_dim) + " " + ini_values_string);

if status == 1 % check for errors
    tt = NaN;
    fun_avg = NaN;
    fun_var = NaN;
    fprintf("Exiting...\n\n");
    return;
end

%% import traces from files
fprintf("\nImporting traces...\n");
tic;
ini_t = toc;

fun_avg_x = importdata("./simulated_traces/sde_sample_path_avg_0.txt");
fun_avg_v = importdata("./simulated_traces/sde_sample_path_avg_1.txt");
fun_var_x = importdata("./simulated_traces/sde_sample_path_var_0.txt");
fun_var_v = importdata("./simulated_traces/sde_sample_path_var_1.txt");

% 1st moment
fun_avg = [fun_avg_x; fun_avg_v];

% 2nd moment (not yet variance)
sec_mom = [fun_var_x; fun_var_v];

% Convert to variance
fun_var = sec_mom - fun_avg.^2;

end_t = toc;
fprintf("Elapsed time to import traces: %f s\n", end_t - ini_t);
fprintf("\n");

%% Plot results
fprintf("Plotting results.\n\n");
% Generate time vector
tt = [t_interval(1):dt*subs_f:t_interval(2)];
tt = tt(1:length(fun_avg_x)); % Check that both vectors are of same length

% Plot x
figure(1);
clf;
hold on;
% Plot simulated variance
nice_plot(tt, fun_avg_x, "time (s)", "$\mathbf{E}[f(x(t))]$", "Simulated variance");

if eq_type == "sde" && num_traces > 1 % add 1 sigma standard error interval
    upper = fun_avg_x + sqrt(fun_var(1,:))/sqrt(num_traces);
    down = fun_avg_x - sqrt(fun_var(1,:))/sqrt(num_traces);
    shade(tt, upper, tt, down, 'FillType',[1 2;2 1], 'FillColor', [0 0 0], 'FillAlpha', 0.1);
end

% Plot v
figure(2);
clf;
hold on;
% Plot simulated variance
nice_plot(tt, fun_avg_v, "time (s)", "$\mathbf{E}[f(v(t))]$", "Simulated variance");

if eq_type == "sde" && num_traces > 1 % add 1 sigma standard error interval
    upper = fun_avg_v + sqrt(fun_var(2,:))/sqrt(num_traces);
    down = fun_avg_v - sqrt(fun_var(2,:))/sqrt(num_traces);
    shade(tt, upper, tt, down, 'FillType',[1 2;2 1], 'FillColor', [0 0 0], 'FillAlpha', 0.1);
% end
end



















