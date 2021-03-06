###########################################################################
# numerical_sde_cpp.py
# 
# Calculate sample path/estimate moments of a dynamical system defined by
# a SDE by using Runge_kutta + averaging on a fast C++ compiled routine.
#
# Alternatively, simulate a deterministic ODE.
#
# Plot output
###########################################################################
# Numerically solve a langevin type SDE dY = a(t,Y)dt + b*dW_t, Y(t0) = Y0
#
# inputs: dt             -- time step size, e.g. 1e-5
#         t_interval     -- vector with [t_ini t_end]
#         num_traces     -- Number of traces for average to estimate moments
#                           e.g. 100
#         eq_type        -- Equation type: "ode" (deterministic) or 
#                           "sde" (stochastic)
#         subs_f         -- Subsampling factor. For every N points generated
#                           by the C++ routine, import only N/subs_f
#         Ito            -- Boolean type. True if Ito equation, false if
#                           Stratonovich equation.
#         n_dim          -- problem dimension, e.g. 2
#         initial_values -- vector with initial values, e.g. [0;0]
#
# output: tt             -- Time vector
#         avg_var        -- E(y^2) estimated from num_traces runs 
#
###########################################################################
# GP Conangla, 09.10.2019
###########################################################################

import subprocess
import os
import time
import numpy as np
import matplotlib.pyplot as plt


# Calculate analytical expression vector for comparison
def plot_theory(tt):
    # Experimental parameters (just an example, free Brownian particle)
    m = 1.495e-17
    gamma = 4.075e-11
    T = 295
    k_B = 1.38065e-23
    sigma_noise = (2*k_B*T*gamma)**0.5
    C = (110e3*2*np.pi)**2
    A = gamma/m
    B = sigma_noise/m
    C2 = (C - A**2/4)

    # Return analytical variance expression
    var_analytical = B**2/(2*A*C)*(1 - np.exp(-A*tt)*(C/C2 - A**2/(4*C2)*np.cos(2*C2**0.5*tt)
           + A/(2*C2**0.5)*np.sin(2*C2**0.5*tt)))
    return var_analytical


# main script to generate traces and plot result
def numerical_sde(dt, t_interval, num_traces, subs_f, eq_type, Ito, n_dim, initial_values):
    ####################################################################
    # Initial checks/allocations
    # Check if folder where trace is going to be exists
    if not os.path.exists("./simulated_traces"):
        # if not, create it
        os.makedirs("./simulated_traces")

    # Put initial values to a string to call ./sde.out routine
    ini_values_string = ""
    for i in range(n_dim):
        ini_values_string = ini_values_string + str(initial_values[i]) + " "

    # Print info on stdout
    # Info on the method
    print("====================================================================")
    print("Estimate moments of SDE with modified Runge_kutta of strong order 1.\n")

    # delete files with old traces, if there are any
    if os.path.exists("./simulated_traces/sde_sample_path_avg_0.txt"):
        for i in range(n_dim):
            os.remove("./simulated_traces/sde_sample_path_avg_" + str(i) + ".txt")
            os.remove("./simulated_traces/sde_sample_path_var_" + str(i) + ".txt")

    ####################################################################
    # Recompile + link drift and diffusion equations
    print("Updating drift and diffusion equations...")
    os.chdir("../C++")
    subprocess.call("g++ -c eq_definitions.cpp", shell=True)
    subprocess.call("g++ -o sde.out sde.o sde_library.o num_vector.o physics.o eq_definitions.o -O3 -pthread", shell=True)
    os.chdir("../Python")

    ####################################################################
    # call compiled C++ routine to generate traces
    print("Calling C++ routine...\n")
    cpp_command = "../C++/sde.out " + str(dt) + " " + str(t_interval[0]) + " " + str(t_interval[1]) + " " + str(num_traces) + " " + str(subs_f) + " " + eq_type + " " + str(Ito) + " " + str(n_dim) + " " + ini_values_string
    status = subprocess.call(cpp_command, shell = True)

    # check that there were no errors
    if status != 0:
        print("Exiting...\n")
        return

    # import traces from files
    print("\nImporting traces...")
    ini_t = time.time()

    avg_f = np.loadtxt("./simulated_traces/sde_sample_path_avg_0.txt")
    avg_f = avg_f.reshape((avg_f.size, 1))

    for i in range(1, n_dim):
        new_row = np.loadtxt("./simulated_traces/sde_sample_path_avg_" + str(i) + ".txt")
        new_row = new_row.reshape((new_row.size, 1))
        avg_f = np.concatenate((avg_f, new_row), axis=1)

    var_f = np.loadtxt("./simulated_traces/sde_sample_path_var_0.txt")
    var_f = var_f.reshape((var_f.size, 1))
    for i in range(1, n_dim):
        new_row = np.loadtxt("./simulated_traces/sde_sample_path_var_" + str(i) + ".txt")
        new_row = new_row.reshape((new_row.size, 1))
        var_f = np.concatenate((var_f, new_row), axis=1)

    # Transpose to maintain coherence with rest of the code
    avg_f = np.transpose(avg_f)
    var_f = np.transpose(var_f)

    # calculate variance
    # the absolute is just a sanity check, it should be always positive
    # but I've found some rounding errors that lead to negative results
    var_f = np.absolute(var_f - np.square(avg_f))
    end_t = time.time()
    print("Elapsed time to import traces: %f s\n" % (end_t - ini_t))

    ####################################################################
    # Plot results
    print("Plotting results.\n")
    # Generate time vector
    tt = np.arange(t_interval[0], t_interval[1], dt*subs_f)

    if len(tt) != len(avg_f[0, :]):  # Check that both vectors are of same length
        tt = tt[0:len(avg_f[0, :])]
        avg_f = avg_f[:, 0:len(tt)]
        var_f = var_f[:, 0:len(tt)]

    x_dim = 0
    v_dim = 1

    # Plot 0
    fig0, ax = plt.subplots()

    # Plot theory comparison
    # yy_theory = plot_theory(tt)    
    # ax.plot(tt, yy_theory)

    # Plot uncertainty bounds if sde
    if eq_type == "sde" and num_traces > 1:
        # standard error bounds
        upper_bound = avg_f[x_dim, :] + np.sqrt(var_f[x_dim, :]/num_traces)
        lower_bound = avg_f[x_dim, :] - np.sqrt(var_f[x_dim, :]/num_traces)
        ax.fill_between(tt, lower_bound, upper_bound, color = "grey", alpha = 0.2)

    # Plot simulated stat. moment
    ax.plot(tt, avg_f[x_dim, :])

    ax.grid(True)

    plt.xlabel("time (s)", fontsize = 18)
    plt.ylabel("$\mathbb{E}[f(x(t))]$", fontsize = 18)
    plt.title("Simulated variance", fontsize = 18)

    ax.set_xlim([min(tt), max(tt)])
    ax.set_ylim([min(avg_f[x_dim, :]), max(avg_f[x_dim, :]*1.1)])
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

    fig0.savefig("./test0.svg")

    # Plot 1
    fig1, ax = plt.subplots()

    # Plot uncertainty bounds if sde
    if eq_type == "sde" and num_traces > 1:
        # standard error bounds
        upper_bound = avg_f[v_dim, :] + np.sqrt(var_f[v_dim, :]/num_traces)
        lower_bound = avg_f[v_dim, :] - np.sqrt(var_f[v_dim, :]/num_traces)
        ax.fill_between(tt, lower_bound, upper_bound, color = "grey", alpha = 0.2)

    # Plot simulated stat. moment
    ax.plot(tt, avg_f[v_dim, :])
    ax.grid(True)

    plt.xlabel("time (s)", fontsize = 18)
    plt.ylabel("$\mathbb{E}[f(v(t))]$", fontsize = 18)
    plt.title("Simulated variance", fontsize = 18)

    ax.set_xlim([min(tt), max(tt)])
    ax.set_ylim([min(avg_f[v_dim, :]), max(avg_f[v_dim, :]*1.1)])
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

    fig1.savefig("./test1.svg")


















