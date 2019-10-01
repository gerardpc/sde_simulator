###########################################################################
# numerical_sde_cpp.py
# 
# Calculate sample path/estimate moments of SDE by using Runge_kutta + 
# averaging on a fast C++ compiled routine.
# (x500 speed increase with respect to equivalent MATLAB routine).
# Plot output with nice_plot
###########################################################################
# Numerically solve a langevin type SDE dY = a(t,Y)dt + b*dW_t, Y(t0) = Y0
#
# inputs: dt             -- time step size, e.g. 1e-5
#         t_interval     -- vector with [t_ini t_end]
#         num_traces     -- Number of traces for average to estimate moments
#                           e.g. 100
#         n_dim          -- problem dimension, e.g. 2
#         initial_values -- vector with initial values, e.g. [0;0]
#
# output: tt             -- Time vector
#         avg_var        -- E(y^2) estimated from num_traces runs 
#
###########################################################################
# GP Conangla, 26.09.2019
###########################################################################

import subprocess
import os
import time
import numpy
import matplotlib.pyplot as plt

# Calculate analytical expression vector for comparison
def plot_theory(tt):
    # Experimental parameters (just an example, free Brownian particle)
    m = 9.2e-18
    gamma = 3.5e-11
    T = 295
    k_B = 1.38065e-23
    sigma_noise = (2*k_B*T*gamma)**(1/2)
    # Return analytical variance expression
    return (sigma_noise/m)**2*tt

# main script to generate traces and plot result
def numerical_sde(dt, t_interval, num_traces, n_dim, initial_values):

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
    if os.path.exists("./simulated_traces/sde_sample_path_0.txt"):
        for i in range(n_dim):
            os.remove("./simulated_traces/sde_sample_path_" + str(i) + ".txt")

    ####################################################################
    # call compiled C++ routine to generate traces
    print("Calling C++ routine...\n")
    cpp_command = "../C++/sde.out " + str(dt) + " " + str(t_interval[0]) + " " + str(t_interval[1]) + " " + str(num_traces) + " " + str(n_dim) + " " + ini_values_string
    status = subprocess.call(cpp_command, shell = True)
    
    # import traces from files
    print("\nImporting traces...")
    ini_t = time.time()
    
    avg_var = numpy.loadtxt("./simulated_traces/sde_sample_path_0.txt")
    for i in range(1, n_dim):
        new_row = numpy.loadtxt("./simulated_traces/sde_sample_path_" + str(i) + ".txt")
        avg_var = numpy.stack((avg_var, new_row))
        
    end_t = time.time()
    print("Elapsed time to import traces: %f s\n" % (end_t - ini_t))
    
    ####################################################################
    # Plot results
    print("Plotting results.\n")
    # Generate time vector
    tt = numpy.arange(t_interval[0], t_interval[1], dt)
    
    if len(tt) != len(avg_var[0,:]): # Check that both vectors are of same length
        tt = tt[0:len(avg_var[0,:])]
        avg_var = avg_var[:, 0:len(tt)]
    
    x_dim = 0;
    v_dim = 1;
    
    # Plot 0
    fig0, ax = plt.subplots()
    
    # Plot theory comparison
    yy_theory = plot_theory(tt)
    
    # Plot simulated variance
    ax.plot(tt, avg_var[x_dim, :])
    ax.grid(True)
    
    plt.xlabel("time (s)", fontsize = 18)
    plt.ylabel("$\mathbb{E}(x^2(t))$", fontsize = 18)
    plt.title("Simulated variance", fontsize = 18)  
    
    ax.set_xlim([min(tt), max(tt)])  
    ax.set_ylim([min(avg_var[x_dim, :]), max(avg_var[x_dim, :]*1.1)])    
    
    fig0.savefig("./test0.png", dpi = 1000)
    
    # Plot 0
    fig1, ax = plt.subplots()
    
    # Plot theory comparison
    yy_theory = plot_theory(tt)
    
    # Plot simulated variance
    ax.plot(tt, yy_theory)
    ax.plot(tt, avg_var[v_dim, :])
    ax.grid(True)
    
    plt.xlabel("time (s)", fontsize = 18)
    plt.ylabel("$\mathbb{E}(v^2(t))$", fontsize = 18)
    plt.title("Simulated variance", fontsize = 18)  
    
    ax.set_xlim([min(tt), max(tt)])  
    ax.set_ylim([min(avg_var[v_dim, :]), max(avg_var[v_dim, :]*1.1)])    
    
    fig1.savefig("./test1.png", dpi = 1000)


















