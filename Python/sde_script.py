###########################################################################
# python script: call numerical_sde_cpp.py (simulate trace + plot result)
###########################################################################
# Numerically solve a langevin type SDE dY = a(t,Y)dt + b*dW_t, Y(t0) = Y0
#
#         dt             -- time step size, e.g. 1e-5
#         t_interval     -- vector with [t_ini t_end]
#         initial_values -- vector with initial values, e.g. [0; 0],
#                           can also be 'rand_ini'
#         subs_f         -- Subsampling factor. For every N points generated
#                           by the C++ routine, import only N/subs_f
#         num_traces     -- Number of traces for average to estimate moments
#                           e.g. 100
#         n_dim          -- problem dimension, e.g. 2
#         initial_values -- vector with initial values, e.g. [0;0]
###########################################################################

import numerical_sde_cpp as sde

dt = 1e-8 # seconds
t_interval = [0, 2e-4]
num_traces = 1000
subs_f = 40
Ito = 1
n_dim = 2 # 2D: x and v
initial_values = [0, 0] # start at origin

# call function
sde.numerical_sde(dt, t_interval, num_traces, subs_f, Ito, n_dim, initial_values)
