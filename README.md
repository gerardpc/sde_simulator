# SDE_SIMULATOR project

Description
--------------
This project includes functions and libraries to numerically simulate a 
dynamical system defined by the SDE
  
dY = a(t,Y)dt + b(t,Y)*dW_t, Y(t0) = Y0,

where a(t,Y) (drift term) and b(t,Y) (diffusion term) are arbitrary (i.e.,
defined by the user) and Y can be of arbitrary (but finite) dimension.

Two simulation options are allowed:
 
       - Deterministic dynamical system (ODE). No noise in the equations of 
         motion is considered. In this case, an Adams predictor-corrector 
         method of 4th order is used.
 
       - Stochastic dynamical system (SDE): a Runge-Kutta type method for SDE
         of strong order 1 and deterministic order 2 is used. 
         Since every run of the method for the SDE will generate a different 
         sample path, the code is optimized to generate many sample paths
         at once using multithreading. When the simulation phase is finished,
         the averages and standard errors are calculated.

Assuming the "sde" simulation type is chosen, a typical application 
would consist of generating n sample paths Y_t, 
computing a certain function of them, f(Y_t), and then averaging the result, to find 
the expected value of f(Y_t). The most common function choice is the 2nd 
moment, the variance of Y(t): <Y^2(t)>.

The code uses a stochastic Runge-Kutta type method of strong order 1 (in 
the presence of noise) that converges to deterministic order 2 (i.e., for zero noise) 
that does not require any non-zero derivatives of b (diffusion term).
Other methods (e.g. the Milstein method) have strong order 1 but reduce
to the Euler-Maruyama method (strong order 0.5) when b is a constant. 
The method also allows choosing between Ito or Stratonovich interpretations
of the SDE.

The code is divided in two parts:

- First: A C++ main() + library functions for the most numerically intensive 
part: in the deterministic case, numerically integrate the ODE; in
the stochastic case, generate the traces Y_t, calculate f(Y_t) and f^2(Y_t) 
(to be able to calculate the standard error) and average the traces. 
Then write the result in a file. Note that in the present form, the C++
main calculates (by default) the VARIANCE, i.e. f(Y_t) = Y_t^2, and not 
any arbitrary function. The function whose expected value is calculated
is defined in num_vector.cpp as "double f(double f)" and can be changed
to any other function the user wants to define. 
The code uses multithreading (detects number of cores and adapts number of
threads).

- Second: MATLAB and Python scripts to load the saved traces and plot the results.

In my personal computer, this code boosted the performance (reduction of
execution time) a factor ~1000 with respect to equivalent code implemented
100% in MATLAB.

Everything is well commented, so you should be able to modify the code 
(whether the C++ or the MATLAB/Python part) easily and perform your own 
numerical analysis.

 The deterministic ODE Adams method can be found in
 > https://en.wikiversity.org/wiki/Adams-Bashforth_and_Adams-Moulton_methods

The RK stochastic method is an implementation of 
> A. J. Roberts. Modify the improved Euler scheme to integrate stochastic differential equations. [1], Oct 2012.
which you can find here: https://arxiv.org/abs/1210.0933

If you use this code (or a modified version of it) in a paper, I ask that you cite 
> Conangla, Gerard P., et al. "Paper in preparation"
in the references section.

For any other comments, suggestions or questions, drop me a message.


Guide:
--------------
Clone the project to a local folder. I suggest using a Linux distribution,
and the guide assumes this is the case.


1- Compile the C++ code
--------------
Open a terminal, go to the C++ folder

    >> cd C++

Compile the code with the following command

    >> g++ sde.cpp sde_library.cpp num_vector.cpp physics.cpp eq_definitions.cpp -o sde.out -O3 -pthread

The option -std=c++11 is in general not needed (your compiler probably 
uses the C++11 standard by default). If you have trouble compiling, try
it again by adding this option.


2- Generate traces and plot the results
--------------
The following is a list of the C++ routine parameters:

1. dt is the discretization time step. Usually it is enough to use a dt 
small enough for the method to be stable (requires experimenting a bit).
2. t_int, an interval (e.g. [0 1]) in seconds to simulate.
3. num_traces, the number of sample paths that will be simulated and averaged.
This parameter is not used in the ODE case: uses 1 by default.
4. subsampling_f, an integer number greater or equal than 1 that imports
only one every subsampling_f points. Makes importing a file much faster
for large traces.
5. eq_type, a string with two possible values ("ode" = deterministic simulation,
"sde" = stochastic simulation).
6. Ito, a boolean that selects an Ito SDE when is true and a Stratonovich
SDE when false. This parameter is not used in the ODE case.
7. n_dim, the equation dimension (i.e., the number of variables or equations).
8. y0, a vector with the initial conditions for each of the n_dim degrees 
of freedom.

If you are working with MATLAB:

Open MATLAB, go to the MATLAB folder. Call the script "sde_script" from
the MATLAB Command Window:
    
    >> sde_script;
    
If you inspect the code, all the parameters are detailed and can be changed.
The results will be plotted in different figures.

If you are working with Python:

Open a terminal, go to the Python folder. Call the function "sde_script.py".
If you inspect the code, all the parameters are detailed and can be changed.
E.g.

    >> python sde_script.py

The results will be saved as svg figures. Check the "test.svg" file for
an example output.


2- Changing the a(t,Y) and b(t,Y) functions of the SDE
--------------
The a(t,Y) and b(t,Y) functions define the dynamical system. Their expression
can be modified in the file eq_definitions.cpp as the "drift_function" and
"diffusion_function", respectively.

The f[i] vector inside each of these functions contains the driving
term for the equations of motion (the drift_function is the "deterministic"
part and the diffusion_function the "stochastic" part). 

For instance, consider a Brownian particle completely driven by thermal 
noise (i.e., not subject to any force potential). If the particle only 
moves in one dimension, the equations will be

    x' = v
    v' = 0 + sigma/m

The first equation is the definition of the velocity, and the 2nd equation 
is Newton's 2nd law. This is precisely the equation included in the repository 
by default: inspect the code to understand the syntax.
Following this example, other, arbitrary, drift and diffusion functions
can be written. Of course, every time any of these functions is changed, the
C++ code has to be recompiled with the command indicated in section 1.

3- Changing the function f(x) that is used in the expected value, <f(Y_t)>
--------------
The f(Y_t) function is defined in the file num_vector.cpp as
"double f(double x)". Replacing the function, which by default returns 
"return x*x;" (to calculate the variance) by any other is straightforward.
Just remember to write it in C++.

4- If you want to use this code to simulate Brownian particles...
--------------
The code is specifically designed to be easy to do that (although of 
course it is not mandatory; any drift/diffusion functions can be used). 
You should use the struct "eq_params", defined in physics.cpp/physics.hpp,
which uses the parameters defined in the folder "eq_params" and the functions 
and struct members from particle_traps.cpp and particle_traps.hpp.
The values from this struct can then be used to define the drift and diffusion,
since eq_params is passed as one of the function parameters.

eq_params has five substructs:

1. Gaussian beam: defines a Gaussian beam parameters. Its values can be
modified in ./eq_params/gaussian_beam.txt.

2. Particle: defines the trapped particle parameters. Its values can be
modified in ./eq_params/particle.txt.

3. Thermodynamics: defines the thermodynamical parameters. Its values can be
modified in ./eq_params/thermodynamics.txt.

4. Optical trap: defines the optical dipole trap parameters. Its values can be
modified in ./eq_params/optical_trap.txt.

5. Paul trap: defines the Paul trap parameters. Its values can be
modified in ./eq_params/paul_trap.txt.

Since the parameters are imported from files, the values can be modified
without the need to recompile the code every time. 

More details can be found by inspecting the struct definitions in 
physics.cpp/physics.hpp

> By: Gerard Planes Conangla









