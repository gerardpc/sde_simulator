# SDE_SIMULATOR project

Description
--------------
Functions + libraries to generate sample paths of a given vectorial
stochastic process, defined by the SDE

dY = a(t,Y)dt + b(t,Y)*dW_t, Y(t0) = Y0

and where Y has arbitrary (but finite!) dimension.

A typical application would be generating n sample paths Y_t, computing
a certain function of them f(Y_t) and then averaging the result, to find 
the expected value of f(Y_t). The most common case is calculating the 2nd 
moment, the variance of Y(t): <Y^2(t)>.

The code uses a stochastic Runge-Kutta type method of strong order 1
that does not require any non-zero derivatives of b (diffusion term).
Other methods (e.g. the Milstein method) have strong order 1 but reduce
to the Euler-Maruyama method (strong order 0.5) when b is a constant.

The code is divided in two parts:

- A C++ main() + library functions for the most numerically intensive 
part: generating the traces Y_t, calculating f(Y_t) and averaging them. 
Then write the result in a file. Note that in the present form, the C++
main calculates the VARIANCE, i.e. f(Y_t) = Y_t^2, and not any arbitrary
function. The code uses multithreading and assumes a computer with 4 cores.

- MATLAB and Python scripts to load the saved trace and plot the results.

In my personal computer, this code boosts the performance (i.e., reduces
execution time) a factor 500 with respect to equivalent code implemented
100% in MATLAB.

Everything is well commented, so you should be able to modify the code 
(whether the C++ or the MATLAB/Python part) easily and perform your own 
numerical analysis.

The RK method is an implementation of 
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

    >> g++ sde.cpp sde_library.cpp -o sde.out -O3 -pthread


2- Generate traces and plot the results
--------------
If you are working with MATLAB:

Open MATLAB, go to the MATLAB folder. Call the function "numerical_sde_cpp".
If you inspect the code, all the inputs and outputs are detailed. E.g.

    >> [tt, avg_var] = numerical_sde_cpp(1e-5, [0 1], 100, 2, [0 0]);

The results will be plotted in different figures.

If you are working with Python:

Open a terminal, go to the Python folder. Call the function "sde_script.py".
If you inspect the code, all the parameters are detailed and can be changed.
E.g.

    >> python sde_script.py

The results will be saved as png figures. Check the "test.png" file for
an example output.


2- Changing the a(t,Y) and b(t,Y) functions of the SDE
--------------
The a(t,Y) and b(t,Y) functions are defined in the file sde_library.cpp 
as the drift_function and diffusion_function, respectively.

The f[i] vector inside each of these functions contains the driving
term for the equations of motion (the drift_function is the "deterministic"
part and the diffusion_function the "stochastic" part). 

For instance, consider a Brownian particle completely driven by thermal 
noise (i.e., not subject to any force potential). If the particle only 
moves in one dimension, the equations will be

 x' = v
 v' = 0 + sigma/m

The first equation is clear, and the 2nd equation is Newton's 2nd law.
This is precisely the equation included in the code by default, and 
should be used as an example for other, arbitrary, drift and diffusion
functions. Of course, every time any of these functions is changed, the
C++ code should be recompiled with the command indicated in section 1.


> By: Gerard Planes Conangla









