/********************************************************************** 
 * DESCRIPTION:
 * 
 * Function + libraries to generate sample paths of a given stochastic 
 * process, defined by a user-defined SDE
 * 
 * dY = a(t,Y)dt + b(t,Y)*dW_t, Y(t0) = Y0,
 * 
 * using a Runge-Kutta type method for SDE of strong order 1.
 * This method doesn't require non-zero derivatives of b 
 * (diffusion term), since many methods (e.g. the Milstein method) 
 * have an effective strong order < 1 when b is a constant.
 * 
 ***********************************************************************  
 * INPUTS:  Call function with
 *          >> ./sde.out inp1 inp2 inp3 inp4 inp5 inp6 inp7...
 *          inp1: dt, time discretization step. Type: double
 *          inp2: t_ini, initial time. Type: double
 *          inp3: t_end, final time. Type: double
 *          inp4: num_traces, number of traces to simulate. Type: int
 *                Obs: needs to be num_traces%4 = 0
 *          inp5: n_dim, number of dimensions; e.g., harmonic oscillator
 *                would be 2. Type: int
 *          inp6: y0[0], 1st initial condition. Type: double
 *          inp7: y0[1], 2nd initial condition. Type: double. More inputs
 *                if n_dim is larger.
 *          E.g:
 *          ./sde.out 1e-5 0 1e-1 100 2 0 0
 * 
 * OUTPUTS: Prints elapsed time on stdout. Prints sample traces as rows 
 *          in "./simulated_traces/sde_sample_path_i.txt", will create 
 *          4 of them (one for every thread, see Observations). 
 *********************************************************************** 
 * OBSERVATIONS:
 * 
 * Uses 4 threads (3 + main), assuming a computer with 4 cores.
 * 
 * Compile with:
 * >> g++ sde.cpp sde_library.cpp -o sde.out -O3 -pthread
 * This assumes that sde.cpp, sde_library.cpp and sde_library.h are on 
 * the same folder.
 *********************************************************************** 
 * REFERENCES:
 * 
 * GIT REPOSITORY
 * https://github.com/gerardpc/sde_simulator
 * 
 * SDE RK method
 * https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_method_%28SDE%29#Variation_of_the_Improved_Euler_is_flexible
 ***********************************************************************  
 * Versions: 
 *  By GP Conangla
 *  01.10.2019
 *      Obs: Working function. Prints on a file estimated <x^2(t)>
 *      Performance, compared with pure MATLAB code is about x500 times 
 *      faster. Using main + 3 more threads (i.e., assumes computer with
 *      4 cores).
 *********************************************************************** 
 */

#include <iostream> // stdin, stdout 
#include <cstdio> // printf family
#include <cmath> // most math functions
#include <vector> // vector library
#include <string> // string library
#include <random> // random number generation
#include <stdlib.h> // includes rand()
#include <chrono> // time related library, for seeeding
#include <thread> // for multithreading

// C++ function that mimics MATLAB linspace
std::vector<double> linspace(std::vector<double> t_interval, int n);

// C++ function that mimics MATLAB a:h:b
std::vector<double> generate_v_from_h(std::vector<double> t_interval, double dt);

// C++ function that mimics MATLAB zeros (but with arbitrary "value")
std::vector<std::vector<double>> ini_matrix(int n_arrays, int n_columns, double value);

// C++ function to print array on stdout
int print_array_asrow(std::vector<std::vector <double>> a, int dim, std::FILE* fp);

// C++ function to print vector on stdout
int print_vector(std::vector<double> v, std::FILE* fp);

// C++ function to print array on stdout
int print_array(std::vector<std::vector<double>> a, std::FILE* fp);

// Print "hola!" in stdout for debugging
int hola();

// sum of 2 vectors OF THE SAME SIZE
std::vector<double> vector_sum(std::vector<double> &a, std::vector<double> &b);

// sum of 2 arrays OF THE SAME SIZE
std::vector<std::vector<double>> array_sum(std::vector<std::vector<double>> &a,
    std::vector<std::vector<double>> &b);

// dot product of 2 vectors OF THE SAME SIZE
std::vector<double> dot_product(std::vector<double> &a, std::vector<double> &b);

// product of array times scalar
std::vector<std::vector<double>> array_scalar_multiplication(std::vector<std::vector<double>> &a, double k);

// dot product of 2 arrays OF THE SAME SIZE
std::vector<std::vector<double>> array_dot_product(std::vector<std::vector<double>> &a,
    std::vector<std::vector<double>> &b);

// function that is applied element by element in function_array
double f(double x);

// apply f(x) on array element by element
std::vector<std::vector<double>> function_array(std::vector<std::vector<double>> &a);

// product of vector times scalar
std::vector<double> scalar_multiplication(std::vector<double> &a, double k);

// drift of process
std::vector<double> drift_function(std::vector<double> y, double t);

// diffusion of process
std::vector<double> diffusion_function(std::vector<double> y, double t);

// Roll a dice between low and high
int roll_dice(int low, int high);

// Runge-Kutta function as programmed in MATLAB
void runge_kutta(std::vector<double> t_interval,
    std::vector<double> y0, double dt, std::vector<std::vector<double>> &y);

// Generate num_traces traces with the RK method and print 
// each of them to file in different rows
void generate_traces(int num_traces, std::string filename, 
    std::vector<double> t_interval, std::vector<double> y0, double dt);
    
// Generate num_traces traces with the RK method and return 
// sum of variance to file. Obs: you need to divide by num_traces
// outside this function to get average!
void generate_avg_trace(int num_traces, std::vector<double> t_interval, 
std::vector<double> y0, std::vector<std::vector<double>> &avg_var, double dt);
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
