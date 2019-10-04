/********************************************************************** 
 * DESCRIPTION:
 * 
 * Library to generate sample paths of a given stochastic 
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
 * OBSERVATIONS:
 * 
 * Function RK_all (arguably the main function of the library) 
 * automatically adapts to number of cores (divides trace generation among
 * them).
 *
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
 *  04.10.2019
 *      Obs: Working library. Prints on a file estimated <x^2(t)>. This
 *      function can be changed, defined as function double f(double x).
 *      Performance, compared with pure MATLAB code is about x500 times 
 *      faster with 4 cores. Automatically adapts to number of cores.
 *********************************************************************** 
 */
 
// STANDARD LIBRARIES
#include <vector> // vector library
#include <string> // string library
#include <stdlib.h> // includes rand()


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
    
// Generate num_traces traces with generate_avg_trace in 1 or 4 different
// threads. Print elapsed time for execution on stdout. Returns average
// of f(Y_t)
std::vector<std::vector<double>> RK_all(int num_traces, bool many_traces, 
std::vector<double> t_interval, std::vector<double> y0, double dt);

// Print results   
// avg trace number i (where i is degree of freedom number i) will be
// printed on file ./simulated_traces/sde_sample_path_i.txt
int print_results(int n_dim, std::vector<std::vector<double>> avg_trace);

// fill problem parameters with inputs, if given (otherwise use default
// values)
int fill_parameters_w_inputs(int argc, char* argv[], double &dt, 
std::vector<double> &t_interval, int &num_traces, bool &many_traces, 
int &n_dim, std::vector<double> &y0);
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
