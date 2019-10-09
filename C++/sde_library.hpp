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
 *  09.10.2019
 *      Obs: Working library. Prints on a file estimated <x^2(t)>. This
 *      function can be changed (defined in num_vector.hpp as
 *      "double f(double x)".
 *********************************************************************** 
 */
 
// STANDARD LIBRARIES
#include <vector> // vector library
#include <string> // string library
#include <stdlib.h> // includes rand()
// MY LIBRARIES
#include "particle_traps.hpp" // optical tweezer, paul trap and gaussian beam

//======================================================================
// EQUATION PARAMETERS STRUCT
//======================================================================
struct eq_params{
    // the members of this struct are specific to my problem of interest,
    // but can be replaced to any other members/parameters to send variables
    // inside the drift and diffusion functions
    opt_tweezer ot;
    paul_trap pt;
    gaussian_beam gb;
    void fill();
};

//======================================================================
// PROBLEM FUNCTIONS
//======================================================================
// SDE FUNCTIONS
// dY = a(t,Y)dt + b(t,Y)*dW_t, Y(t0) = Y0,

// drift of process: a(t,Y)
std::vector<double> drift_function(std::vector<double> y, double t, const eq_params &args);

// diffusion of process: b(t, Y)
std::vector<double> diffusion_function(std::vector<double> y, double t, const eq_params &args);

//======================================================================
// LIBRARY FUNCTIONS
//======================================================================

// Runge-Kutta function as programmed in MATLAB
void runge_kutta(std::vector<double> t_interval,
    std::vector<double> y0, double dt, std::vector<std::vector<double>> &y, 
    bool Ito, const eq_params &args);

// Generate num_traces traces with the RK method and print 
// each of them to file in different rows
void generate_traces(unsigned int num_traces, std::string filename, 
std::vector<double> t_interval, std::vector<double> y0, double dt, 
unsigned int subsampling_f, bool Ito, const eq_params &args);
    
// Generate num_traces traces with the RK method and return 
// sum of variance to file. Obs: you need to divide by num_traces
// outside this function to get average!
void generate_avg_trace(unsigned int num_traces, std::vector<double> t_interval, 
std::vector<double> y0, std::vector<std::vector<double>> &avg_var, double dt, 
bool Ito, const eq_params &args);
    
// Generate num_traces traces with generate_avg_trace in 1 or 4 different
// threads. Print elapsed time for execution on stdout. Returns average
// of f(Y_t)
std::vector<std::vector<double>> RK_all(unsigned int num_traces, bool many_traces, 
std::vector<double> t_interval, std::vector<double> y0, double dt, bool Ito, const eq_params &args);

// Print results   
// avg trace number i (where i is degree of freedom number i) will be
// printed on file ./simulated_traces/sde_sample_path_i.txt
int print_results(unsigned int n_dim, const std::vector<std::vector<double>> avg_trace,
unsigned int subsampling_f);

// fill problem parameters with inputs, if given (otherwise use default
// values)
int fill_parameters_w_inputs(int argc, char* argv[], double &dt, 
std::vector<double> &t_interval, unsigned int &num_traces, bool &many_traces, 
unsigned int &subsampling_f, bool &Ito, unsigned int &n_dim, std::vector<double> &y0);
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
