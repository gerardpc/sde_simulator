/********************************************************************** 
 * DESCRIPTION:
 * 
 * Library with functions to simulate a dynamical system defined by an
 * arbitrary (user-defined) SDE
 * 
 * dY = a(t,Y)dt + b(t,Y)*dW_t, Y(t0) = Y0,
 * 
 * Two simulation options are allowed: ODE type (deterministic and,
 * therefore, only drift and no diffusion) and SDE type (stochastic, 
 * considering both terms). The numerical methods used differ depending 
 * on the case: 4th order Adams predictor-corrector for ODEs, Runge-Kutta
 * type method of strong order 1 and deterministic order 2 for SDEs.
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
 *  14.10.2019
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
#include "eq_definitions.hpp"
#include "num_vector.hpp" // basic vector manipulation library

//======================================================================
// LIBRARY FUNCTIONS
//======================================================================
// DETERMINISTIC
// 4th order Runge-Kutta
// Used to calculate first steps in Adams predictor-corrector method
void runge4(std::vector<double> t_interval, std::vector<double> y0, double dt,
std::vector<std::vector<double>> &y, const eq_params &args);

// 4-step Adams predictor-corrector
void adams_pc(std::vector<double> t_interval, std::vector<double> y0, double dt,
std::vector<std::vector<double>> &y, const eq_params &args);

int adams_all(std::vector<double> t_interval, std::vector<double> y0, 
std::vector<std::vector<double>> &y, double dt, const eq_params &args);

//======================================================================
// STOCHASTIC
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
std::vector<double> y0, std::vector<std::vector<double>> &fun_avg, 
std::vector<std::vector<double>> &fun_var, double dt, bool Ito, 
const eq_params &args, unsigned int thread_id);
    
// Generate num_traces traces with generate_avg_trace in 1 or 4 different
// threads. Print elapsed time for execution on stdout. Returns average
// of f(Y_t)
int RK_all(unsigned int num_traces, bool many_traces, 
std::vector<double> t_interval, std::vector<double> y0, 
std::vector<std::vector<double>> &fun_avg, std::vector<std::vector<double>> &fun_var,
double dt, bool Ito, const eq_params &args);

// Print results   
// avg trace number i (where i is degree of freedom number i) will be
// printed on file ./simulated_traces/sde_sample_path_i.txt
int print_results(unsigned int n_dim, const std::vector<std::vector<double>> fun_avg,
unsigned int subsampling_f, std::string statistic);

// fill problem parameters with inputs, if given (otherwise use default
// values)
int fill_parameters_w_inputs(int argc, char* argv[], double &dt, 
std::vector<double> &t_interval, unsigned int &num_traces, bool &many_traces, 
unsigned int &subsampling_f, std::string &eq_type, bool &Ito, unsigned int &n_dim, 
std::vector<double> &y0);
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
