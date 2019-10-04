/********************************************************************** 
 * DESCRIPTION:
 * 
 * Function to generate sample paths of a given stochastic 
 * process, defined by the SDE
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
 * OUTPUTS: Prints elapsed time and other details on stdout. 
 *          Prints sample traces as rows in files
 *          "./simulated_traces/sde_sample_path_i.txt"
 *          Every file contains one of the degrees of freedom
 *********************************************************************** 
 * OBSERVATIONS:
 * 
 * Detects number of cores and divides traces generation among threads.
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
 *  04.10.2019
 *      Obs: Working function. Prints on a file estimated <x^2(t)>
 *      Performance, compared with pure MATLAB code is about x500 times 
 *      faster with 4 cores. Automatically adapts to number of cores.
 *********************************************************************** 
 */

// STANDARD LIBRARIES
#include <vector> // vector library

// MY LIBRARIES
#include "sde_library.h"

using namespace std; // avoid writing std::function every time

// main function: read problem inputs, generate traces Y_i, calculate
// <f(Y_i)> and print results in files
int main(int argc, char* argv[]){    
    //==================================================================
    // Declare problem parameters
    // dt: discretization step
    double dt;
    // time interval
    vector<double> t_interval;
    // number of traces to simulate
    int num_traces;
    bool many_traces; // true if num_traces >= 4
    // problem dimension (e.g., for a harm. oscillator = 2) 
    int n_dim;
    // initial conditions, if given
    vector<double> y0;
    
    // Fill problem parameters with inputs, if given.
    // Otherwise use default values: dt = 1e-5, t_interval = [0 1e-1],
    // num_traces = 100
    fill_parameters_w_inputs(argc, argv, dt, t_interval, num_traces,
    many_traces, n_dim, y0);

    //==================================================================
    // RK method. Put <f(Y_t)> on avg_trace
    vector<vector<double>> avg_trace = RK_all(num_traces, many_traces, t_interval, y0, dt);
    
    //==================================================================    
    // Print results   
    // avg trace number i (where i is degree of freedom number i) will be
    // printed on file ./simulated_traces/sde_sample_path_i.txt
    print_results(n_dim, avg_trace);
    
    //==================================================================    
    // Finish execution    
    return 0;
}







































