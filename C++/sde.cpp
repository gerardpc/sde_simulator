/********************************************************************** 
 * DESCRIPTION:
 * 
 * Function to numerically simulate a dynamical system defined by the SDE
 * 
 * dY = a(t,Y)dt + b(t,Y)*dW_t, Y(t0) = Y0.
 * 
 * Two simulation options are allowed:
 *
 *      - Deterministic dynamical system (ODE). No noise in the equations of 
 *        motion, an Adams predictor-corrector method of 4th order is used.
 *
 *      - Stochastic dynamical system (SDE): a Runge-Kutta type method for SDE
 *        of strong order 1 and deterministic order 2 is used. The method
 *        maintains the order even when the diffusion term is a constant
 *        (as opposed to, for instance, Milstein method). 
 *        Since every run of the method for an SDE will be a different 
 *        sample path, the code is optimized to generate many sample paths
 *        at once using multithreading. When the simulation phase is finished,
 *        the averages and standard errors are calculated.
 * 
 ***********************************************************************  
 * INPUTS:  Call function with
 *          >> ./sde.out inp1 inp2 inp3 inp4 inp5 inp6 inp7...
 *          inp1:  dt, time discretization step. Type: double
 *          inp2:  t_ini, initial time. Type: double
 *          inp3:  t_end, final time. Type: double
 *          inp4:  num_traces, number of traces to simulate. Type: int
 *                 Obs: needs to be num_traces%num_cores = 0 (otherwise 
 *                 program will round to last %num_cores = 0 number). This
 *                 parameter is not used in the ODE case: only one trace
 *                 is generated.
 *          inp5:  subsampling_f: even though N points will be generated
 *                 per trace, print only N/10 points. Helps when small dt 
 *                 is required for numerical stability, but the equations 
 *                 can be well represented with less points (i.e., big
 *                 speed-up in importing data to MATLAB/Python). Type: int
 *          inp6:  eq_type: selects ode or sde equation type. Type: string
 *          inp7:  Ito: selects Ito or Stratonovich equation. Type: bool.
 *                 This parameter is not used in the ODE case.
 *          inp8:  n_dim, number of dimensions; e.g., harmonic oscillator
 *                 would be 2. Type: int
 *          inp9:  y0[0], 1st initial condition. Type: double
 *          inp10: y0[1], 2nd initial condition. Type: double. 
 *                 More inputs for y0 if n_dim is larger.
 *          E.g:
 *          ./sde.out 1e-5 0 1e-1 100 1 2 0 0
 * 
 * OUTPUTS: Prints elapsed time and other details on stdout. 
 *          Prints sample traces as rows in files
 *          "./simulated_traces/sde_sample_path_i.txt"
 *          Every file contains one of the degrees of freedom
 *********************************************************************** 
 * OBSERVATIONS:
 * 
 * Detects number of cores and divides traces generation among threads
 * (only in the case of SDE simulation; for ODEs it doesn't make any sense
 * since all of them are going to be the same).
 * 
 * Compile with:
 * >> g++ sde.cpp sde_library.cpp num_vector.cpp particle_traps.cpp -o sde.out -O3 -pthread
 *
 * This assumes that sde.cpp, sde_library.cpp and sde_library.h are on 
 * the same folder.
 *********************************************************************** 
 * REFERENCES:
 * 
 * GIT REPOSITORY
 * https://github.com/gerardpc/sde_simulator
 * 
 * ODE Adams method
 * https://en.wikiversity.org/wiki/Adams-Bashforth_and_Adams-Moulton_methods
 * 
 * SDE RK method
 * https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_method_(SDE)
 ***********************************************************************  
 * Versions: 
 *  By GP Conangla
 *  14.10.2019
 *      Obs: Working function. Prints on a file estimated <x^2(t)>
 *      Performance, compared with pure MATLAB code is about x500 times 
 *      faster with 4 cores. Automatically adapts to number of cores.
 *********************************************************************** 
 */

// STANDARD LIBRARIES
#include <iostream> // stdin, stdout 
#include <vector> // vector library
#include <string> // string library
#include <stdlib.h> // includes "exit"

// MY LIBRARIES
#include "num_vector.hpp"
#include "sde_library.hpp"

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
    unsigned int num_traces;
    bool many_traces; // true if num_traces >= 4
    // Subsampling factor
    unsigned int subs_f;
    // Equation type: ode or sde?
    string eq_type;
    // Ito (true) or Stratonovich (false)?
    bool Ito;
    // problem dimension (e.g., for a harm. oscillator = 2) 
    unsigned int n_dim;
    // initial conditions, if given
    vector<double> y0;
    
    // Fill problem parameters with inputs, if given.
    // Otherwise use default values (see sde_library::fill_parameters)
    fill_parameters_w_inputs(argc, argv, dt, t_interval, num_traces, many_traces, subs_f, eq_type, Ito, n_dim, y0);
    
    //==================================================================
    // import trap equations from files in eq_params folder
    // (only necessary if particle trap equations are going to be used)
    eq_params trap;
    trap.fill();
    //trap.print();
    
    // Define average trace vectors
    vector<vector<double>> fun_avg;
    vector<vector<double>> fun_var;
    
    // Call numerical method
    if(eq_type == "ode"){ // Deterministic case    
        // call Adams method    
        adams_all(t_interval, y0, fun_avg, dt, trap);
        fun_var = ini_matrix(fun_avg.size(), y0.size(), 0);
    } else if(eq_type == "sde"){ // Stochastic case
        // call RK method. Put <f(Y_t)> on fun_avg and <f^2(Y_t)> on fun_var    
        RK_all(num_traces, many_traces, t_interval, y0, fun_avg, fun_var, dt, Ito, trap);
    } else { // unknown equation type
        std::cout << "Unknown equation type. Please, specify one of the following: 'ode', 'sde'.\n";
        return 1;
    }
    
    //==================================================================    
    // Print results   
    // avg trace number i (where i is degree of freedom number i) will be
    // printed on file ./simulated_traces/sde_sample_path_i.txt
    std::cout << "-> Average of f(x(t))\n";
    print_results(n_dim, fun_avg, subs_f, "avg");    
    std::cout << "\n-> Variance of f(x(t))\n";
    print_results(n_dim, fun_var, subs_f, "var");
    
    //==================================================================    
    // Finish execution    
    return 0;
}

































