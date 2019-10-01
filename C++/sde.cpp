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

// STANDARD LIBRARIES
#include <iostream> // stdin, stdout 
#include <cstdio> // printf family
#include <cmath> // most math functions
#include <vector> // vector library
#include <string> // string library
#include <random> // random number generation
#include <stdlib.h> // includes rand()
#include <chrono> // time related library, for seeeding
#include <thread> // for multithreading

// MY LIBRARIES
#include "sde_library.h"

using namespace std; // avoid writing std::function every time

int main(int argc, char* argv[]){    
    //==================================================================
    // Define/declare parameters (use default or inputs, if any)
    // dt: discretization step
    double dt;
    if(argc < 2){
        dt = 1e-5; // default value
    } else {
        dt = atof(argv[1]); // input value
    }
    // time interval
    std::vector<double> t_interval;
    if(argc < 3){ 
        t_interval = {0, 1e-1}; // default value
    } else {
        t_interval = {atof(argv[2]), atof(argv[3])}; // input value
    }
    // number of traces to simulate
    int num_traces;
    bool many_traces;
    if(argc < 5){ 
        num_traces = 25; // default value, 25*4 = 100 traces
    } else {
        // case of few traces: not multithreading
        if(atoi(argv[4]) < 4){
            num_traces = atoi(argv[4]);
            many_traces = false;
            // print number of traces on stdout
            printf("Number of traces that will be generated: %d\n", num_traces);
        // case of many traces
        } else {
            // input value of num_traces
            // notice it gets rounded to 0 mod(4)
            num_traces = atoi(argv[4])/4;             
            many_traces = true;
            // print number of traces on stdout
            printf("Number of traces that will be generated: %d\n", num_traces*4);
        }
    }
    // problem dimension
    int n_dim;
    if(argc < 6){ 
        n_dim = 2; // 2D problem, think of x and v.
    } else {
        n_dim = atoi(argv[5]); // input value
    }
    // initial conditions, if given
    std::vector<double> y0(n_dim);
    if(argc < 6){ 
        // default: initial conditions = 0
        y0 = {0, 0};
    } else {
        // read initial conditions
        for(int i = 0; i < n_dim; i++){
            y0[i] = atof(argv[6 + i]);
        }
    }

    //==================================================================
    // RK method    
    // measure initial time
    auto start = std::chrono::high_resolution_clock::now(); 
    
    // preallocate average of f(x) trace vectors
    std::vector<std::vector<double>> avg_var1;
    std::vector<std::vector<double>> avg_var2;
    std::vector<std::vector<double>> avg_var3;
    std::vector<std::vector<double>> avg_var4; 
    
    // Run 4 parallel RK methods to speed code x4 
    // run 3 simulations in thread t1, t2, t3
    if(many_traces){
        // multithreading if num_traces >= 4        
        thread t1(generate_avg_trace, num_traces, t_interval, y0, ref(avg_var1), dt);
        thread t2(generate_avg_trace, num_traces, t_interval, y0, ref(avg_var2), dt);
        thread t3(generate_avg_trace, num_traces, t_interval, y0, ref(avg_var3), dt);
    
        // run 1 more simulation in main
        generate_avg_trace(num_traces, t_interval, y0, avg_var4, dt);  
    
        // join threads when they finish
        t1.join();
        t2.join();
        t3.join();
    
        // Calculate sum of variances and put it to avg_var1
        avg_var1 = array_sum(avg_var1, avg_var2);
        avg_var1 = array_sum(avg_var1, avg_var3);
        avg_var1 = array_sum(avg_var1, avg_var4);
    
        // divide by num_traces to estimate <f(y_t)_i(t)>
        double num_tracesf = num_traces;
        double ntraces_factor = 1/(num_tracesf*4);
        avg_var1 = array_scalar_multiplication(avg_var1, ntraces_factor);
        
        // calculate and print execution time
        auto stop = std::chrono::high_resolution_clock::now(); 
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count();
        printf("Execution time to simulate %d x %g s trace with dt = %g: %g s\n", 
                num_traces*4, t_interval[1] - t_interval[0], dt, ((float) duration)/1e6);
    } else {
        // generate single thread
        generate_avg_trace(num_traces, t_interval, y0, avg_var1, dt);  
        // divide by num_traces to estimate <f(y_t)_i(t)>
        double num_tracesf = num_traces;
        double ntraces_factor = 1/(num_tracesf);
        avg_var1 = array_scalar_multiplication(avg_var1, ntraces_factor);
        
        // calculate and print execution time
        auto stop = std::chrono::high_resolution_clock::now(); 
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count();
        printf("Execution time to simulate %d x %g s trace with dt = %g: %g s\n", 
                num_traces, t_interval[1] - t_interval[0], dt, ((float) duration)/1e6);
    }
    
    //==================================================================    
    // Print results   
    // filenames vector of strings.
    // avg trace number i (where i is degree of freedom number i) will be
    // printed on file filename[i]
    vector<string> filename(n_dim);
    cout << "Average traces saved in files:\n";
    for(int i = 0; i < n_dim; i++){
        // generate new filename
        filename[i] = "./simulated_traces/sde_sample_path_" + to_string(i) + ".txt";
        cout << filename[i] << "\n";
    }
    
    // Print results of dimension dim to file
    // Notice that, by default, this APPENDS a new row to the file
    for(int i = 0; i < n_dim; i++){
        // open file
        FILE* fp = fopen(filename[i].c_str(), "a"); // <-- a for "append"
        if (fp == NULL) {
          fprintf(stderr, "Can't open output file\n");
          exit(1);
        }
        // print to file
        print_array_asrow(avg_var1, i, fp);    
        // Close file
        fclose(fp);
    }
    
    //==================================================================    
    // Finish execution    
    return 0;
}







































