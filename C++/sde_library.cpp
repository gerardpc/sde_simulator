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
#include <iostream> // stdin, stdout 
#include <cstdio> // printf family
#include <cmath> // most math functions
#include <vector> // vector library
#include <string> // string library
#include <random> // random number generation
#include <chrono> // time related library, for seeeding
#include <thread> // for multithreading

// MY LIBRARIES
#include "num_vector.hpp" // basic vector manipulation library
#include "sde_library.hpp"


//======================================================================
// LIBRARY FUNCTIONS
//======================================================================
// DETERMINISITC
// 4th order Runge-Kutta
// Used to calculate first steps in Adams predictor-corrector method
void runge4(std::vector<double> t_interval, std::vector<double> y0, double dt,
std::vector<std::vector<double>> &y, const eq_params &args){

    // PARAMETERS AND INITIALIZATION
    // generate time vector    
    std::vector<double> t = generate_v_from_h(t_interval, dt);
    unsigned int n_time = t.size();

    // initialize y with initial conditions as y(0,:) = y0
    y[0] = y0;
    
    // initialize temp vectors
    std::vector<double> k1(y0.size());
    std::vector<double> k2(y0.size());
    std::vector<double> k3(y0.size());
    std::vector<double> k4(y0.size());
    std::vector<double> tmp1(y0.size());
    std::vector<double> tmp2(y0.size());
    std::vector<double> tmp3(y0.size());
    
    // RK method
    for(unsigned int n = 1; n < n_time; n++){
        // k1
        k1 = drift_function(y[n - 1], t[n - 1], args);
        
        // k2
        tmp1 = scalar_multiplication(k1, dt/2);
        tmp1 = vector_sum(y[n - 1], tmp1);
        k2 = drift_function(tmp1, t[n - 1] + dt/2, args);
        
        // k3
        tmp1 = scalar_multiplication(k2, dt/2);
        tmp1 = vector_sum(y[n - 1], tmp1);
        k3 = drift_function(tmp1, t[n - 1] + dt/2, args);
        
        // k4
        tmp1 = scalar_multiplication(k3, dt);
        tmp1 = vector_sum(y[n - 1], tmp1);
        k4 = drift_function(y[n - 1], t[n - 1] + dt, args);
        
        // next y value
        tmp1 = vector_sum(k1, k4);
        tmp2 = vector_sum(k2, k3);
        tmp2 = scalar_multiplication(tmp2, 2);
        tmp3 = vector_sum(tmp1, tmp2);
        tmp3 = scalar_multiplication(tmp3, dt/6);
        y[n] = vector_sum(y[n - 1], tmp3);
    }
    // solution is in y. End
}

// 4-step Adams predictor-corrector
void adams_pc(std::vector<double> t_interval, std::vector<double> y0, double dt,
std::vector<std::vector<double>> &y, const eq_params &args){

    // PARAMETERS AND INITIALIZATION
    // generate time vector    
    std::vector<double> t = generate_v_from_h(t_interval, dt);
    unsigned int n_time = t.size();
    
    // preallocate solution vector y, full of zeros
    y = ini_matrix(n_time, y0.size(), 0);

    // use Runge-Kutta method to get 4 initial values in y
    std::vector<double> t_interval_rk{t[0], t[3]};
    runge4(t_interval_rk, y0, dt, y, args);
    
    // Temporary vectors
    std::vector<double> tmp1(y0.size());
    std::vector<double> tmp2(y0.size());
    std::vector<double> tmp3(y0.size());
    std::vector<double> tmp4(y0.size());
    
    // The Adams method itself
    for(unsigned int n = 4; n < n_time; n++){
        // Predictor: this y[n] is the predicted value
        // previous steps
        tmp1 = scalar_multiplication(drift_function(y[n - 1], t[n - 1], args), 55);
        tmp2 = scalar_multiplication(drift_function(y[n - 2], t[n - 2], args), -59);
        tmp3 = scalar_multiplication(drift_function(y[n - 3], t[n - 3], args), 37);
        tmp4 = scalar_multiplication(drift_function(y[n - 4], t[n - 4], args), -9);
        // combine steps
        tmp1 = vector_sum(tmp1, tmp2);
        tmp1 = vector_sum(tmp1, tmp3);
        tmp1 = vector_sum(tmp1, tmp4);
        tmp1 = scalar_multiplication(tmp1, dt/24);
        y[n] = vector_sum(y[n - 1], tmp1);
            
        // Corrector: now the value y[n] is updated with the corrected value
        tmp1 = scalar_multiplication(drift_function(y[n], t[n], args), 9);
        tmp2 = scalar_multiplication(drift_function(y[n - 1], t[n - 1], args), 19);
        tmp3 = scalar_multiplication(drift_function(y[n - 2], t[n - 2], args), -5);
        tmp4 = scalar_multiplication(drift_function(y[n - 3], t[n - 3], args), 1);
        // combine steps
        tmp1 = vector_sum(tmp1, tmp2);
        tmp1 = vector_sum(tmp1, tmp3);
        tmp1 = vector_sum(tmp1, tmp4);
        tmp1 = scalar_multiplication(tmp1, dt/24);
        y[n] = vector_sum(y[n - 1], tmp1);
    }
    // solution is in y. End
}

// Generate 1 trace with the Adams method (ODE case)
int adams_all(std::vector<double> t_interval, std::vector<double> y0, 
std::vector<std::vector<double>> &y, double dt, const eq_params &args){
    // measure initial time
    auto start = std::chrono::high_resolution_clock::now();     
    printf("------------------------------------------\n");
    printf("A single (deterministic) trace will be generated.\n");
    
    // generate single thread
    adams_pc(t_interval, y0, dt, y, args);
    
    // calculate <f(y_t)>
    y = function_array(y);
    
    // calculate and print execution time
    auto stop = std::chrono::high_resolution_clock::now(); 
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count();
    printf("\nExecution time to simulate one x %g s trace with dt = %.1e: %g s\n", 
            t_interval[1] - t_interval[0], dt, ((float) duration)/1e6);
    return 0;
}

//======================================================================
// STOCHASTIC
// Stochastic Runge-Kutta function: numerical code.
// output in y
void runge_kutta(std::vector<double> t_interval,
std::vector<double> y0, double dt, std::vector<std::vector<double>> &y, 
bool Ito, const eq_params &args){
    
    // PARAMETERS AND INITIALIZATION
    // generate time vector    
    std::vector<double> t = generate_v_from_h(t_interval, dt);
    unsigned int n_time = t.size();
    // preallocate solution vector y, full of zeros
    y = ini_matrix(n_time, y0.size(), 0);

    // initialize y with initial conditions as y(0,:) = y0
    y[0] = y0;
    
    srand (time(NULL)); // initialize random seed:
    double S_k;
    
    // k_1, k_2 for the RK method
    std::vector<double> k_1(y0.size());
    std::vector<double> k_2(y0.size());
    std::vector<double> tmp1(y0.size());
    std::vector<double> tmp2(y0.size());
    std::vector<double> tmp3(y0.size());
    
    // dW_t, Wiener process differential
    std::vector<double> dW_t(y0.size());
    // random engine generator
    unsigned seed1 = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed1);
    // normal distribution type N(0,1)
    std::normal_distribution<double> n_distribution(0, 1);
    
    // the RUNGE - KUTTA method itself
    for(unsigned int n = 1; n < n_time; n++){
        // calculate dW_t
        for(unsigned int i = 0; i < dW_t.size(); i++){
            dW_t[i] = n_distribution(generator)*sqrt(dt);
        }
        y[n] = vector_sum(y[n - 1], tmp1);
        // calculate S_k
        if(Ito){ // Ito case
            if(n_distribution(generator) >= 0){
                S_k = 1;
            } else {
                S_k = -1;
            }
        } else { // Stratonovich case
            S_k = 0;
        }
        
        // Calculate k_1
        tmp1 = drift_function(y[n - 1], t[n - 1], args);
        tmp1 = scalar_multiplication(tmp1, dt);

        for(unsigned int i = 0; i < dW_t.size(); i++){
            tmp2[i] = dW_t[i] - S_k*sqrt(dt);
        }

        tmp3 = diffusion_function(y[n - 1], t[n - 1], args);
        tmp2 = dot_product(tmp2, tmp3);

        k_1 = vector_sum(tmp1, tmp2);

        // Calculate k_2
        tmp1 = drift_function(vector_sum(y[n - 1], k_1), t[n], args);
        tmp1 = scalar_multiplication(tmp1, dt);
        
        for(unsigned int i = 0; i < dW_t.size(); i++){
            tmp2[i] = dW_t[i] + S_k*sqrt(dt);
        }
        
        tmp3 = diffusion_function(y[n - 1], t[n], args);      
        tmp2 = dot_product(tmp2, tmp3);  
        
        k_2 = vector_sum(tmp1, tmp2);        
        
        // update solution matrix
        tmp1 = vector_sum(k_1, k_2);
        tmp1 = scalar_multiplication(tmp1, 0.5);
        y[n] = vector_sum(y[n - 1], tmp1);
    }
    // solution is in y. End
}

// Generate num_traces traces with the RK method and print 
// each of them to file in different rows
void generate_traces(unsigned int num_traces, std::string filename, 
std::vector<double> t_interval, std::vector<double> y0, double dt, 
unsigned int subsampling_f, bool Ito, const eq_params &args){
    // open to write to file
    FILE* fp = fopen(filename.c_str(), "a");
    if (fp == NULL) {
      fprintf(stderr, "Can't open output file\n");
      exit(1);
    }
    // preallocate solution vector y
    std::vector<std::vector<double>> y;
    // call RK method num_traces times
    for(unsigned int i = 0; i < num_traces; i++){
        runge_kutta(t_interval, y0, dt, y, Ito, args);  
        // Print results to file 
        print_array_asrow(y, 0, subsampling_f, fp);
    }
    // Close file and exit
    fclose(fp);
}

// Generate num_traces traces with the RK method and return 
// sum of function <f(Y_t)> (Obs: remains to divide by num_traces
// to estimate average)
void generate_avg_trace(unsigned int num_traces, std::vector<double> t_interval, 
std::vector<double> y0, std::vector<std::vector<double>> &fun_avg, 
std::vector<std::vector<double>> &fun_var, double dt, bool Ito, 
const eq_params &args, unsigned int thread_id){
    // preallocate solution vector y and tmp1 vector
    std::vector<std::vector<double>> y;
    std::vector<std::vector<double>> tmp1;
    
    // call RK method num_traces times
    runge_kutta(t_interval, y0, dt, y, Ito, args); 
    fun_avg = function_array(y);
    fun_var = array_dot_product(fun_avg, fun_avg);
    if(thread_id == 0){ // main thread, print progress
        unsigned int progress;
        unsigned int last_progress = 0;
        // Print progress bar
        std::cout << "\nProgress:\n";
        std::cout << "0 " << std::flush;
        // Calculate rest of traces
        for(unsigned int i = 1; i < num_traces; i++){
            runge_kutta(t_interval, y0, dt, y, Ito, args); 
            // calculate <f(y_t)>
            tmp1 = function_array(y);
            fun_avg = array_sum(fun_avg, tmp1);
            // calculate Var(f(y_t)) (for standard error)
            tmp1 = array_dot_product(tmp1, tmp1);
            fun_var = array_sum(fun_var, tmp1);            
            progress = ((i + 1)*50)/num_traces; 
            if(progress > last_progress){
                for(unsigned int diff = 0; diff < progress - last_progress; diff++){
                    std::cout << "#";
                }
                std::cout << " " << progress*2 << "%" << std::flush;
                if(progress*2 < 10){
                    std::cout << "\b\b\b";
                } else if(progress*2 < 100) {                        
                    std::cout << "\b\b\b\b";
                }
                last_progress = progress;
            }
        }
    } else { // not main thread, do not print
        for(unsigned int i = 1; i < num_traces; i++){
            runge_kutta(t_interval, y0, dt, y, Ito, args); 
            // calculate <f(y_t)>
            tmp1 = function_array(y);
            fun_avg = array_sum(fun_avg, tmp1);
            // calculate Var(f(y_t)) (for standard error)
            tmp1 = array_dot_product(tmp1, tmp1);
            fun_var = array_sum(fun_var, tmp1);        
        }
    }
}

// Generate num_traces traces with generate_avg_trace using different
// threads (automatic core detection. 
// Print elapsed time for execution on stdout. Returns average
// of f(Y_t)
int RK_all(unsigned int num_traces, bool many_traces, 
std::vector<double> t_interval, std::vector<double> y0, 
std::vector<std::vector<double>> &fun_avg, std::vector<std::vector<double>> &fun_var,
double dt, bool Ito, const eq_params &args){
    // measure initial time
    auto start = std::chrono::high_resolution_clock::now(); 
    
    // get number of cores
    unsigned num_cores = std::thread::hardware_concurrency();
    // in case it is not able to detect any core, set to 1
    if(num_cores == 0){
        num_cores = 1;
    }
    unsigned num_threads = num_cores - 1;
    printf("------------------------------------------\n");
    printf("Detected number of cores: %d\n", num_cores);
    
    if(Ito){
        printf("Equation of type: Ito\n");
    } else {
        printf("Equation of type: Stratonovich\n");
    }
    
    // Run parallel RK methods to speed code x num_threads 
    if(many_traces and num_threads >= 1){ // if multithreading
        // calculate traces per core
        num_traces = num_traces - num_traces % num_cores;
        printf("Number of traces that will be generated: %d\n", num_traces);
        int traces_per_core = num_traces/num_cores;
        printf("Number of traces per core: %d\n", traces_per_core);
        
        // vector of RK results and threads
        std::vector<std::vector<std::vector<double>>> thread_fun_avg(num_threads);
        std::vector<std::vector<std::vector<double>>> thread_fun_var(num_threads);
        std::vector<std::thread> t;
        
        // initiate simulations in threads
        for(unsigned int i = 0; i < num_threads; i++){
            t.push_back(std::thread(generate_avg_trace, traces_per_core, t_interval, y0, ref(thread_fun_avg[i]), 
            ref(thread_fun_var[i]), dt, Ito, args, i + 1));
        }
        
        // run 1 more simulation in main
        generate_avg_trace(traces_per_core, t_interval, y0, fun_avg, fun_var, dt, Ito, args, 0);  

        // join threads when they finish
        for(unsigned int i = 0; i < num_threads; i++){
            t[i].join();
        }
    
        // Calculate sum of averages and put it on fun_avg
        for(unsigned int i = 0; i < num_threads; i++){
            fun_avg = array_sum(thread_fun_avg[i], fun_avg);
        }
        // Calculate sum of variances and put it on fun_avg
        for(unsigned int i = 0; i < num_threads; i++){
            fun_var = array_sum(thread_fun_var[i], fun_var);
        }
        
        // divide by num_traces to estimate <f(y_t)_i>
        double num_tracesf = num_traces;
        double ntraces_factor = 1/(num_tracesf);
        fun_avg = array_scalar_multiplication(fun_avg, ntraces_factor);
        // divide by num_traces to estimate Var(f(y_t)_i)
        fun_var = array_scalar_multiplication(fun_var, ntraces_factor);
        
        // calculate and print execution time
        auto stop = std::chrono::high_resolution_clock::now(); 
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count();
        printf("\n\nExecution time to simulate %d x %g s trace with dt = %.1e: %g s\n\n", 
                num_traces, t_interval[1] - t_interval[0], dt, ((float) duration)/1e6);

    } else {  // not multithreading
        printf("Number of traces that will be generated: %d\n", num_traces);
        // generate single thread
        generate_avg_trace(num_traces, t_interval, y0, fun_avg, fun_var, dt, Ito, args, 0);  
        // divide by num_traces to estimate <f(y_t)_i(t)>
        double num_tracesf = num_traces;
        double ntraces_factor = 1/(num_tracesf);
        fun_avg = array_scalar_multiplication(fun_avg, ntraces_factor);
        // divide by num_traces to estimate Var(f(y_t)_i)
        fun_var = array_scalar_multiplication(fun_var, ntraces_factor);
        
        // calculate and print execution time
        auto stop = std::chrono::high_resolution_clock::now(); 
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count();
        printf("\nExecution time to simulate %d x %g s trace with dt = %.1e: %g s\n", 
                num_traces, t_interval[1] - t_interval[0], dt, ((float) duration)/1e6);
    }
    return 0;
}

// Print results   
// avg trace number i (where i is degree of freedom number i) will be
// printed (APPENDED) on file ./simulated_traces/sde_sample_path_statistic_i.txt
int print_results(unsigned int n_dim, const std::vector<std::vector<double>> fun_avg,
unsigned int subsampling_f, std::string statistic){
    std::vector<std::string> filename(n_dim);
    std::cout << "Traces saved in files:\n";
    for(unsigned int i = 0; i < n_dim; i++){
        // generate new filename
        filename[i] = "./simulated_traces/sde_sample_path_" + statistic + "_" + std::to_string(i) + ".txt";
        std::cout << filename[i] << "\n";
    }
    
    // Print results of dimension dim to file
    // Notice that, by default, this APPENDS a new row to the file
    for(unsigned int i = 0; i < n_dim; i++){
        // open file
        FILE* fp = fopen(filename[i].c_str(), "a"); // <-- a for "append"
        if (fp == NULL) {
          fprintf(stderr, "Can't open output file\n");
          exit(1);
        }
        // print to file
        print_array_asrow(fun_avg, i, subsampling_f, fp);    
        // Close file
        fclose(fp);
    }
    printf("------------------------------------------\n");
    return 0;
}

// fill problem parameters with inputs, if given (otherwise use default
// values)
int fill_parameters_w_inputs(int argc, char* argv[], double &dt, 
std::vector<double> &t_interval, unsigned int &num_traces, bool &many_traces, 
unsigned int &subsampling_f, std::string &eq_type, bool &Ito, unsigned int &n_dim, 
std::vector<double> &y0){
    // dt: time step
    if(argc < 2){
        dt = 1e-5; // default value
    } else {
        dt = atof(argv[1]); // input value
    }
    
    // time interval
    if(argc < 3){ 
        t_interval = {0, 1e-1}; // default value
    } else {
        t_interval = {atof(argv[2]), atof(argv[3])}; // input value
    }

    // number of traces
    if(argc < 5){ 
        num_traces = 100; // default value
    } else {
        num_traces = atoi(argv[4]);
        // case of few traces: not multithreading
        if(atoi(argv[4]) < 4){ 
            many_traces = false;
        // case of many traces: multithreading
        } else {   
            many_traces = true;
        }
    }

    // subsampling factor
    if(argc < 6){ 
        subsampling_f = 1; // don't do any subsampling when printing results
    } else {
        subsampling_f = atoi(argv[5]); // input value
    }
    
    // ODE or SDE
    if(argc < 7){ 
        eq_type = "sde"; // SDE by default
    } else {
        eq_type = argv[6]; // input value
    }
    
    // Ito or Stratonovich
    if(argc < 8){ 
        Ito = true; // Ito by default
    } else {
        Ito = atoi(argv[7]); // input value
    }
    
    // problem dimension
    if(argc < 9){ 
        n_dim = 2; // 2D problem, think of x and v.
    } else {
        n_dim = atoi(argv[8]); // input value
    }    
        
    // initial conditions, if given
    std::vector<double> v_zeros(n_dim);
    y0 = v_zeros;
    if(argc < 9){ 
        // default: initial conditions = 0
        y0 = {0, 0};
    } else {
        // read initial conditions
        for(unsigned int i = 0; i < n_dim; i++){
            y0[i] = atof(argv[9 + i]);
        }
    }
    return 0;
}
