/********************************************************************** 
 * DESCRIPTION:
 * 
 * Library to generate sample paths of a given stochastic 
 * process, defined by a user-defined SDE
 * 
 * dY = a(t,Y)dt + b(t,Y)*dW_t, Y(t0) = Y0,
 * 
 * using a Runge-Kutta type method for SDE of strong order 1 and
 * deterministic order 2.
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
// MATH PROBLEM DEFINITION FUNCTIONS
//======================================================================
// SDE FUNCTIONS
// dY = a(t,Y)dt + b(t,Y)*dW_t, Y(t0) = Y0,
// drift of process: a(t,Y)

std::vector<double> drift_function(std::vector<double> y, double t, const eq_params &eq){
    // output vector
    std::vector<double> f(y.size(), 0);
    
    // use physical names
    double x = y[0];
    double v = y[1];    
    
    // Definitions: experimental parameters
    
    // x' = v
    f[0] =  v;
    // v' = ...
    //f[1] = force_r_gb(x, 0, eq.ot.alpha, eq.gb, 1e-5*eq.ot.r)/eq.ot.m - eq.ot.g_norm*v;
    //f[1] = -pow(eq.ot.w_gb_r, 2)*x - eq.ot.g_norm*v; 
    //f[1] = - eq.ot.g_norm*v; 
    f[1] = -eq.th.g_norm*v + eq.pt.eps*cos(eq.pt.w_dr*t)*x/eq.part.m;
        
    // return output
    return f;
}

// diffusion of process: b(t, Y)
std::vector<double> diffusion_function(std::vector<double> y, double t, const eq_params &eq){
    // output vector
    std::vector<double> f(y.size(), 0);
    
    // use physical names
    double x = y[0];
    double v = y[1];
    
    // Calculate force
    // No noise in x:
    f[0] = 0;
    // Stochastic force in momentum:
    f[1] = eq.th.sigma/eq.part.m;
    
    // return output
    return f;
}


//======================================================================
// Equation parameters struct constructor
//======================================================================
void eq_params::fill(){
    gb.fill();
    part.fill();
    th.fill(part);
    ot.fill();
    ot.fill_gb_w(part, gb);
    pt.fill(part, th);
}

void eq_params::print(){
    printf("Equation parameters:\n\n");
    printf("Gaussian beam\n");
    gb.print();
    printf("\nParticle\n");
    part.print();
    printf("\nThermodynamics\n");
    th.print();
    printf("\nOptical tweezer\n");
    ot.print();
    printf("\nPaul trap\n");
    pt.print();
}

// Dipole force f_r(r,z), Gaussian beam 
// for a Rayleigh particle (a << lambda) or atom
double eq_params::force_r(double r, double z){
    return force_r_gb(r, z, part.alpha, gb, part.h);
}

// Dipole force f_z(r,z), Gaussian beam 
// for a Rayleigh particle (a << lambda) or atom
double eq_params::force_z(double r, double z){
    return force_z_gb(r, z, part.alpha, gb, part.h);
}

//======================================================================
// LIBRARY FUNCTIONS
//======================================================================
// Runge-Kutta function: numerical code.
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
std::vector<double> y0, std::vector<std::vector<double>> &avg_var, double dt, 
bool Ito, const eq_params &args, unsigned int thread_id){
    // preallocate solution vector y and tmp1 vector
    std::vector<std::vector<double>> y;
    std::vector<std::vector<double>> tmp1;
    
    // call RK method num_traces times
    runge_kutta(t_interval, y0, dt, y, Ito, args); 
    avg_var = function_array(y);
    if(thread_id == 0){ // main thread, print progress
        unsigned int progress;
        unsigned int last_progress = 0;
        // Print progress bar
        std::cout << "\nProgress:\n";
        std::cout << "0 " << std::flush;
        // Calculate rest of traces
        for(unsigned int i = 1; i < num_traces; i++){
            runge_kutta(t_interval, y0, dt, y, Ito, args); 
            tmp1 = function_array(y);
            avg_var = array_sum(avg_var, tmp1);
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
            tmp1 = function_array(y);
            avg_var = array_sum(avg_var, tmp1);
        }
    }
}

// Generate num_traces traces with generate_avg_trace using different
// threads (automatic core detection. 
// Print elapsed time for execution on stdout. Returns average
// of f(Y_t)
std::vector<std::vector<double>> RK_all(unsigned int num_traces, bool many_traces, 
std::vector<double> t_interval, std::vector<double> y0, double dt, bool Ito, const eq_params &args){
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
    
    // preallocate average of f(x) trace vectors 
    std::vector<std::vector<double>> avg_trace;
    
    // Run parallel RK methods to speed code x num_threads 
    if(many_traces and num_threads >= 1){ // if multithreading
        // calculate traces per core
        num_traces = num_traces - num_traces % num_cores;
        printf("Number of traces that will be generated: %d\n", num_traces);
        int traces_per_core = num_traces/num_cores;
        printf("Number of traces per core: %d\n", traces_per_core);
        
        // vector of RK results and threads
        std::vector<std::vector<std::vector<double>>> thread_avg_trace(num_threads);
        std::vector<std::thread> t;
        
        // initiate simulations in threads
        for(unsigned int i = 0; i < num_threads; i++){
            t.push_back(std::thread(generate_avg_trace, traces_per_core, t_interval, y0, ref(thread_avg_trace[i]), dt, Ito, args, i + 1));
        }
        
        // run 1 more simulation in main
        generate_avg_trace(traces_per_core, t_interval, y0, avg_trace, dt, Ito, args, 0);  
    
        // join threads when they finish
        for(unsigned int i = 0; i < num_threads; i++){
            t[i].join();
        }
    
        // Calculate sum of variances and put it on avg_trace
        for(unsigned int i = 0; i < num_threads; i++){
            avg_trace = array_sum(thread_avg_trace[i], avg_trace);
        }
        
        // divide by num_traces to estimate <f(y_t)_i>
        double num_tracesf = num_traces;
        double ntraces_factor = 1/(num_tracesf);
        avg_trace = array_scalar_multiplication(avg_trace, ntraces_factor);
        
        // calculate and print execution time
        auto stop = std::chrono::high_resolution_clock::now(); 
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count();
        printf("\n\nExecution time to simulate %d x %g s trace with dt = %.1e: %g s\n", 
                num_traces, t_interval[1] - t_interval[0], dt, ((float) duration)/1e6);

    } else {  // not multithreading
        printf("Number of traces that will be generated: %d\n", num_traces);
        // generate single thread
        generate_avg_trace(num_traces, t_interval, y0, avg_trace, dt, Ito, args, 0);  
        // divide by num_traces to estimate <f(y_t)_i(t)>
        double num_tracesf = num_traces;
        double ntraces_factor = 1/(num_tracesf);
        avg_trace = array_scalar_multiplication(avg_trace, ntraces_factor);
        
        // calculate and print execution time
        auto stop = std::chrono::high_resolution_clock::now(); 
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count();
        printf("\nExecution time to simulate %d x %g s trace with dt = %.1e: %g s\n", 
                num_traces, t_interval[1] - t_interval[0], dt, ((float) duration)/1e6);
    }        
    return avg_trace;
}

// Print results   
// avg trace number i (where i is degree of freedom number i) will be
// printed on file ./simulated_traces/sde_sample_path_i.txt
int print_results(unsigned int n_dim, const std::vector<std::vector<double>> avg_trace,
unsigned int subsampling_f){
    std::vector<std::string> filename(n_dim);
    std::cout << "Average traces saved in files:\n";
    for(unsigned int i = 0; i < n_dim; i++){
        // generate new filename
        filename[i] = "./simulated_traces/sde_sample_path_" + std::to_string(i) + ".txt";
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
        print_array_asrow(avg_trace, i, subsampling_f, fp);    
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
unsigned int &subsampling_f, bool &Ito, unsigned int &n_dim, std::vector<double> &y0){
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
    
    // Ito or Stratonovich
    if(argc < 7){ 
        Ito = true; // Ito by default
    } else {
        Ito = atoi(argv[6]); // input value
    }
    
    // problem dimension
    if(argc < 8){ 
        n_dim = 2; // 2D problem, think of x and v.
    } else {
        n_dim = atoi(argv[7]); // input value
    }    
        
    // initial conditions, if given
    std::vector<double> v_zeros(n_dim);
    y0 = v_zeros;
    if(argc < 8){ 
        // default: initial conditions = 0
        y0 = {0, 0};
    } else {
        // read initial conditions
        for(unsigned int i = 0; i < n_dim; i++){
            y0[i] = atof(argv[8 + i]);
        }
    }
    return 0;
}
