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
#include "sde_library.h" // SDE library

// C++ function that mimics MATLAB linspace
std::vector<double> linspace(std::vector<double> interval, int n) {
    // preallocate output vector
    std::vector<double> v_linspace(n);
    // calculate h
    double h = (interval[1] - interval[0])/(n-1);
    // fill vector
    for(int i = 0; i < n + 1; i++) {
        v_linspace[i] = interval[0] + h*i;
    }
    return v_linspace;
}

// C++ function that mimics MATLAB a:h:b
std::vector<double> generate_v_from_h(std::vector<double> interval, double h) {
    // output vector
    std::vector<double> v;
    double current_value = interval[0];
    while(current_value <= interval[1]) {
        v.push_back(current_value);
        current_value = current_value + h;
    }
    return v;
}

// C++ function that mimics MATLAB zeros (but with arbitrary "value")
std::vector<std::vector<double>> ini_matrix(int n_arrays, int n_columns, double value){
    std::vector<std::vector<double>> y(n_arrays, std::vector<double>(n_columns, value));
    return y;
}

// C++ function to print array on stdout
int print_array(std::vector<std::vector <double>> a, std::FILE* fp) {
    for(int i = 0; i < a.size(); i++){
        for(int j = 0; j < a[0].size(); j++){
            fprintf(fp, "%.10f ", a[i][j]);
        }
        fprintf(fp, "\n");
    }
    return 0;
}

// C++ function to print array on stdout
int print_array_asrow(std::vector<std::vector <double>> a, int dim, std::FILE* fp) {
    for(int i = 0; i < a.size(); i++){
        fprintf(fp, "%.10f ", a[i][dim]);
    }
    fprintf(fp, "\n");
    return 0;
}

// C++ function to print vector on stdout
int print_vector(std::vector<double> v, std::FILE* fp) {
    for(int i = 0; i < v.size(); i++){
        fprintf(fp, "%f\n", v[i]);
    }
    return 0;
}

// Print "hola!" in stdout for debugging
int hola(){
    printf("hola!\n");
}

// sum of 2 vectors OF THE SAME SIZE
std::vector<double> vector_sum(std::vector<double> &a, std::vector<double> &b){
    std::vector<double> c(a.size());
    for(int i = 0; i < a.size(); i++){
        c[i] = a[i] + b[i];
    }
    return c;
}

// sum of 2 arrays OF THE SAME SIZE
std::vector<std::vector<double>> array_sum(std::vector<std::vector<double>> &a,
std::vector<std::vector<double>> &b){
        
    std::vector<std::vector<double>> c = ini_matrix(a.size(), a[0].size(), 0);
    for(int i = 0; i < a.size(); i++){
        for(int j = 0; j < a[0].size(); j++){
            c[i][j] = a[i][j] + b[i][j];
        }
    }
    return c;
}

// dot product of 2 vectors OF THE SAME SIZE
std::vector<double> dot_product(std::vector<double> &a, std::vector<double> &b){
    std::vector<double> c(a.size());
    for(int i = 0; i < a.size(); i++){
        c[i] = a[i]*b[i];
    }
    return c;
}

// dot product of 2 arrays OF THE SAME SIZE
std::vector<std::vector<double>> array_dot_product(std::vector<std::vector<double>> &a,
std::vector<std::vector<double>> &b){
            
    std::vector<std::vector<double>> c = ini_matrix(a.size(), a[0].size(), 0);
    for(int i = 0; i < a.size(); i++){
        for(int j = 0; j < a[0].size(); j++){
            c[i][j] = a[i][j]*b[i][j];
        }
    }
    return c;
}

// function that is applied element by element in function_array
double f(double x){
    return x*x;
}

// apply f(x) on array element by element
std::vector<std::vector<double>> function_array(std::vector<std::vector<double>> &a){
            
    std::vector<std::vector<double>> c = ini_matrix(a.size(), a[0].size(), 0);
    for(int i = 0; i < a.size(); i++){
        for(int j = 0; j < a[0].size(); j++){
            c[i][j] = f(a[i][j]);
        }
    }
    return c;
}


// product of vector times scalar
std::vector<double> scalar_multiplication(std::vector<double> &a, double k){
    std::vector<double> result(a.size());
    for(int i = 0; i < a.size(); i++){
        result[i] = k*a[i];
    }
    return result;
}

// product of array times scalar
std::vector<std::vector<double>> array_scalar_multiplication(std::vector<std::vector<double>> &a, double k){
    std::vector<std::vector<double>> result = ini_matrix(a.size(), a[0].size(), 0);
    for(int i = 0; i < a.size(); i++){
        for(int j = 0; j < a[0].size(); j++){
            result[i][j] = k*a[i][j];
        }
    }
    return result;
}

// drift of process
std::vector<double> drift_function(std::vector<double> y, double t){
    // output vector
    std::vector<double> f(y.size(), 0);
    
    // use physical names
    double x = y[0];
    double v = y[1];    
    
    // Definitions: experimental parameters
    double m = 9.2e-18;
    double gamma = 3.5e-11;
    double w = 2*M_PI*2000;
    double eps = 6.3e-5;
    double T = 295;
    double k_B = 1.38065e-23;
    double sigma = sqrt(2*k_B*T*gamma);
    double gamma_m = gamma/m;
    double eps_m = eps/m;
    
    // i.e., x' = v
    f[0] = v; 
    f[1] = 0; //(-gamma_m*v + eps_m*std::cos(w*t)*x);//1/m*(-gamma*v + eps*std::cos(w*t)*x);
        
    // return output
    return f;
}

// diffusion of process
std::vector<double> diffusion_function(std::vector<double> y, double t){
    // output vector
    std::vector<double> f(y.size(), 0);
    
    // use physical names
    double x = y[0];
    double v = y[1];
    
    // Definitions: experimental parameters
    double m = 9.2e-18;
    double gamma = 3.5e-11;
    double w = 2*M_PI*2e4;
    double eps = 6.3e-9;
    double T = 295;
    double k_B = 1.38065e-23;
    double sigma = sqrt(2*k_B*T*gamma);
    
    // Calculate force
    // No noise in x:
    // f[0] = 0; -> but already 0 by initialization
    // Stochastic force in momentum:
    f[1] = sigma/m;
    
    // return output
    return f;
}

// Roll a dice between low and high
int roll_dice(int low, int high){
    int dice_result = rand() % (high - low + 1) + low;
    return dice_result;
}

// Runge-Kutta function as programmed in MATLAB
// output in y
void runge_kutta(std::vector<double> t_interval,
std::vector<double> y0, double dt, std::vector<std::vector<double>> &y){
    
    // PARAMETERS AND INITIALIZATION
    // generate time vector
    
    std::vector<double> t = generate_v_from_h(t_interval, dt);
    int n_time = t.size();
    // preallocate solution vector y, full of zeros
    y = ini_matrix(n_time, y0.size(), 0);

    // initialize y with initial conditions as y(0,:) = y0
    y[0] = y0;
    
    // S_k for the runge_kutta method
    std::vector<double> sk_vec{-1, 1};
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
    for(int n = 1; n < n_time; n++){
        // calculate dW_t
        for(int i = 0; i < dW_t.size(); i++){
            dW_t[i] = n_distribution(generator)*sqrt(dt);
        }
        
        // calculate S_k
        S_k = sk_vec[roll_dice(0,1)];
        
        // Calculate k_1
        tmp1 = drift_function(y.back(), t[n - 1]);
        tmp1 = scalar_multiplication(tmp1, dt);

        for(int i = 0; i < dW_t.size(); i++){
            tmp2[i] = dW_t[i] - S_k*sqrt(dt);
        }

        tmp3 = diffusion_function(y[n - 1], t[n - 1]);
        tmp2 = dot_product(tmp2, tmp3);

        k_1 = vector_sum(tmp1, tmp2);

        // Calculate k_2
        tmp1 = drift_function(vector_sum(y[n - 1], k_1), t[n]);
        tmp1 = scalar_multiplication(tmp1, dt);
        
        for(int i = 0; i < dW_t.size(); i++){
            tmp2[i] = dW_t[i] + S_k*sqrt(dt);
        }
        
        tmp3 = diffusion_function(y[n - 1], t[n]);      
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
void generate_traces(int num_traces, std::string filename, 
std::vector<double> t_interval, std::vector<double> y0, double dt){
    // open to write to file
    FILE* fp = fopen(filename.c_str(), "a");
    if (fp == NULL) {
      fprintf(stderr, "Can't open output file\n");
      exit(1);
    }
    // preallocate solution vector y
    std::vector<std::vector<double>> y;
    // call RK method num_traces times
    for(int i = 0; i < num_traces; i++){
        runge_kutta(t_interval, y0, dt, y);  
        // Print results to file 
        print_array_asrow(y, 0, fp);
    }
    // Close file and exit
    fclose(fp);
}

// Generate num_traces traces with the RK method and return 
// sum of function <f(Y_t)> (Obs: remains to divide by num_traces
// to estimate average)
void generate_avg_trace(int num_traces, std::vector<double> t_interval, 
std::vector<double> y0, std::vector<std::vector<double>> &avg_var, double dt){

    // preallocate solution vector y and tmp1 vector
    std::vector<std::vector<double>> y;
    std::vector<std::vector<double>> tmp1;
    
    // call RK method num_traces times
    runge_kutta(t_interval, y0, dt, y); 
    avg_var = array_dot_product(y, y);
    for(int i = 1; i < num_traces; i++){
        runge_kutta(t_interval, y0, dt, y); 
        tmp1 = function_array(y);
        avg_var = array_sum(avg_var, tmp1);
    }    
}
