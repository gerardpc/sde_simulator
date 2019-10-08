/********************************************************************** 
 * DESCRIPTION:
 * 
 * Basic library with elementary vector operations and manipulations
 *********************************************************************** 
 * OBSERVATIONS:
 * 
 *
 *********************************************************************** 
 * REFERENCES:
 * 
 * GIT REPOSITORY
 * https://github.com/gerardpc/sde_simulator
 *
 * **********************************************************************  
 * Versions: 
 *  By GP Conangla
 *  08.10.2019
 *          Working
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

// MY LIBRARIES
#include "num_vector.hpp" 

// Function f(Y_t) that is applied to every trace and averaged to
// estimate <f(Y_t)>
double f(double x){
    return x*x; // square function by default, to estimate variance
}

// C++ function that mimics MATLAB linspace
std::vector<double> linspace(std::vector<double> interval, unsigned int n) {
    // preallocate output vector
    std::vector<double> v_linspace(n);
    // calculate h
    double h = (interval[1] - interval[0])/(n-1);
    // fill vector
    for(unsigned int i = 0; i < n + 1; i++) {
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
int print_array(const std::vector<std::vector <double>> a, std::FILE* fp) {
    for(unsigned int i = 0; i < a.size(); i++){
        for(unsigned int j = 0; j < a[0].size(); j++){
            fprintf(fp, "%.10f ", a[i][j]);
        }
        fprintf(fp, "\n");
    }
    return 0;
}

// C++ function to print array on stdout
int print_array_asrow(const std::vector<std::vector <double>> a, unsigned int dim, std::FILE* fp) {
    for(unsigned int i = 0; i < a.size(); i++){
        fprintf(fp, "%.10f ", a[i][dim]);
    }
    fprintf(fp, "\n");
    return 0;
}

// C++ function to print vector on stdout
int print_vector(const std::vector<double> v, std::FILE* fp) {
    for(unsigned int i = 0; i < v.size(); i++){
        fprintf(fp, "%f\n", v[i]);
    }
    return 0;
}

// Print "hola!" in stdout for debugging
void hola(){
    printf("hola!\n");
}

// sum of 2 vectors OF THE SAME SIZE
std::vector<double> vector_sum(const std::vector<double> &a, const std::vector<double> &b){
    std::vector<double> c(a.size());
    for(unsigned int i = 0; i < a.size(); i++){
        c[i] = a[i] + b[i];
    }
    return c;
}

// sum of 2 arrays OF THE SAME SIZE
std::vector<std::vector<double>> array_sum(const std::vector<std::vector<double>> &a,
const std::vector<std::vector<double>> &b){
        
    std::vector<std::vector<double>> c = ini_matrix(a.size(), a[0].size(), 0);
    for(unsigned int i = 0; i < a.size(); i++){
        for(unsigned int j = 0; j < a[0].size(); j++){
            c[i][j] = a[i][j] + b[i][j];
        }
    }
    return c;
}

// dot product of 2 vectors OF THE SAME SIZE
std::vector<double> dot_product(const std::vector<double> &a, const std::vector<double> &b){
    std::vector<double> c(a.size());
    for(unsigned int i = 0; i < a.size(); i++){
        c[i] = a[i]*b[i];
    }
    return c;
}

// dot product of 2 arrays OF THE SAME SIZE
std::vector<std::vector<double>> array_dot_product(const std::vector<std::vector<double>> &a,
std::vector<std::vector<double>> &b){
            
    std::vector<std::vector<double>> c = ini_matrix(a.size(), a[0].size(), 0);
    for(unsigned int i = 0; i < a.size(); i++){
        for(unsigned int j = 0; j < a[0].size(); j++){
            c[i][j] = a[i][j]*b[i][j];
        }
    }
    return c;
}

// apply f(x) on array element by element
std::vector<std::vector<double>> function_array(const std::vector<std::vector<double>> &a){
            
    std::vector<std::vector<double>> c = ini_matrix(a.size(), a[0].size(), 0);
    for(unsigned int i = 0; i < a.size(); i++){
        for(unsigned int j = 0; j < a[0].size(); j++){
            c[i][j] = f(a[i][j]);
        }
    }
    return c;
}


// product of vector times scalar
std::vector<double> scalar_multiplication(const std::vector<double> &a, double k){
    std::vector<double> result(a.size());
    for(unsigned int i = 0; i < a.size(); i++){
        result[i] = k*a[i];
    }
    return result;
}

// product of array times scalar
std::vector<std::vector<double>> array_scalar_multiplication(const std::vector<std::vector<double>> &a, double k){
    std::vector<std::vector<double>> result = ini_matrix(a.size(), a[0].size(), 0);
    for(unsigned int i = 0; i < a.size(); i++){
        for(unsigned int j = 0; j < a[0].size(); j++){
            result[i][j] = k*a[i][j];
        }
    }
    return result;
}

// Roll a dice between low and high
int roll_dice(int low, int high){
    int dice_result = rand() % (high - low + 1) + low;
    return dice_result;
}
