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
#include <vector> // vector library
#include <cstdio> // printf family
 
// Function f(Y_t) that is applied to every trace and averaged to
// estimate <f(Y_t)>
double f(double x);

// C++ function that mimics MATLAB linspace
std::vector<double> linspace(std::vector<double> interval, unsigned int n);

// C++ function that mimics MATLAB a:h:b
std::vector<double> generate_v_from_h(std::vector<double> interval, double h);

// C++ function that mimics MATLAB zeros (but with arbitrary "value")
std::vector<std::vector<double>> ini_matrix(int n_arrays, int n_columns, double value);

// C++ function to print array on stdout
int print_array(const std::vector<std::vector <double>> a, std::FILE* fp);

// C++ function to print array on stdout
int print_array_asrow(const std::vector<std::vector <double>> a, unsigned int dim, std::FILE* fp);

// C++ function to print vector on stdout
int print_vector(const std::vector<double> v, std::FILE* fp);

// Print "hola!" in stdout for debugging
void hola();

// sum of 2 vectors OF THE SAME SIZE
std::vector<double> vector_sum(const std::vector<double> &a, const std::vector<double> &b);

// sum of 2 arrays OF THE SAME SIZE
std::vector<std::vector<double>> array_sum(const std::vector<std::vector<double>> &a,
const std::vector<std::vector<double>> &b);

// dot product of 2 vectors OF THE SAME SIZE
std::vector<double> dot_product(const std::vector<double> &a, const std::vector<double> &b);

// dot product of 2 arrays OF THE SAME SIZE
std::vector<std::vector<double>> array_dot_product(const std::vector<std::vector<double>> &a,
std::vector<std::vector<double>> &b);

// apply f(x) on array element by element
std::vector<std::vector<double>> function_array(const std::vector<std::vector<double>> &a);


// product of vector times scalar
std::vector<double> scalar_multiplication(const std::vector<double> &a, double k);

// product of array times scalar
std::vector<std::vector<double>> array_scalar_multiplication(const std::vector<std::vector<double>> &a, double k);

// Roll a dice between low and high
int roll_dice(int low, int high);

