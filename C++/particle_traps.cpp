/********************************************************************** 
 * DESCRIPTION:
 * 
 * Particle trap libraries: optical tweezers (i.e. dipole trap),  
 * Paul trap, Gaussian beams, etc.
 * 
 *********************************************************************** 
 * OBSERVATIONS:
 * 
 * 
 *********************************************************************** 
 * REFERENCES:
 * 
 * GIT REPOSITORY
 * https://github.com/gerardpc/sde_simulator
 ***********************************************************************  
 * Versions: 
 *  By GP Conangla
 *  04.10.2019
 *      Parts are working, field functions not yet tested.
 *********************************************************************** 
 */

// STANDARD LIBRARIES
#include <iostream> // stdin, stdout 
#include <cstdio> // printf family
#include <cmath> // most math functions
#include <vector> // vector library
#include <string> // string library
#include <stdlib.h> // includes rand()
#include <fstream> // for reading file

// MY LIBRARIES
#include "particle_traps.h"


//======================================================================
// Struct member functions
//======================================================================
// Print functions
void gaussian_beam::print(){
    printf("length_unit: %.4e\n", length_unit);
    printf("h_bar: %.4e\n", h_bar);
    printf("c: %.4e\n", c);
    printf("eps_0: %.4e\n", eps_0);
    printf("lambda: %.4e\n", lambda);
    printf("P_0: %.4e\n", P_0);
    printf("NA: %.4e\n", NA);
    printf("w_0: %.4e\n", w_0);
    printf("z_R: %.4e\n", z_R);
    printf("I_0: %.4e\n", I_0);
}

void paul_trap::print(){
    printf("r: %.4e\n", r);
    printf("density: %.4e\n", density);
    printf("m: %.4e\n", m);
    printf("n: %.4e\n", n);
    printf("alpha: %.4e\n", alpha);
    printf("f_dr: %.4e\n", f_dr);
    printf("w_dr: %.4e\n", w_dr);
    printf("Q: %.4e\n", Q);
    printf("V: %.4e\n", V);
    printf("d: %.4e\n", d);
    printf("eps: %.4e\n", eps);
    printf("T: %.4e\n", T);
    printf("pressure: %.4e\n", pressure);
    printf("gamma_ambient: %.4e\n", gamma_ambient);
    printf("gamma: %.4e\n", gamma);
    printf("g_norm: %.4e\n", g_norm);
    printf("sigma: %.4e\n", sigma);
}

void opt_tweezer::print(){
    printf("r: %.4e\n", r);
    printf("density: %.4e\n", density);
    printf("m: %.4e\n", m);
    printf("n: %.4e\n", n);
    printf("alpha: %.4e\n", alpha);
    printf("freq: %.4e\n", freq);
    printf("w: %.4e\n", w);
    printf("T: %.4e\n", T);
    printf("pressure: %.4e\n", pressure);
    printf("gamma_ambient: %.4e\n", gamma_ambient);
    printf("gamma: %.4e\n", gamma);
    printf("g_norm: %.4e\n", g_norm);
    printf("sigma: %.4e\n", sigma);
}

// Constructor functions
void gaussian_beam::fill(){        
    // gaussian beam file
    std::ifstream input_file("./gaussian_beam.txt");
    
    // vector of doubles containing values from file
    std::vector<double> file_values;
    
    // read from file and put values in file_values
    std::string member, description;
    double value;
    while(input_file >> member >> value >> description){
        file_values.push_back(value);
    }
    input_file.close();
    
    // Physical constants
    length_unit = 1e-9; 
    h_bar = 1.0546e-34;
    c = 3e8;
    eps_0 = 8.854e-12;
    
    // laser beam
    lambda = file_values[0]*length_unit; // laser wavelength 
    P_0 = file_values[1]*1e-3; // laser power in W
    NA = file_values[2]; // Gaussian beam NA
    
    // numerical details    
    double c_factor = 1.5; // correction factor for high NA in Gaussian beam (set to "1" to eliminate correction)
    
    // Gaussian beam
    w_0 = lambda/(M_PI*NA)*c_factor;
    z_R = M_PI*pow(w_0, 2)/lambda;
    I_0 = P_0*2/(M_PI*pow(w_0, 2));
}

void paul_trap::fill(){
    // paul trap file
    std::ifstream input_file("./paul_trap.txt");
    
    // vector of doubles containing values from file
    std::vector<double> file_values;
    
    // read from file and put values in file_values
    std::string member, description;
    double value;
    while(input_file >> member >> value >> description){
        file_values.push_back(value);
    }
    input_file.close();
    
    // Write values on struct members
    // Physical constants
    double k_B = 1.38065e-23;
    double eps_0 = 8.854e-12;
    // particle parameters
    r = file_values[0];
    density = file_values[1];
    m = 4./3.*M_PI*pow(r, 3)*density;
    n = file_values[2];
    alpha = 4*M_PI*pow(r, 3)*eps_0*(pow(n, 2) - 1)/(pow(n, 2) + 2); // induced bulk dipole
    // trap parameters
    f_dr = file_values[3];
    w_dr = 2*M_PI*f_dr;
    Q = file_values[6];
    V = file_values[7];
    d = file_values[8];
    eps = 1.602176620e-19*Q*V/pow(d,2);
    // ambient parameters
    T = file_values[4];
    pressure = file_values[5];
    gamma_ambient = 6*M_PI*r*1.84e-5;
    gamma = gamma_ambient*pressure/1010;
    g_norm = gamma/m;
    sigma = sqrt(2*k_B*T*gamma);
    
    // print value of beta, as defined in Iz. et al PRE 1995
    double beta = 4*Q*V*1.602e-19/(m*pow(d*w_dr, 2));
    printf("Value of Beta: %.4e\n", beta);
    printf("Value of eps: %.4e\n", eps);
}

void opt_tweezer::fill(){
    // optical tweezer file
    std::ifstream input_file("./optical_tweezer.txt");
    
    // vector of doubles containing values from file
    std::vector<double> file_values;
    
    // read from file and put values in file_values
    std::string member, description;
    double value;
    while(input_file >> member >> value >> description){
        file_values.push_back(value);
    }
    input_file.close();
    
    // Write values on struct members
    // Physical constants
    double k_B = 1.38065e-23;
    double eps_0 = 8.854e-12;
    // particle parameters
    r = file_values[0];
    density = file_values[1];
    m = 4./3.*M_PI*pow(r, 3)*density;
    n = file_values[2];
    alpha = 4*M_PI*pow(r, 3)*eps_0*(pow(n, 2) - 1)/(pow(n, 2) + 2); // induced bulk dipole
    // trap parameters
    freq = file_values[3];
    w = 2*M_PI*freq;
    // ambient parameters
    T = file_values[4];
    pressure = file_values[5];
    gamma_ambient = 6*M_PI*r*1.84e-5;
    gamma = gamma_ambient*pressure/1010;
    g_norm = gamma/m;
    sigma = sqrt(2*k_B*T*gamma);
}

//======================================================================
// FIELD FUNCTIONS
//======================================================================
// calculate particle mass of an assumed spherical silica particle
// r in meters, density in kg/m^3
double particle_mass(double r, double density){
    return 4./3.*M_PI*pow(r, 3)*density;
}

// Gaussian beam width
double gb_w(double z, gaussian_beam &gb){
    return gb.w_0*sqrt(1 + pow((z/gb.z_R), 2));
}

// Gaussiain beam intensity
double gb_I(double r, double z, gaussian_beam &gb){
    return gb.I_0*pow((gb.w_0/gb_w(z, gb)),2)*exp(-(2*pow(r/gb_w(z, gb), 2)));
}

// Intensity gradient in the r direction
double grad_I_r(double r, double z, gaussian_beam &gb, double h){
    return (1./12.*gb_I(r-2*h, z, gb) - 2./3.*gb_I(r-h, z, gb) 
    + 2./3.*gb_I(r+h, z, gb) - 1./12.*gb_I(r+2*h, z, gb))/h;
}

// Intensity gradient in the z direction
double grad_I_z(double r, double z, gaussian_beam &gb, double h){
    return (1./12.*gb_I(r, z-2*h, gb) - 2./3.*gb_I(r, z-h, gb)
    + 2./3.*gb_I(r, z+h, gb) - 1./12.*gb_I(r, z+2*h, gb))/h;
}

// Dipole force f_r(r,z), Gaussian beam 
// for a Rayleigh particle (a << lambda) or atom
double force_r_gb(double r, double z, double alpha, gaussian_beam &gb, double h){
    return 1./2.*alpha*grad_I_r(r,z,gb,h);
}

// Dipole force f_z(r,z), Gaussian beam 
// for a Rayleigh particle (a << lambda) or atom
double force_z_gb(double r, double z, double alpha, gaussian_beam &gb, double h){
    return 1./2.*alpha*grad_I_z(r,z,gb,h);
}

// paul trap force field
double force_paul_trap(double x, double t, paul_trap &pt){
    return (pt.eps*cos(pt.w_dr*t)*x);
}

//======================================================================
// FUNCTIONS HERE
//======================================================================
// Main function for testing
int main(){
    gaussian_beam gb;
    paul_trap pt;
    opt_tweezer ot;
    
    gb.fill();
    gb.print();
    printf("\n");
    
    pt.fill();
    pt.print();
    printf("\n");
    
    ot.fill();
    ot.print();
    
    double a = 150e-9;
    
    double h = 1e-5*a; // discretization size
    return 0;
}




























