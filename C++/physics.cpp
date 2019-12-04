/********************************************************************** 
 * DESCRIPTION:
 * 
 * Particle trap libraries: equation parameters, 
 * optical tweezers (i.e. dipole trap),  Paul trap, Gaussian beams, etc.
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
 *  14.10.2019
 *      Working.
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
#include "physics.hpp"

//======================================================================
// Struct member functions
//======================================================================
// Print functions
void gaussian_beam::print(){
    printf("length_unit: %.4e\n", length_unit);
    printf("h_bar: %.4e\n", h_bar);
    printf("c: %.4e\n", c);
    printf("eps_0: %.4e\n", eps_0);
    printf("c_factor: %.4e\n", c_factor);
    printf("lambda: %.4e\n", lambda);
    printf("P_0: %.4e\n", P_0);
    printf("NA: %.4e\n", NA);
    printf("w_0: %.4e\n", w_0);
    printf("z_R: %.4e\n", z_R);
    printf("I_0: %.4e\n", I_0);
}

void particle::print(){
    printf("r: %.4e\n", r);
    printf("density: %.4e\n", density);
    printf("m: %.4e\n", m);
    printf("n: %.4e\n", n);
    printf("alpha: %.4e\n", alpha);
    printf("Q: %.4e\n", Q);
    printf("h: %.4e\n", h);
}

void thermodynamics::print(){    
    printf("k_B: %.4e\n", k_B); 
    printf("T: %.4e\n", T);
    printf("pressure: %.4e\n", pressure);
    printf("gamma_ambient: %.4e\n", gamma_ambient);
    printf("gamma: %.4e\n", gamma);
    printf("g_norm: %.4e\n", g_norm);
    printf("sigma: %.4e\n", sigma);
}

void opt_trap::print(){
    printf("freq: %.4e\n", freq);
    printf("w: %.4e\n", w);    
    printf("f_r GB: %.4e\n", f_gb_r);
    printf("w_r GB: %.4e\n", w_gb_r);
    printf("f_z GB: %.4e\n", f_gb_z);
    printf("w_z GB: %.4e\n", w_gb_z);
}

void paul_trap::print(){    
    printf("f_dr: %.4e\n", f_dr);
    printf("w_dr: %.4e\n", w_dr);
    printf("V: %.4e\n", V);
    printf("d: %.4e\n", d);
    printf("eps: %.4e\n", eps);
    printf("alpha_iz: %.4e\n", alpha_iz);
    printf("beta_iz: %.4e\n", beta_iz);
}

// Constructor functions
void gaussian_beam::fill(){        
    // gaussian beam file
    std::ifstream input_file("../eq_params/gaussian_beam.txt");    
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
    
    // Write values on struct members   
    // laser beam
    lambda = file_values[0]*length_unit; // laser wavelength 
    P_0 = file_values[1]*1e-3; // laser power in W
    NA = file_values[2]; // Gaussian beam NA    
    // numerical details    
    c_factor = 1.5; // correction factor for high NA in Gaussian beam (set to "1" to eliminate correction)
    // Gaussian beam
    w_0 = lambda/(M_PI*NA)*c_factor;
    z_R = M_PI*pow(w_0, 2)/lambda;
    I_0 = P_0*2/(M_PI*pow(w_0, 2));
}

void particle::fill(){    
    // particle file
    std::ifstream input_file("../eq_params/particle.txt");    
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
    double eps_0 = 8.854e-12;
    
    // Write values on struct members
    r = file_values[0];
    density = file_values[1];
    m = 4./3.*M_PI*pow(r, 3)*density;
    n = file_values[2];
    alpha = 4*M_PI*pow(r, 3)*eps_0*(pow(n, 2) - 1)/(pow(n, 2) + 2); // induced bulk dipole
    Q = file_values[3];
    h = r*1e-4;
}

void thermodynamics::fill(particle &pt){
    // thermodynamics file
    std::ifstream input_file("../eq_params/thermodynamics.txt");    
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
    T = file_values[0];
    pressure = file_values[1];
    gamma_ambient = 6*M_PI*pt.r*1.84e-5; // Stokes drag
    gamma = gamma_ambient*pressure/1010; // gamma is proport. to pressure
    g_norm = gamma/pt.m;
    sigma = sqrt(2*k_B*T*gamma);
}

void opt_trap::fill(){
    // optical tweezer file
    std::ifstream input_file("../eq_params/optical_trap.txt");    
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
    // trap parameters
    freq = file_values[0];
    w = 2*M_PI*freq;
}

// calculate frequencies from Gaussian Beam
// ref: Gerard's theory notes
void opt_trap::fill_gb_w(particle &part, gaussian_beam &gb){        
    double k_r = 2*M_PI*pow(part.r,3)/gb.c*(pow(part.n,2) - 1)/(pow(part.n,2)+ 2)*8*gb.P_0/M_PI*pow(M_PI*gb.NA/(gb.lambda*gb.c_factor), 4);
    double k_z = 2*M_PI*pow(part.r,3)/gb.c*(pow(part.n,2) - 1)/(pow(part.n,2) + 2)*4*gb.P_0*pow(gb.lambda*gb.c_factor, 2)/pow(M_PI,3)*pow(M_PI*gb.NA/(gb.lambda*gb.c_factor), 6);
    w_gb_r = pow(k_r/part.m, 0.5);
    w_gb_z = pow(k_z/part.m, 0.5);
    f_gb_r = w_gb_r/(2*M_PI);
    f_gb_z = w_gb_z/(2*M_PI);
}

void paul_trap::fill(particle &part, thermodynamics &th){
    // paul trap file
    std::ifstream input_file("../eq_params/paul_trap.txt");    
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
    double e_charge = 1.602176620e-19; // electron charge (C)
    // trap parameters
    f_dr = file_values[0];
    w_dr = 2*M_PI*f_dr;
    V = file_values[1];
    d = file_values[2];
    eps = e_charge*part.Q*V/pow(d,2);    
    // value of alpha, beta, as defined in Iz. et al PRE 1995
    alpha_iz = 2*th.gamma/(part.m*w_dr);
    beta_iz = 4*part.Q*V*e_charge/(part.m*pow(d*w_dr, 2));
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

//======================================================================
// FIELD FUNCTIONS
//======================================================================
// calculate particle mass of an assumed spherical silica particle
// r in meters, density in kg/m^3
double particle_mass(double r, double density){
    return 4./3.*M_PI*pow(r, 3)*density;
}

// Gaussian beam width
double gb_w(double z, const gaussian_beam &gb){
    return gb.w_0*sqrt(1 + pow((z/gb.z_R), 2));
}

// Gaussian beam intensity
double gb_I(double r, double z, const gaussian_beam &gb){
    return gb.I_0*pow((gb.w_0/gb_w(z, gb)),2)*exp(-(2*pow(r/gb_w(z, gb), 2)));
}

// Intensity gradient in the r direction
double grad_I_r(double r, double z, const gaussian_beam &gb, double h){
    return (1./12.*gb_I(r-2*h, z, gb) - 2./3.*gb_I(r-h, z, gb) 
    + 2./3.*gb_I(r+h, z, gb) - 1./12.*gb_I(r+2*h, z, gb))/h;
}

// Intensity gradient in the z direction
double grad_I_z(double r, double z, const gaussian_beam &gb, double h){
    return (1./12.*gb_I(r, z-2*h, gb) - 2./3.*gb_I(r, z-h, gb)
    + 2./3.*gb_I(r, z+h, gb) - 1./12.*gb_I(r, z+2*h, gb))/h;
}

// Field square gradient in the r direction
double grad_E2_r(double r, double z, const gaussian_beam &gb, double h){
    return grad_I_r(r,z,gb,h)/(gb.c*gb.eps_0);
}

// Field square gradient in the z direction
double grad_E2_z(double r, double z, const gaussian_beam &gb, double h){
    return grad_I_z(r,z,gb,h)/(gb.c*gb.eps_0);
}

// Dipole force f_r(r,z), Gaussian beam 
// for a Rayleigh particle (a << lambda) or atom
double force_r_gb(double r, double z, double alpha, const gaussian_beam &gb, double h){
    return 1./2.*alpha*grad_E2_r(r,z,gb,h);
}

// Dipole force f_z(r,z), Gaussian beam 
// for a Rayleigh particle (a << lambda) or atom
double force_z_gb(double r, double z, double alpha, const gaussian_beam &gb, double h){
    return 1./2.*alpha*grad_E2_z(r,z,gb,h);
}

// paul trap force field
double force_paul_trap(double x, double t, const paul_trap &pt){
    return (pt.eps*cos(pt.w_dr*t)*x);
}

// Dipole force f_r(r,z), takes parameters from eq_params
double force_r(double r, double z, const eq_params &eq){
    return force_r_gb(r, z, eq.part.alpha, eq.gb, eq.part.h);
}

// Dipole force f_z(r,z), takes parameters from eq_params
double force_z(double r, double z, const eq_params &eq){
    return force_z_gb(r, z, eq.part.alpha, eq.gb, eq.part.h);
}


























