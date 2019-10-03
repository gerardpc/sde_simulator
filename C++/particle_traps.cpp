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
 *  03.10.2019
 *      Not yet working.
 *********************************************************************** 
 */

#include <iostream> // stdin, stdout 
#include <cstdio> // printf family
#include <cmath> // most math functions
#include <vector> // vector library
#include <string> // string library
#include <stdlib.h> // includes rand()

//======================================================================
// STRUCT DEFINITIONS HERE
//======================================================================
// Gaussian beam struct
struct gaussian_beam{
    double w_0;
    double z_R;
    double I_0;
};

// Paul trap struct
struct paul_trap{
    // particle parameters
    double r;
    double density;
    double m;
    // trap parameters
    double f_dr;
    double w_dr;
    double V;
    double d;
    double eps;
    // ambient parameters
    double T; 
    double pressure;
    double gamma_air;
    double gamma;
    double g_norm;
    double sigma;
};

// Optical tweezer struct
struct opt_tweezer{
    // particle parameters
    double r;
    double density;
    double m;
    // trap parameters
    double freq;
    double w;
    double T; 
    // ambient parameters
    double pressure;
    double gamma_air;
    double gamma;
    double g_norm;
    double sigma;
};

//======================================================================
// FUNCTIONS HERE
//======================================================================
// calculate particle mass of an assumed spherical silica particle
// r in meters
double particle_mass(double r){
    double silica_dens = 2200; // kg/m^3
    return 4./3.*M_PI*pow(r, 3)*silica_dens;
}

// Gaussian beam width
double gb_w(double z, gaussian_beam gb){
    return gb.w_0*sqrt(1 + pow((z/gb.z_R), 2));
}

// Gaussiain beam intensity
double gb_I(double r, double z, gaussian_beam gb){
    return gb.I_0*pow((gb.w_0/gb_w(z, gb)),2)*exp(-(2*pow(r/gb_w(z, gb), 2)));
}

// Intensity gradient in the r direction
double grad_I_r(double r, double z, gaussian_beam gb, double h){
    return (1./12.*gb_I(r-2*h, z, gb) - 2./3.*gb_I(r-h, z, gb) 
    + 2./3.*gb_I(r+h, z, gb) - 1./12.*gb_I(r+2*h, z, gb))/h;
}

// Intensity gradient in the z direction
double grad_I_z(double r, double z, gaussian_beam gb, double h){
    return (1./12.*gb_I(r, z-2*h, gb) - 2./3.*gb_I(r, z-h, gb)
    + 2./3.*gb_I(r, z+h, gb) - 1./12.*gb_I(r, z+2*h, gb))/h;
}

// Dipole force, Gaussian beam on r direction
double force_r_gb(double r, double z, double alpha, gaussian_beam gb, double h){
    return 1./2.*alpha*grad_I_r(r,z,gb,h);
}

// Dipole force, Gaussian beam on z direction
double force_z_gb(double r, double z, double alpha, gaussian_beam gb, double h){
    return 1./2.*alpha*grad_I_z(r,z,gb,h);
}

// Fill a Gaussian beam
int fill_gaussian_beam(gaussian_beam gb){
}

// Fill an optical tweezer
int fill_optical_tweezer(){
}

// Fill a Paul trap
int fill_paul_trap(){
}

// Calculate vector force for a Gaussian Beam dipole optical trap
// F = (f_r, f_z)
std::vector<double> force_gaussian_beam(double r, double z){
    // General parameters (i.e. constants)
    double length_unit = 1e-9; // nm
    double h_bar = 1.0546e-34;
    double c = 3e8;
    double eps_0 = 8.854e-12;
    double kB;
    
    // laser beam
    double lambda = 1064*unit; // laser wavelength 
    double P_0 = 70e-3; // laser power in W
    double NA = 0.8; // Gaussian beam NA
    
    // particle parameters
    double n = 1.43; // refractive index of particle
    double a = 150./2.*unit; // particle radius
    double m_particle = particle_mass(a); // particle mass
    double alpha = 4*M_PI*pow(a, 3)*eps_0*(pow(n, 2) - 1)/(pow(n, 2) + 2); // induced bulk dipole
    
    // numerical details
    double h = 1e-5*a; // discretization size
    double c_factor = 1.5; // correction factor for high NA in Gaussian beam (set to "1" to eliminate correction)
    
    // Gaussian beam
    gaussian_beam gb;
    gb.w_0 = lambda/(M_PI*NA)*c_factor;
    gb.z_R = M_PI*pow(gb.w_0, 2)/lambda;
    gb.I_0 = P_0*2/(M_PI*pow(gb.w_0, 2));

    // dipole force (on a Rayleigh particle (a << lambda) or atom)
    return {force_r_gb(r, z, alpha, gb, h), force_z_gb(r, z, alpha, gb, h)};
}

double force_paul_trap(double x){
    return (-trap.gamma*v + trap.eps*cos(trap.w_0*t)*x);
}

int main(){
    return 0;
}
