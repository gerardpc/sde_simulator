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
 *  08.10.2019
 *      Parts are working, field functions not yet tested.
 *********************************************************************** 
 */
 
// STANDARD LIBRARIES
#include <vector> // vector library

//======================================================================
// STRUCT DEFINITIONS HERE
//
// Stack experimental parameters together and
// pass them between functions in a more elegant way
//======================================================================
// Gaussian beam struct
struct gaussian_beam{
    // Physical constants
    double length_unit = 1e-9; // nm
    double h_bar = 1.0546e-34; // Planck's constant
    double c = 3e8;  // speed of light
    double eps_0 = 8.854e-12; // vacuum diel. permittivity
    // laser/lens parameters
    double c_factor;
    double lambda;
    double P_0;
    double NA;
    // proper Gaussian beam parameters
    double w_0;
    double z_R;
    double I_0;
    void print(); // print values
    void fill(); // constructor
};

// Paul trap struct
struct paul_trap{
    // particle parameters
    double r;
    double density;
    double m;
    double n;
    double alpha;
    // trap parameters
    double f_dr;
    double w_dr;
    double Q;
    double V;
    double d;
    double eps;
    double alpha_iz;
    double beta_iz;
    // ambient parameters
    double T; 
    double pressure;
    double gamma_ambient;
    double gamma;
    double g_norm;
    double sigma;
    void print(); // print values
    void fill(); // constructor
};

// Optical tweezer struct
struct opt_tweezer{
    // particle parameters
    double r;
    double density;
    double m;
    double n;
    double alpha;
    // trap parameters from file
    double freq;
    double w;
    // trap parameters from gaussian beam
    double f_gb_r;
    double w_gb_r;
    double f_gb_z;
    double w_gb_z;
    // ambient parameters
    double T; 
    double pressure;
    double gamma_ambient;
    double gamma;
    double g_norm;
    double sigma;
    void print(); // print values
    void fill(); // constructor
    void fill_gb_w(gaussian_beam &gb); // calculate frequencies from GB
};

//======================================================================
// FIELD FUNCTIONS
//======================================================================
// calculate particle mass of an assumed spherical silica particle
// r in meters, density in kg/m^3
double particle_mass(double r, double density);

// Gaussian beam width
double gb_w(double z, const gaussian_beam &gb);

// Gaussiain beam intensity
double gb_I(double r, double z, const gaussian_beam &gb);

// Intensity gradient in the r direction
double grad_I_r(double r, double z, const gaussian_beam &gb, double h);

// Intensity gradient in the z direction
double grad_I_z(double r, double z, const gaussian_beam &gb, double h);

// Field square gradient in the r direction
double grad_E2_r(double r, double z, const gaussian_beam &gb, double h);

// Field square gradient in the z direction
double grad_E2_z(double r, double z, const gaussian_beam &gb, double h);

// Dipole force f_r(r,z), Gaussian beam 
// for a Rayleigh particle (a << lambda) or atom
double force_r_gb(double r, double z, double alpha, const gaussian_beam &gb, double h);

// Dipole force f_z(r,z), Gaussian beam 
// for a Rayleigh particle (a << lambda) or atom
double force_z_gb(double r, double z, double alpha, const gaussian_beam &gb, double h);

// paul trap force field
double force_paul_trap(double x, double t, const paul_trap &pt);
























