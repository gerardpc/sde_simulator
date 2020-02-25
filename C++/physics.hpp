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
    double k; 
    double P_0;
    double NA;
    // proper Gaussian beam parameters
    double w_0;
    double z_R;
    double I_0;
    void print(); // print values
    void fill(); // constructor
};

// Particle parameters struct
struct particle{
    double r;
    double density;
    double m;
    double n;
    double alpha;
    double alpha_im;
    double Q;
    double h;
    void print();
    void fill(gaussian_beam &gb);
};

// Thermodynamics parameters struct
struct thermodynamics{
    double k_B = 1.38065e-23;
    double T; 
    double pressure;
    double gamma_ambient;
    double gamma;
    double g_norm;
    double sigma;
    void print();
    void fill(particle &part);
};

// Optical tweezer struct
struct opt_trap{
    // freq. from file
    double freq;
    double w;
    // freq. derived from gaussian beam
    double f_gb_r;
    double w_gb_r;
    double f_gb_z;
    double w_gb_z;
    void print(); // print values
    void fill(); // constructor
    void fill_gb_w(particle &part, gaussian_beam &gb); // calculate frequencies from GB
};

// Paul trap struct
struct paul_trap{
    // trap parameters
    double f_dr;
    double w_dr;
    double V;
    double d;
    // simplified parameters
    double eps;
    double alpha_iz;
    double beta_iz;
    void print(); // print values
    void fill(particle &part, thermodynamics &th); // constructor
};

//======================================================================
// EQUATION PARAMETERS STRUCT
//======================================================================
struct eq_params{
    // the members of this struct are specific to my problem of interest,
    // but can be replaced to any other members/parameters to send variables
    // inside the drift and diffusion functions
    gaussian_beam gb;
    particle part;
    thermodynamics th;
    opt_trap ot;
    paul_trap pt;
    void fill();
    void print();
};

//======================================================================
// FIELD FUNCTIONS
//======================================================================
// calculate particle mass of an assumed spherical silica particle
// r in meters, density in kg/m^3
double particle_mass(double r, double density);

// Gaussian beam width
double gb_w(double z, const gaussian_beam &gb);

// Gaussian beam intensity
double gb_I(double r, double z, const gaussian_beam &gb);

// Intensity gradient in the r direction
double grad_I_r(double r, double z, const gaussian_beam &gb, double h);

// Intensity gradient in the z direction
double grad_I_z(double r, double z, const gaussian_beam &gb, double h);
    
// Field square gradient in the r direction
double grad_E2_r(double r, double z, const gaussian_beam &gb, double h);

// Field square gradient in the z direction
double grad_E2_z(double r, double z, const gaussian_beam &gb, double h);

// Gradient force f_r(r,z), Gaussian beam 
// for a Rayleigh particle (a << lambda) or atom
double force_r_gb(double r, double z, double alpha, const gaussian_beam &gb, double h);

// Gradient force f_z(r,z), Gaussian beam 
// for a Rayleigh particle (a << lambda) or atom
double force_z_gb(double r, double z, double alpha, const gaussian_beam &gb, double h);

// Scattering force f_scat(r,z), takes parameters from eq_params
double scat_force_z(double r, double z, const eq_params &eq);

//======================================================================
// FULL FORCES
//======================================================================
// OPTICAL FORCE:
// Dipole force f_r(r,z), takes parameters from eq_params
double force_r(double r, double z, const eq_params &eq);

// Dipole force f_z(r,z), takes parameters from eq_params
double force_z(double r, double z, const eq_params &eq);

// PAUL TRAP FORCE:
// paul trap force field
double force_paul_trap(double x, double t, const paul_trap &pt);




















