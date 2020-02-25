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
 *  4.12.2019
 *********************************************************************** 
 */
 
 // STANDARD LIBRARIES
#include <cmath> // most math functions
#include <vector> // vector library

// MY LIBRARIES
#include "eq_definitions.hpp" 

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
    // x' = v
    f[0] =  v;
    
    //~ // IRENE & GIACOMO
    //~ double T1 = y[0];
    //~ double T2 = y[1];
    
    //~ double T0 = 295;
    //~ double l_ov_T1 = 1;
    //~ double l_ov_T2 = 2;
    //~ double l_res = l_0 + l_ov_T1*(T1 - T0) + l_ov_T2*(T2 - T0);
    
    // Definitions: experimental parameters

    // PAUL TRAP
    // f[1] = -eq.th.g_norm*v + eq.pt.eps*cos(eq.pt.w_dr*t)*x/eq.part.m; 
    // f[2] = x*v;
    
    // OTHER EXAMPLES:
    //
    // HARMONIC OSCILLATOR
    //f[1] = -eq.th.g_norm*v - pow(eq.ot.w, 2)*x;
    
    // OPT. TWEEZER TRAP
    f[1] = -eq.th.g_norm*v + force_z(0, x, eq)/eq.part.m;

    // HYBRID TRAP
    //f[1] = -eq.th.g_norm*v + eq.pt.eps*cos(eq.pt.w_dr*t)*x/eq.part.m + force_r(x, 0, eq)/eq.part.m;
   
    
    // MARC & JAN 
    // double w_0 = 2*M_PI*120e3;
    // double w_d = 2*M_PI*240e3;
    // double xi = -10e12;
    // double eta = 2e12; 
    // double eps = 0.02;
    // f[1] = -eq.th.g_norm*v - pow(w_0,2)*(1 + eps*cos(w_d*t) + eta*x*v/w_0 + xi*pow(x,2))*x;
        
    // return output
    return f;
}

// diffusion of process: b(t, Y)
std::vector<double> diffusion_function(std::vector<double> y, double t, const eq_params &eq){
    // output vector
    std::vector<double> f(y.size(), 0);
    // x = y[0], v = y[1]
    
    // Calculate force
    // No noise in x:
    f[0] = 0;
    // Stochastic force in momentum:
    f[1] = eq.th.sigma/eq.part.m;
    f[2] = 0;
    
    // return output
    return f;
}
