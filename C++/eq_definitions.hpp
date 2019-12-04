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
#include <vector> // vector library

// MY LIBRARIES
#include "physics.hpp" // optical tweezer, paul trap and gaussian beam

//======================================================================
// PROBLEM FUNCTIONS
//======================================================================
// SDE FUNCTIONS
// dY = a(t,Y)dt + b(t,Y)*dW_t, Y(t0) = Y0,

// drift of process: a(t,Y)
std::vector<double> drift_function(std::vector<double> y, double t, const eq_params &args);

// diffusion of process: b(t, Y)
std::vector<double> diffusion_function(std::vector<double> y, double t, const eq_params &args);
