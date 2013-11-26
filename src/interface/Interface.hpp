#ifndef INTERFACE_HPP
#define INTERFACE_HPP
/*

  Interface functions prototypes.

*/
#include <vector>

#include "Config.hpp"

// Disable obnoxious warnings from Eigen headers
#if defined (__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wall"
#include <Eigen/Dense>
#pragma GCC diagnostic pop
#elif (__INTEL_COMPILER)
#pragma warning push
#pragma warning disable "-Wall"
#include <Eigen/Dense>
#pragma warning pop
#endif

#include "Sphere.hpp"
#include "WaveletCavity.hpp"

/*

	Functions visible to host program 

*/

extern "C" void hello_pcm_(int * a, double * b);

extern "C" void init_pcm_();

extern "C" void tear_down_pcm_();

extern "C" void comp_chg_pcm_(char* potString, char* chgString);

extern "C" void comp_pol_ene_pcm_(double * energy);

extern "C" void dot_surface_functions_(double * result, const char * potString, const char * chgString);

extern "C" void collect_nctot_(int * nuclei);

extern "C" void collect_atoms_(double * charges, double * centers);

extern "C" void get_cavity_size_(int * nts);

extern "C" void get_tess_centers_(double * centers);

extern "C" void get_tess_cent_coord_(int * its, double * center);

extern "C" void print_pcm_();

extern "C" void print_gepol_cavity_();

extern "C" void set_surface_function_(int * nts, double * values, char * name);

extern "C" void get_surface_function_(int * nts, double * values, char * name);

extern "C" void add_surface_function_(char * result, double * coeff, char * part);

extern "C" void print_surface_function_(char * name);

extern "C" void clear_surf_func_(char * name);

extern "C" void append_surf_func_(char * name);

/*

	Functions not visible to host program	

 */

void setupInput();

void initCavity(); 

void initSolver();

WaveletCavity * initWaveletCavity();

void init_pwcsolver_();

void init_pwlsolver_();

void initAtoms(Eigen::VectorXd & charges_, Eigen::Matrix3Xd & sphereCenter_);

void initSpheresImplicit(const Eigen::VectorXd & charges_, const Eigen::Matrix3Xd & sphereCenter_, std::vector<Sphere> & spheres_);

void initSpheresAtoms(const Eigen::VectorXd & charges_, const Eigen::Matrix3Xd & sphereCenter_, std::vector<Sphere> & spheres_);

bool surfaceFunctionExists(const std::string & name);

template<typename T> void safe_delete( T *& ptr ); 

#endif // INTERFACE_HPP

