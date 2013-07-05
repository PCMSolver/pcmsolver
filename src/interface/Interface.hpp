#ifndef INTERFACE_H_
#define INTERFACE_H_
/*

  Interface functions prototypes.

*/
#include <vector>

#include "Config.hpp"

#include <Eigen/Dense>

#include "Sphere.hpp"
#include "WaveletCavity.hpp"

/*

	Functions visible to host program 

*/

extern "C" void hello_pcm_(int * a, double * b);

extern "C" void init_pcm_();

extern "C" void comp_chg_pcm_(char* potString, char* chgString);

extern "C" void comp_pol_ene_pcm_(double * energy, int * separate_or_total);

extern "C" void get_epsilon_static_(double * epsilon);

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

extern "C" bool surf_func_exists_(char * name);

extern "C" void clear_surf_func_(char * name);

extern "C" void append_surf_func_(char * name);

/*

	Functions not visible to host program	

 */

void setupInput();

// 1. Declare a global Cavity * _cavity; 2. use the factory inside here; 3. _cavity = _theCavityYouWant
void initCavity(); 

// 1. Declare a global PCMSolver * _solver; 2. use the factory inside here; 3. _solver = _theSolverYouWant
void initSolver(); // The GreensFunctionFactory will be used here to generate the inside & outside Green's Functions

WaveletCavity * initWaveletCavity();

void init_wavelet_cavity_();

void init_iefsolver_();

void init_cpcmsolver_();

void init_pwcsolver_();

void init_pwlsolver_();

void build_isotropic_matrix_();

void build_anisotropic_matrix_();

void initAtoms(Eigen::VectorXd & _charges, Eigen::Matrix3Xd & _sphereCenter);

std::vector<Sphere> initSpheresImplicit(const Eigen::VectorXd & _charges, const Eigen::Matrix3Xd & _sphereCenter);

std::vector<Sphere> initSpheresAtoms(const Eigen::VectorXd & _charges, const Eigen::Matrix3Xd & _sphereCenter);

void init_atoms_(Eigen::VectorXd & charges, 
                 Eigen::Matrix3Xd & sphereCenter);

void init_spheres_implicit_(Eigen::VectorXd & charges, 
                 Eigen::Matrix3Xd & centers);

void init_spheres_atoms_(Eigen::VectorXd & charges, 
                                    Eigen::Matrix3Xd & centers);


#endif

