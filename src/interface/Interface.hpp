#ifndef INTERFACE_HPP
#define INTERFACE_HPP
/*

  Interface functions prototypes.

*/
#include <vector>

#include "Config.hpp"
#include "FCMangle.hpp"

// Disable obnoxious warnings from Eigen headers
#if defined (__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wall" 
#pragma GCC diagnostic ignored "-Weffc++" 
#pragma GCC diagnostic ignored "-Wextra"
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

#define hello_pcm \
	FortranCInterface_GLOBAL_(hello_pcm, HELLO_PCM)
extern "C" void hello_pcm(int * a, double * b);

#define init_pcm \
	FortranCInterface_GLOBAL_(init_pcm, INIT_PCM)
extern "C" void init_pcm();

#define tear_down_pcm \
	FortranCInterface_GLOBAL_(tear_down_pcm, TEAR_DOWN_PCM)
extern "C" void tear_down_pcm();

#define comp_chg_pcm \
	FortranCInterface_GLOBAL_(comp_chg_pcm, COMP_CHG_PCM)
extern "C" void comp_chg_pcm(char* potString, char* chgString);

#define comp_pol_ene_pcm \
	FortranCInterface_GLOBAL_(comp_pol_ene_pcm, COMP_POL_ENE_PCM)
extern "C" void comp_pol_ene_pcm(double * energy);

#define dot_surface_functions \
	FortranCInterface_GLOBAL_(dot_surface_functions, DOT_SURFACE_FUNCTIONS)
extern "C" void dot_surface_functions(double * result, const char * potString, const char * chgString);

extern "C" void collect_nctot(int * nuclei);

extern "C" void collect_atoms(double * charges, double * centers);

#define get_cavity_size \
	FortranCInterface_GLOBAL_(get_cavity_size, GET_CAVITY_SIZE)
extern "C" void get_cavity_size(int * nts);

#define get_tess_centers \
	FortranCInterface_GLOBAL_(get_tess_centers, GET_TESS_CENTERS)
extern "C" void get_tess_centers(double * centers);

#define get_tess_cent_coord \
	FortranCInterface_GLOBAL_(get_tess_cent_coord, GET_TESS_CENT_COORD)
extern "C" void get_tess_cent_coord(int * its, double * center);

#define print_pcm \
	FortranCInterface_GLOBAL_(print_pcm, PRINT_PCM)
extern "C" void print_pcm();

#define print_gepol_cavity \
	FortranCInterface_GLOBAL_(print_gepol_cavity, PRINT_GEPOL_CAVITY)
extern "C" void print_gepol_cavity();

#define set_surface_function \
	FortranCInterface_GLOBAL_(set_surface_function, SET_SURFACE_FUNCTION)
extern "C" void set_surface_function(int * nts, double * values, char * name);

#define get_surface_function \
	FortranCInterface_GLOBAL_(get_surface_function, GET_SURFACE_FUNCTION)
extern "C" void get_surface_function(int * nts, double * values, char * name);

#define add_surface_function \
	FortranCInterface_GLOBAL_(add_surface_function, ADD_SURFACE_FUNCTION)
extern "C" void add_surface_function(char * result, double * coeff, char * part);

#define print_surface_function \
	FortranCInterface_GLOBAL_(print_surface_function, PRINT_SURFACE_FUNCTION)
extern "C" void print_surface_function(char * name);

#define clear_surf_func \
	FortranCInterface_GLOBAL_(clear_surf_func, CLEAR_SURF_FUNC)
extern "C" void clear_surf_func(char * name);

#define append_surf_func \
	FortranCInterface_GLOBAL_(append_surf_func, APPEND_SURF_FUNC)
extern "C" void append_surf_func(char * name);

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

