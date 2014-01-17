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

#define set_up_pcm \
	FortranCInterface_GLOBAL_(set_up_pcm, SET_UP_PCM)
extern "C" void set_up_pcm();

#define tear_down_pcm \
	FortranCInterface_GLOBAL_(tear_down_pcm, TEAR_DOWN_PCM)
extern "C" void tear_down_pcm();

#define compute_asc \
	FortranCInterface_GLOBAL_(compute_asc, COMPUTE_ASC)
extern "C" void compute_asc(char * potString, char * chgString);

#define compute_polarization_energy \
	FortranCInterface_GLOBAL_(compute_polarization_energy, COMPUTE_POLARIZATION_ENERGY)
extern "C" void compute_polarization_energy(double * energy);

#define dot_surface_functions \
	FortranCInterface_GLOBAL_(dot_surface_functions, DOT_SURFACE_FUNCTIONS)
extern "C" void dot_surface_functions(double * result, const char * potString, const char * chgString);

extern "C" void collect_nctot(int * nuclei);

extern "C" void collect_atoms(double * charges, double * centers);

extern "C" void host_writer(const char * message, size_t * message_length);

#define get_cavity_size \
	FortranCInterface_GLOBAL_(get_cavity_size, GET_CAVITY_SIZE)
extern "C" void get_cavity_size(int * nts);

#define get_tesserae \
	FortranCInterface_GLOBAL_(get_tesserae, GET_TESSERAE)
extern "C" void get_tesserae(double * centers);

#define get_tesserae_centers \
	FortranCInterface_GLOBAL_(get_tesserae_centers, GET_TESSERAE_CENTERS)
extern "C" void get_tesserae_centers(int * its, double * center);

#define print_citation \
    FortranCInterface_GLOBAL_(print_citation, PRINT_CITATION)
extern "C" void print_citation();

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

#define clear_surface_function \
	FortranCInterface_GLOBAL_(clear_surface_function, CLEAR_SURFACE_FUNCTION)
extern "C" void clear_surface_function(char * name);

#define append_surface_function \
	FortranCInterface_GLOBAL_(append_surface_function, APPEND_SURFACE_FUNCTION)
extern "C" void append_surface_function(char * name);

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

void initSpheresAtoms(const Eigen::Matrix3Xd & sphereCenter_, std::vector<Sphere> & spheres_);

bool surfaceFunctionExists(const std::string & name);

template<typename T> 
void safe_delete( T *& ptr ); 

#endif // INTERFACE_HPP

