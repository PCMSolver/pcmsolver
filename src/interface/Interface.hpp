#ifndef INTERFACE_HPP
#define INTERFACE_HPP
/*

  Interface functions prototypes.

*/

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
/*! \fn extern "C" void compute_asc(char * potString, char * chgString, int * irrep)
 *  \param[in] potString name of the potential SurfaceFunction
 *  \param[in] chgString name of the charge SurfaceFunction
 *  \param[in] irrep     the irreducible representation the potential/charge pair belongs to
 */
extern "C" void compute_asc(char * potString, char * chgString, int * irrep);

#define compute_polarization_energy \
	FortranCInterface_GLOBAL_(compute_polarization_energy, COMPUTE_POLARIZATION_ENERGY)
extern "C" void compute_polarization_energy(double * energy);

#define dot_surface_functions \
	FortranCInterface_GLOBAL_(dot_surface_functions, DOT_SURFACE_FUNCTIONS)
extern "C" void dot_surface_functions(double * result, const char * potString, const char * chgString);

extern "C" void collect_nctot(int * nuclei);

extern "C" void collect_atoms(double * charges, double * centers);

extern "C" void host_writer(const char * message, size_t * message_length);

extern "C" void set_point_group(int * pg);

#define get_cavity_size \
	FortranCInterface_GLOBAL_(get_cavity_size, GET_CAVITY_SIZE)
/*! \fn extern "C" void get_cavity_size(int * nts, int * ntsirr)
 *  \param nts    the total size of the cavity
 *  \param ntsirr the irreducible size of the cavity
 */
extern "C" void get_cavity_size(int * nts, int * ntsirr);

#define get_tesserae \
	FortranCInterface_GLOBAL_(get_tesserae, GET_TESSERAE)
/*! \fn extern "C" void get_tesserae(double * centers)
 *  \param centers the grid of points provided by discretization of the cavity
 */
extern "C" void get_tesserae(double * centers);

#define get_tesserae_centers \
	FortranCInterface_GLOBAL_(get_tesserae_centers, GET_TESSERAE_CENTERS)
/*! \fn extern "C" void get_tesserae_centers(int * its, double * center)
 *  \param[in] its the index of the tesserae
 *  \param  center the coordinates of the its-th tessera
 */
extern "C" void get_tesserae_centers(int * its, double * center);

#define print_citation \
    FortranCInterface_GLOBAL_(print_citation, PRINT_CITATION)
extern "C" void print_citation();

#define print_pcm \
	FortranCInterface_GLOBAL_(print_pcm, PRINT_PCM)
extern "C" void print_pcm();

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

/*! \fn extern "C" void scale_surface_function(char * func, double * coeff)
 *  \param[in] func  the name of the SurfaceFunction to be scaled
 *  \param[in] coeff the scaling coefficient
 */
#define scale_surface_function \
	FortranCInterface_GLOBAL_(scale_surface_function, SCALE_SURFACE_FUNCTION)
extern "C" void scale_surface_function(char * func, double * coeff);

/*

	Functions not visible to host program	

 */

/*! \fn void setupInput()
 *
 *  Sets up Input object
 */
void setupInput();

/*! \fn void initCavity()
 *  
 *  Creates Cavity object
 */
void initCavity(); 

/*! \fn void initSolver()
 *
 *  Creates PCMSolver object
 */
void initSolver();

/*! \fn void initWaveletCavity()
 * 
 *  Initializes the _waveletCavity global object
 */
void initWaveletCavity();

/*! \fn void initAtoms(Eigen::VectorXd & charges_, Eigen::Matrix3Xd & sphereCenter_)
 *  \param[in] charges_      contains charges of atomic centers
 *  \param[in] sphereCenter_ contains coordinates of atomic centers
 *
 *  Initializes list of atoms and atoms positions
 */
void initAtoms(Eigen::VectorXd & charges_, Eigen::Matrix3Xd & sphereCenter_);

/*! \fn void initSpheresImplicit(const Eigen::VectorXd & charges_, const Eigen::Matrix3Xd & sphereCenter_, std::vector<Sphere> & spheres_)
 *  \param[in] charges_      contains charges of atomic centers
 *  \param[in] sphereCenter_ contains coordinates of atomic centers
 *  \param[in] spheres_      contains list of spheres
 *
 *  Generates the list of spheres needed to build the cavity from list of atoms and atoms centers
 *  Uses one of the predefined sets of radii
 */
void initSpheresImplicit(const Eigen::VectorXd & charges_, const Eigen::Matrix3Xd & sphereCenter_, std::vector<Sphere> & spheres_);

/*! \fn void initSpheresAtoms(const Eigen::Matrix3Xd & sphereCenter_, std::vector<Sphere> & spheres_)
 *  \param[in] sphereCenter_ contains coordinates of atomic centers
 *  \param[in] spheres_      contains list of spheres
 *
 *  Generates the list of spheres needed to build the cavity from list of atoms centers
 *  Substitures radii on specified atoms with custom values from input
 */
void initSpheresAtoms(const Eigen::Matrix3Xd & sphereCenter_, std::vector<Sphere> & spheres_);

/*! \fn bool surfaceFunctionExists(const std::string & name)
 *  \param[in] name name of the SurfaceFunction
 *
 *  Checks if SurfaceFunction exists
 */
bool surfaceFunctionExists(const std::string & name);

/*! \fn template<typename T> void safe_delete( T *& ptr );
 *  \param[in] ptr the pointer to be deleted
 *  \tparam    T   the type of the pointed-to object
 *
 *  Implements safe deletion of pointers
 */
template<typename T> 
void safe_delete( T *& ptr ); 

/*! \fn inline void printer(const std::string & message)
 *  \param[in] message the message to be printed
 *
 *  Flushes message to host program printing unit
 */
inline void printer(const std::string & message);

/*! \fn inline void printer(const std::ostringstream & stream)
 *  \param[in] stream the message to be printed
 *
 *  Flushes stream to host program printing unit
 */
inline void printer(std::ostringstream & stream);

#endif // INTERFACE_HPP

