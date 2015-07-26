/* pcmsolver_copyright_start */
/*
 *     PCMSolver, an API for the Polarizable Continuum Model
 *     Copyright (C) 2013 Roberto Di Remigio, Luca Frediani and contributors
 *
 *     This file is part of PCMSolver.
 *
 *     PCMSolver is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU Lesser General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *
 *     PCMSolver is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU Lesser General Public License for more details.
 *
 *     You should have received a copy of the GNU Lesser General Public License
 *     along with PCMSolver.  If not, see <http://www.gnu.org/licenses/>.
 *
 *     For information on the complete list of contributors to the
 *     PCMSolver API, see: <http://pcmsolver.github.io/pcmsolver-doc>
 */
/* pcmsolver_copyright_end */

#ifndef INTERFACE_HPP
#define INTERFACE_HPP
/*

  Interface functions prototypes.

*/

#include <map>
#include <vector>

#include "Config.hpp"

#include "FCMangle.hpp"

#include <Eigen/Core>

#include "InputManager.hpp"
#include "Molecule.hpp"
#include "Sphere.hpp"
#include "SurfaceFunction.hpp"

/*

	Functions visible to host program

*/

#ifdef __cplusplus
extern "C" {
#endif

#define hello_pcm \
	FortranCInterface_GLOBAL_(hello_pcm, HELLO_PCM)
void hello_pcm(int * a, double * b);

#define set_up_pcm \
	FortranCInterface_GLOBAL_(set_up_pcm, SET_UP_PCM)
/*! \fn void set_up_pcm(int * host_provides_input)
 *  \param[in] host_provides_input whether the host does syntactic parsing for API
 */
void set_up_pcm(int * host_provides_input);

#define tear_down_pcm \
	FortranCInterface_GLOBAL_(tear_down_pcm, TEAR_DOWN_PCM)
void tear_down_pcm();

#define write_timings \
	FortranCInterface_GLOBAL_(write_timings, WRITE_TIMINGS)
/*! \fn void write_timings()
 *  \brief Write API timings to pcmsolver.timer.dat
 */
void write_timings();

#define compute_asc \
	FortranCInterface_GLOBAL_(compute_asc, COMPUTE_ASC)
/*! \fn void compute_asc(char * potString, char * chgString, int * irrep)
 *  \brief Computes the ASC corresponding to the MEP and irreducible representation, using the equilibrium PCM matrix
 *  \param[in] potString name of the potential SurfaceFunction
 *  \param[in] chgString name of the charge SurfaceFunction
 *  \param[in] irrep     the irreducible representation the potential/charge pair belongs to
 *
 *  The equilibrium PCM matrix is used in the computation of the ASC
 */
void compute_asc(char * potString, char * chgString, int * irrep);

#define compute_nonequilibrium_asc \
	FortranCInterface_GLOBAL_(compute_nonequilibrium_asc, COMPUTE_NONEQUILIBRIUM_ASC)
/*! \fn void compute_nonequilibrium_asc(char * potString, char * chgString, int * irrep)
 *  \brief Computes the ASC correspoding to the MEP and irreducible representation, using the nonequilibrium PCM matrix
 *  \param[in] potString name of the potential SurfaceFunction
 *  \param[in] chgString name of the charge SurfaceFunction
 *  \param[in] irrep     the irreducible representation the potential/charge pair belongs to
 *
 *  The nonequilibrium PCM matrix is used in the computation of the ASC
 */
void compute_nonequilibrium_asc(char * potString, char * chgString, int * irrep);

#define compute_polarization_energy \
	FortranCInterface_GLOBAL_(compute_polarization_energy, COMPUTE_POLARIZATION_ENERGY)
/*! \fn void compute_polarization_energy(double * energy)
 *  \param[out] energy polarization energy \f$ U_\mathrm{pol} = \frac{1}{2}\mathbf{q}\cdot\mathbf{v} \f$
 */
void compute_polarization_energy(double * energy);

#define save_surface_functions \
	FortranCInterface_GLOBAL_(save_surface_functions, SAVE_SURFACE_FUNCTIONS)
/*! \fn void save_surface_functions()
 *  \brief Dumps the surface functions in the SharedSurfaceFunctionMap functions to a .npy file
 *
 *  The surface functions are saved as NumPy arrays. This function should be called **before** the
 *  call to tear_down_pcm() so that the converged MEP and ASC are saved to file for a possible later
 *  use. For example, they might be needed in a response calculation.
 *  The filename is decided according to the name of the surface function.
 */
void save_surface_functions();

#define save_surface_function \
	FortranCInterface_GLOBAL_(save_surface_function, SAVE_SURFACE_FUNCTION)
/*! \fn void save_surface_function(const char * name)
 *  \brief Dumps the surface function to a .npy file
 *  \param[in] name the name of the SurfaceFunction to be saved to file
 */
void save_surface_function(const char * name);

#define load_surface_function \
	FortranCInterface_GLOBAL_(load_surface_function, LOAD_SURFACE_FUNCTION)
/*! \fn void load_surface_function(const char * name)
 *  \brief Loads the surface function from a .npy file
 *  \param[in] name the name of the surface function to be loaded from file
 *
 *  The name has to be consistent with that of the .npy file saved using save_surface_functions()
 */
void load_surface_function(const char * name);

#define dot_surface_functions \
	FortranCInterface_GLOBAL_(dot_surface_functions, DOT_SURFACE_FUNCTIONS)
/*! \fn void dot_surface_functions(double * result, const char * potString, const char * chgString)
 *  \param[out]   result dot product of the surface function pair
 *  \param[in] potString name of the potential SurfaceFunction
 *  \param[in] chgString name of the chareg SurfaceFunction
 */
void dot_surface_functions(double * result, const char * potString,
                                      const char * chgString);

/*! \fn void collect_nctot(int * nuclei)
 *  \param[out] nuclei number of nuclei in the molecule
 */
void collect_nctot(int * nuclei);

/*! \fn void collect_atoms(double * charges, double * centers)
 *  \param[out] charges atomic charge of the nuclei
 *  \param[out] centers coordinate of the nuclei (in column-major order)
 */
void collect_atoms(double * charges, double * centers);

/*! \fn void host_writer(const char * message, size_t * message_length)
 *  \param[in]        message message to be flushed to host program output
 *  \param[in] message_length length of the message
 */
void host_writer(const char * message, size_t * message_length);

/*! \fn void set_point_group(int * nr_generators, int * gen1, int * gen2, int * gen3)
 *  \param[out] nr_generators number of generators in the molecular point group
 *  \param[out] gen1 first group generator
 *  \param[out] gen2 second group generator
 *  \param[out] gen3 third group generator
 */
void set_point_group(int * nr_generators, int * gen1, int * gen2,
                                int * gen3);

/*! \fn void host_input(cavityInput * cav, solverInput * solv, greenInput * green)
 *  \brief Pushes input parameters defined host-side to the API
 *  \param[in] cav   cavity input data
 *  \param[in] solv  solver input data
 *  \param[in] green Green's function input data
 */
void host_input(cavityInput * cav, solverInput * solv, greenInput * green);

#define push_input_string \
	FortranCInterface_GLOBAL_(push_input_string, PUSH_INPUT_STRING)
/*! \fn void push_input_string(char * s)
 *  \brief pushes a string input parameter to the API
 *  \param[in] s the string to be passed
 *
 *  The string is saved inside a std::vector<std::string>
 *  This function is to be used in conjunction with input management done
 *  host-side.
 */
void push_input_string(char * s);

#define get_cavity_size \
	FortranCInterface_GLOBAL_(get_cavity_size, GET_CAVITY_SIZE)
/*! \fn void get_cavity_size(int * nts, int * ntsirr)
 *  \param nts    the total size of the cavity
 *  \param ntsirr the irreducible size of the cavity
 */
void get_cavity_size(int * nts, int * ntsirr);

#define get_tesserae \
	FortranCInterface_GLOBAL_(get_tesserae, GET_TESSERAE)
/*! \fn void get_tesserae(double * centers)
 *  \param centers the grid of points provided by discretization of the cavity
 */
void get_tesserae(double * centers);

#define get_tesserae_centers \
	FortranCInterface_GLOBAL_(get_tesserae_centers, GET_TESSERAE_CENTERS)
/*! \fn void get_tesserae_centers(int * its, double * center)
 *  \param[in] its the index of the tesserae
 *  \param  center the coordinates of the its-th tessera
 */
void get_tesserae_centers(int * its, double * center);

#define print_citation \
    FortranCInterface_GLOBAL_(print_citation, PRINT_CITATION)
void print_citation();

#define print_pcm \
	FortranCInterface_GLOBAL_(print_pcm, PRINT_PCM)
void print_pcm();

#define set_surface_function \
	FortranCInterface_GLOBAL_(set_surface_function, SET_SURFACE_FUNCTION)
/*! \fn void set_surface_function(int * nts, double * values, char * name)
 *  \brief Wraps a raw array of doubles into a SurfaceFunction with the given label
 *  \param[in]    nts length of the values vector
 *  \param[in] values raw array containing the values of a function at the cavity surface
 *  \param[in]   name label for the SurfaceFunction
 */
void set_surface_function(int * nts, double * values, char * name);

#define get_surface_function \
	FortranCInterface_GLOBAL_(get_surface_function, GET_SURFACE_FUNCTION)
/*! \fn void get_surface_function(int * nts, double * values, char * name)
 *  \brief Retrieves a raw array of doubles from the SurfaceFunction with the given label
 *  \param[in]     nts length of the values vector
 *  \param[out] values raw array containing the values of a function at the cavity surface
 *  \param[in]    name label for the SurfaceFunction
 */
void get_surface_function(int * nts, double * values, char * name);

#define add_surface_function \
	FortranCInterface_GLOBAL_(add_surface_function, ADD_SURFACE_FUNCTION)
/*! \fn void add_surface_function(char * result, double * coeff, char * part)
 *  \brief DAXPY-like utility for SurfaceFunction \f$ f_\mathrm{result} := f_\mathrm{result} + c f_\mathrm{part} \f$
 *  \param[in] result label for the resulting SurfaceFunction
 *  \param[in]  coeff coefficient
 *  \param[in]   part label for the LHS SurfaceFunction
 */
void add_surface_function(char * result, double * coeff, char * part);

#define print_surface_function \
	FortranCInterface_GLOBAL_(print_surface_function, PRINT_SURFACE_FUNCTION)
/*! \fn void print_surface_function(char * name)
 *  \param[in] name label for the SurfaceFunction to be printed
 */
void print_surface_function(char * name);

#define clear_surface_function \
	FortranCInterface_GLOBAL_(clear_surface_function, CLEAR_SURFACE_FUNCTION)
/*! \fn void clear_surface_function(char * name)
 *  \brief Clears contents of SurfaceFunction
 *  \param[in] name label for the SurfaceFunction to be cleared
 */
void clear_surface_function(char * name);

#define append_surface_function \
	FortranCInterface_GLOBAL_(append_surface_function, APPEND_SURFACE_FUNCTION)
/*! \fn void append_surface_function(char * name)
 *  \brief Appends a new SurfaceFunction to the global map
 *  \param[in] name label for the new SurfaceFunction
 */
void append_surface_function(char * name);

/*! \fn void scale_surface_function(char * func, double * coeff)
 *  \param[in] func  the name of the SurfaceFunction to be scaled
 *  \param[in] coeff the scaling coefficient
 */
#define scale_surface_function \
	FortranCInterface_GLOBAL_(scale_surface_function, SCALE_SURFACE_FUNCTION)
void scale_surface_function(char * func, double * coeff);

#ifdef __cplusplus
}
#endif

/*

	Functions not visible to host program

 */

/*! \fn void setupInput()
 *
 *  Sets up Input object
 *  \param[in] from_host whether input parameters are passed from host
 */
void setupInput(bool from_host);

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

/*! \fn void initNonEqSolver()
 *
 *  Creates PCMSolver object, nonequilibrium version
 */
void initNonEqSolver();

/*! \fn void initAtoms(Eigen::VectorXd & charges_, Eigen::Matrix3Xd & sphereCenter_)
 *  \param[in] charges_      contains charges of atomic centers
 *  \param[in] sphereCenter_ contains coordinates of atomic centers
 *
 *  Initializes list of atoms and atoms positions
 */
void initAtoms(Eigen::VectorXd & charges_, Eigen::Matrix3Xd & sphereCenter_);

/*! \fn void initMolecule(Molecule & molecule_, Eigen::Matrix3Xd & sphereCenter_)
 *  \param[out] molecule_  a Molecule object
 *
 *  Initializes Molecule object with all the information on charges, masses and positions of atomic centers
 */
void initMolecule(Molecule & molecule_);

/*! \fn void initSpheresImplicit(const Eigen::VectorXd & charges_, const Eigen::Matrix3Xd & sphereCenter_, std::vector<Sphere> & spheres_)
 *  \param[in] charges_      contains charges of atomic centers
 *  \param[in] sphereCenter_ contains coordinates of atomic centers
 *  \param[in] spheres_      contains list of spheres
 *
 *  Generates the list of spheres needed to build the cavity from list of atoms and atoms centers
 *  Uses one of the predefined sets of radii
 */
void initSpheresImplicit(const Eigen::VectorXd & charges_,
                         const Eigen::Matrix3Xd & sphereCenter_, std::vector<Sphere> & spheres_);

/*! \fn void initSpheresAtoms(const Eigen::Matrix3Xd & sphereCenter_, std::vector<Sphere> & spheres_)
 *  \param[in] sphereCenter_ contains coordinates of atomic centers
 *  \param[in] spheres_      contains list of spheres
 *
 *  Generates the list of spheres needed to build the cavity from list of atoms centers
 *  Substitures radii on specified atoms with custom values from input
 */
void initSpheresAtoms(const Eigen::Matrix3Xd & sphereCenter_,
                      std::vector<Sphere> & spheres_);

/*! \fn bool surfaceFunctionExists(const std::string & name)
 *  \param[in] sf_map the map to be searched
 *  \param[in] name name of the SurfaceFunction
 *
 *  Checks if SurfaceFunction exists
 */
bool surfaceFunctionExists(const std::map<std::string, SurfaceFunction> & sf_map, const std::string & name);

/*! \fn void insertSurfaceFunction(std::map<std::string, SurfaceFunction> & sf_map, const SurfaceFunction & sf)
 *  \param[in] sf_map the map to be searched
 *  \param[in] sf the SurfaceFunction
 *
 *  No checks on existence are performed by this function
 */
void insertSurfaceFunction(std::map<std::string, SurfaceFunction> & sf_map, const SurfaceFunction & sf);

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

