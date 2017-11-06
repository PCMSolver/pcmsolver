/**
 * PCMSolver, an API for the Polarizable Continuum Model
 * Copyright (C) 2017 Roberto Di Remigio, Luca Frediani and collaborators.
 *
 * This file is part of PCMSolver.
 *
 * PCMSolver is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * PCMSolver is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with PCMSolver.  If not, see <http://www.gnu.org/licenses/>.
 *
 * For information on the complete list of contributors to the
 * PCMSolver API, see: <http://pcmsolver.readthedocs.io/>
 */

#pragma once

#include <stddef.h>

#include "PCMInput.h"
#include "PCMSolverExport.h"

#define PCMSolver_API PCMSolver_EXPORT

// To cope with the fact that C doesn't have bool as primitive type
#ifndef pcmsolver_bool_t_DEFINED
#define pcmsolver_bool_t_DEFINED
#if (defined(__STDC__) && (__STDC_VERSION__ < 199901L)) && !defined(__cplusplus)
typedef enum { pcmsolver_false, pcmsolver_true } pcmsolver_bool_t;
#else /* (defined(__STDC__) || (__STDC_VERSION__ < 199901L)) &&                     \
         !defined(__cplusplus) */
#include <stdbool.h>
typedef bool pcmsolver_bool_t;
#endif /* (defined(__STDC__) || (__STDC_VERSION__ < 199901L)) &&                    \
          !defined(__cplusplus) */
#endif /* pcmsolver_bool_t_DEFINED */

/*! \file pcmsolver.h
 *  \brief C API to PCMSolver
 *  \author Roberto Di Remigio
 *  \date 2015
 */

#ifdef __cplusplus
extern "C" {
#endif

/*! \struct pcmsolver_context_s
 *  Forward-declare opaque handle to a PCM context
 */
struct pcmsolver_context_s;
/*! \typedef pcmsolver_context_t
 *  Workaround to have pcmsolver_context_t available to C
 */
typedef struct pcmsolver_context_s pcmsolver_context_t;

/*! \struct PCMInput
 *  Forward-declare PCMInput input wrapping struct
 */
struct PCMInput;

/*! \enum pcmsolver_reader_t
 *  \brief Input processing strategies
 */
typedef enum {
  PCMSOLVER_READER_OWN, /*!< Module reads input on its own */
  PCMSOLVER_READER_HOST /*!< Module receives input from host */
} pcmsolver_reader_t;

/*! \typedef HostWriter
 *  Flushes module output to host program
 *  \param[in,out] message contents of the module output
 */
typedef void (*HostWriter)(const char * message);

/*! \brief Creates a new PCM context object
 *  \param[in] input_reading input processing strategy
 *  \param[in] nr_nuclei     number of atoms in the molecule
 *  \param[in] charges       atomic charges
 *  \param[in] coordinates   atomic coordinates
 *  \param[in] symmetry_info molecular point group information
 *  \param[in] host_input    input to the module, as read by the host
 *
 *  The molecular point group information is passed as an array
 *  of 4 integers: number of generators, first, second and third generator
 *  respectively. Generators map to integers as in table :ref:`symmetry-ops`
 */
PCMSolver_API pcmsolver_context_t * pcmsolver_new(pcmsolver_reader_t input_reading,
                                                  int nr_nuclei,
                                                  double charges[],
                                                  double coordinates[],
                                                  int symmetry_info[],
                                                  struct PCMInput * host_input,
                                                  HostWriter writer);

/*! \brief Deletes a PCM context object
 *  \param[in, out] context the PCM context object to be deleted
 */
PCMSolver_API void pcmsolver_delete(pcmsolver_context_t * context);

/*! \brief Whether the library is compatible with the header file
 *  Checks that the compiled library and header file version match.
 *  Host should abort when that is not the case.
 *  \warning This function should be called **before** instantiating
 *  any PCM context objects.
 */
PCMSolver_API pcmsolver_bool_t pcmsolver_is_compatible_library(void);

/*! \brief Prints citation and set up information
 *  \param[in, out] context the PCM context object
 */
PCMSolver_API void pcmsolver_print(pcmsolver_context_t * context);

/*! \brief Getter for the number of finite elements composing the molecular cavity
 *  \param[in, out] context the PCM context object
 *  \return the size of the cavity
 */
PCMSolver_API int pcmsolver_get_cavity_size(pcmsolver_context_t * context);

/*! \brief Getter for the number of irreducible finite elements composing the
 * molecular cavity
 *  \param[in, out] context the PCM context object
 *  \return the number of irreducible finite elements
 */
PCMSolver_API int pcmsolver_get_irreducible_cavity_size(
    pcmsolver_context_t * context);

/*! \brief Getter for the centers of the finite elements composing the molecular
 * cavity
 *  \param[in, out] context the PCM context object
 *  \param[out] centers array holding the coordinates of the finite elements centers
 */
PCMSolver_API void pcmsolver_get_centers(pcmsolver_context_t * context,
                                         double centers[]);

/*! \brief Getter for the center of the i-th finite element
 *  \param[in, out] context the PCM context object
 *  \param[in] its index of the finite element
 *  \param[out] center array holding the coordinates of the finite element center
 */
PCMSolver_API void pcmsolver_get_center(pcmsolver_context_t * context,
                                        int its,
                                        double center[]);

/*! \brief Getter for the areas/weights of the finite elements
 *  \param[in, out] context the PCM context object
 *  \param[out] areas array holding the weights/areas of the finite elements
 */
PCMSolver_API void pcmsolver_get_areas(pcmsolver_context_t * context,
                                       double areas[]);

/*! \brief Computes ASC given a MEP and the desired irreducible representation
 *  \param[in, out] context the PCM context object
 *  \param[in] mep_name label of the MEP surface function
 *  \param[in] asc_name label of the ASC surface function
 *  \param[in] irrep index of the desired irreducible representation
 *  The module uses the surface function concept to handle potentials
 *  and charges. Given labels for each, this function retrieves the MEP
 *  and computes the corresponding ASC.
 */
PCMSolver_API void pcmsolver_compute_asc(pcmsolver_context_t * context,
                                         const char * mep_name,
                                         const char * asc_name,
                                         int irrep);

/*! \brief Computes response ASC given a MEP and the desired irreducible
 * representation
 *  \param[in, out] context the PCM context object
 *  \param[in] mep_name label of the MEP surface function
 *  \param[in] asc_name label of the ASC surface function
 *  \param[in] irrep index of the desired irreducible representation
 *  If `Nonequilibrium = True` in the input, calculates a response
 *  ASC using the dynamic permittivity. Falls back to the solver with static
 * permittivity
 *  otherwise.
 */
PCMSolver_API void pcmsolver_compute_response_asc(pcmsolver_context_t * context,
                                                  const char * mep_name,
                                                  const char * asc_name,
                                                  int irrep);

/*! \brief Computes the polarization energy
 *  \param[in, out] context the PCM context object
 *  \param[in] mep_name label of the MEP surface function
 *  \param[in] asc_name label of the ASC surface function
 *  \return the polarization energy
 *  This function calculates the dot product of the given MEP and ASC vectors.
 */
PCMSolver_API double pcmsolver_compute_polarization_energy(
    pcmsolver_context_t * context,
    const char * mep_name,
    const char * asc_name);

/*! \brief Getter for the ASC dipole
 *  \param[in, out] context the PCM context object
 *  \param[in] asc_name label of the ASC surface function
 *  \param[out] dipole  the Cartesian components of the ASC dipole moment
 *  \return the ASC dipole, i.e. \sqrt{\sum_i \mu_i^2}
 */
PCMSolver_API double pcmsolver_get_asc_dipole(pcmsolver_context_t * context,
                                              const char * asc_name,
                                              double dipole[]);

/*! \brief Retrieves data wrapped in a given surface function
 *  \param[in, out] context the PCM context object
 *  \param[in] size the size of the surface function
 *  \param[in] values the values wrapped in the surface function
 *  \param[in] name label of the surface function
 */
PCMSolver_API void pcmsolver_get_surface_function(pcmsolver_context_t * context,
                                                  int size,
                                                  double values[],
                                                  const char * name);

/*! \brief Sets a surface function given data and label
 *  \param[in, out] context the PCM context object
 *  \param[in] size the size of the surface function
 *  \param[in] values the values to be wrapped in the surface function
 *  \param[in] name label of the surface function
 */
PCMSolver_API void pcmsolver_set_surface_function(pcmsolver_context_t * context,
                                                  int size,
                                                  double values[],
                                                  const char * name);

/*! \brief Prints surface function contents to host output
 *  \param[in, out] context the PCM context object
 *  \param[in] name label of the surface function
 */
PCMSolver_API void pcmsolver_print_surface_function(pcmsolver_context_t * context,
                                                    const char * name);

/*! \brief Dumps all currently saved surface functions to NumPy arrays in .npy files
 *  \param[in, out] context the PCM context object
 */
PCMSolver_API void pcmsolver_save_surface_functions(pcmsolver_context_t * context);

/*! \brief Dumps a surface function to NumPy array in .npy file
 *  \param[in, out] context the PCM context object
 *  \param[in] name label of the surface function
 *
 *  \note The name parameter is the name of the NumPy array file
 *  **without** .npy extension
 */
PCMSolver_API void pcmsolver_save_surface_function(pcmsolver_context_t * context,
                                                   const char * name);

/*! \brief Loads a surface function from a .npy file
 *  \param[in, out] context the PCM context object
 *  \param[in] name label of the surface function
 *
 *  \note The name parameter is the name of the NumPy array file
 *  **without** .npy extension
 */
PCMSolver_API void pcmsolver_load_surface_function(pcmsolver_context_t * context,
                                                   const char * name);

/*! \brief Writes timing results for the API
 *  \param[in, out] context the PCM context object
 */
PCMSolver_API void pcmsolver_write_timings(pcmsolver_context_t * context);

#ifdef __cplusplus
}
#endif
