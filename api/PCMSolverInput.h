/*
 * PCMSolver, an API for the Polarizable Continuum Model
 * Copyright (C) 2018 Roberto Di Remigio, Luca Frediani and contributors.
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

/*! \typedef OutputWriter
 *  Flushes module output to host program
 *  \param[in,out] message contents of the module output
 */
typedef void (*OutputWriter)(const char * message);

/*! \struct PCMSolverInput
 *  \brief Data structure for host-API input communication.
 *
 *  The input is expected to be in **atomic units**
 */
struct PCMSolverInput {
  /// Number of atomic centers
  int nr_nuclei;
  /// Array of atomic charges
  double * charges;
  /*! Array of atomic coordinates
   *  These are expected in the format 3*nr_nuclei
   */
  double * coordinates;
  /*! Array specifying the generators of the molecular point group.
   *  The molecular point group information is passed as an array
   *  of 4 integers: number of generators, first, second and third generator
   *  respectively. Generators map to integers as in table :ref:`symmetry-ops`
   */
  int * symmetry_info;
  /// Function pointer for module output flushing
  OutputWriter writer;
  /// Type of cavity requested.
  const char * cavity_type;
  /// Name of the .npz file for GePol cavity restart.
  const char * restart_name;
  /// Average tesserae area.
  double area;
  /// Whether to scale or not the atomic radii.
  pcmsolver_bool_t scaling;
  /// The built-in radii set to be used.
  const char * radii_set;
  /// Minimal radius for the added spheres.
  double min_radius;
  const char * mode;
  int * atoms;
  double * radii;
  double * spheres;

  /// Type of solver requested.
  const char * solver_type;
  pcmsolver_bool_t nonequilibrium;
  /// Name of the solvent.
  const char * solvent;
  pcmsolver_bool_t matrix_symmetry;
  /// Correction in the CPCM apparent surface charge scaling factor.
  double correction;
  const char * diagonal_integrator;
  double diagonal_scaling;
  /// Radius of the spherical probe mimicking the solvent.
  double probe_radius;

  /// Type of Green's function requested inside the cavity.
  const char * inside_type;
  const char * inside_derivative;
  /// Value of the static permittivity outside the cavity.
  double outside_epsilon;
  /// Type of Green's function requested outside the cavity.
  const char * outside_type;
};
