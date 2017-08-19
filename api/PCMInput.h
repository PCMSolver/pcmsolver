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

/*! @struct PCMInput
 *  @brief Data structure for host-API input communication.
 */
typedef struct PCMInput {
  /// Type of cavity requested.
  char cavity_type[8];
  /// Wavelet cavity mesh patch level.
  int patch_level;
  /// Wavelet cavity mesh coarsity.
  double coarsity;
  /// Average tesserae area.
  double area;
  /// The built-in radii set to be used.
  char radii_set[8];
  /// Minimal distance between sampling points.
  double min_distance;
  /// Derivative order for the switching function.
  int der_order;
  /// Whether to scale or not the atomic radii.
  pcmsolver_bool_t scaling;
  /// Name of the .npz file for GePol cavity restart.
  char restart_name[20];
  /// Minimal radius for the added spheres.
  double min_radius;
  /// Type of solver requested.
  char solver_type[7];
  /// Correction in the CPCM apparent surface charge scaling factor.
  double correction;
  /// Name of the solvent.
  char solvent[16];
  /// Radius of the spherical probe mimicking the solvent.
  double probe_radius;
  /// Type of the integral equation to be used.
  char equation_type[11];
  /// Type of Green's function requested inside the cavity.
  char inside_type[7];
  /// Value of the static permittivity outside the cavity.
  double outside_epsilon;
  /// Type of Green's function requested outside the cavity.
  char outside_type[22];
} PCMInput;
