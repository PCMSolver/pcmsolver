/**
 * PCMSolver, an API for the Polarizable Continuum Model
 * Copyright (C) 2016 Roberto Di Remigio, Luca Frediani and collaborators.
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

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "PCMInput.h"
#include "pcmsolver.h"

#include "C_host-functions.h"

#define NR_NUCLEI 6

FILE * output;

void host_writer(const char * message) { fprintf(output, "%s\n", message); }

int main() {

  output = fopen("C_host.out", "w+");
  if (!pcmsolver_is_compatible_library()) {
    fprintf(stderr, "%s\n", "PCMSolver library not compatible");
    exit(EXIT_FAILURE);
  }

  fprintf(output, "%s\n", "Starting a PCMSolver calculation");
  // Use C2H4 in D2h symmetry
  double charges[NR_NUCLEI] = {6.0, 1.0, 1.0, 6.0, 1.0, 1.0};
  double coordinates[3 * NR_NUCLEI] = {0.0,
                                       0.000000,
                                       1.257892,
                                       0.0,
                                       1.745462,
                                       2.342716,
                                       0.0,
                                       -1.745462,
                                       2.342716,
                                       0.0,
                                       0.000000,
                                       -1.257892,
                                       0.0,
                                       1.745462,
                                       -2.342716,
                                       0.0,
                                       -1.745462,
                                       -2.342716};
  // This means the molecular point group has three generators:
  // the Oxy, Oxz and Oyz planes
  int symmetry_info[4] = {3, 4, 2, 1};
  struct PCMInput host_input = pcmsolver_input();

  pcmsolver_context_t * pcm_context = pcmsolver_new(PCMSOLVER_READER_HOST,
                                                    NR_NUCLEI,
                                                    charges,
                                                    coordinates,
                                                    symmetry_info,
                                                    &host_input,
                                                    host_writer);

  pcmsolver_print(pcm_context);

  int grid_size = pcmsolver_get_cavity_size(pcm_context);
  int irr_grid_size = pcmsolver_get_irreducible_cavity_size(pcm_context);
  double * grid = (double *)calloc(3 * grid_size, sizeof(double));
  pcmsolver_get_centers(pcm_context, grid);
  double * areas = (double *)calloc(grid_size, sizeof(double));
  pcmsolver_get_areas(pcm_context, areas);

  double * mep = nuclear_mep(NR_NUCLEI, charges, coordinates, grid_size, grid);
  const char * mep_lbl = {"NucMEP"};
  pcmsolver_set_surface_function(pcm_context, grid_size, mep, mep_lbl);
  const char * asc_lbl = {"NucASC"};
  // This is the Ag irreducible representation (totally symmetric)
  int irrep = 0;
  pcmsolver_compute_asc(pcm_context, mep_lbl, asc_lbl, irrep);
  double * asc_Ag = (double *)calloc(grid_size, sizeof(double));
  pcmsolver_get_surface_function(pcm_context, grid_size, asc_Ag, asc_lbl);

  double energy =
      pcmsolver_compute_polarization_energy(pcm_context, mep_lbl, asc_lbl);

  fprintf(output, "Polarization energy: %20.12f\n", energy);

  double * asc_neq_B3g = (double *)calloc(grid_size, sizeof(double));
  const char * asc_neq_B3g_lbl = {"OITASC"};
  // This is the B3g irreducible representation
  irrep = 3;
  pcmsolver_compute_response_asc(pcm_context, mep_lbl, asc_neq_B3g_lbl, irrep);
  pcmsolver_get_surface_function(
      pcm_context, grid_size, asc_neq_B3g, asc_neq_B3g_lbl);

  // Equilibrium ASC in B3g symmetry.
  // This is an internal check: the relevant segment of the vector
  // should be the same as the one calculated using pcmsolver_compute_response_asc
  double * asc_B3g = (double *)calloc(grid_size, sizeof(double));
  const char * asc_B3g_lbl = {"ASCB3g"};
  pcmsolver_compute_asc(pcm_context, mep_lbl, asc_B3g_lbl, irrep);
  pcmsolver_get_surface_function(pcm_context, grid_size, asc_B3g, asc_B3g_lbl);

  // Check that everything calculated is OK
  // Cavity size
  const int ref_size = 576;
  if (grid_size != ref_size) {
    fprintf(stderr,
            "%s\n",
            "Error in the cavity size, please file an issue on: "
            "https://github.com/PCMSolver/pcmsolver");
    exit(EXIT_FAILURE);
  } else {
    fprintf(output, "%s\n", "Test on cavity size: PASSED");
  }
  // Irreducible cavity size
  const int ref_irr_size = 72;
  if (irr_grid_size != ref_irr_size) {
    fprintf(stderr,
            "%s\n",
            "Error in the irreducible cavity size, please file an "
            "issue on: https://github.com/PCMSolver/pcmsolver");
    exit(EXIT_FAILURE);
  } else {
    fprintf(output, "%s\n", "Test on irreducible cavity size: PASSED");
  }
  // Polarization energy
  const double ref_energy = -0.437960027982;
  if (!check_unsigned_error(energy, ref_energy, 1.0e-7)) {
    fprintf(stderr,
            "%s\n",
            "Error in the polarization energy, please file an issue "
            "on: https://github.com/PCMSolver/pcmsolver");
    exit(EXIT_FAILURE);
  } else {
    fprintf(output, "%s\n", "Test on polarization energy: PASSED");
  }
  // Surface functions
  test_surface_functions(
      output, grid_size, mep, asc_Ag, asc_B3g, asc_neq_B3g, areas);

  pcmsolver_save_surface_functions(pcm_context);
  pcmsolver_save_surface_function(pcm_context, asc_lbl);
  pcmsolver_load_surface_function(pcm_context, mep_lbl);

  pcmsolver_write_timings(pcm_context);

  pcmsolver_delete(pcm_context);

  free(grid);
  free(mep);
  free(asc_Ag);
  free(asc_B3g);
  free(asc_neq_B3g);
  free(areas);

  fclose(output);

  return EXIT_SUCCESS;
}
