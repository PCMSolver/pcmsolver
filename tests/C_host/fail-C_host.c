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

  output = fopen("fail-C_host.out", "w+");
  if (!pcmsolver_is_compatible_library()) {
    fprintf(stderr, "%s\n", "PCMSolver library not compatible");
    exit(EXIT_FAILURE);
  }

  fprintf(output, "%s\n", "Starting a PCMSolver calculation");
  // Use C2H4 in D2h symmetry
  double charges[NR_NUCLEI] = {42.0, 1.0, 1.0, 6.0, 1.0, 1.0};
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
  int symmetry_info[4] = {0, 0, 0, 0};
  struct PCMInput host_input = pcmsolver_input();

  pcmsolver_context_t * pcm_context = pcmsolver_new(PCMSOLVER_READER_HOST,
                                                    NR_NUCLEI,
                                                    charges,
                                                    coordinates,
                                                    symmetry_info,
                                                    &host_input,
                                                    host_writer);

  pcmsolver_delete(pcm_context);

  fclose(output);

  return EXIT_SUCCESS;
}
