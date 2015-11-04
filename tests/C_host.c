#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "pcmsolver.h"

#define NR_NUCLEI 6

#ifdef __GNUC__
#  define UNUSED(x) UNUSED_ ## x __attribute__((__unused__))
#else
#  define UNUSED(x) UNUSED_ ## x
#endif

#ifdef __GNUC__
#  define UNUSED_FUNCTION(x) __attribute__((__unused__)) UNUSED_ ## x
#else
#  define UNUSED_FUNCTION(x) UNUSED_ ## x
#endif

FILE * output;

PCMInput_t pcmsolver_input();

/*! \brief calculates nuclear molecular electrostatic potential (MEP) at cavity points
 *  \author Roberto Di Remigio
 *  \date 2015
 *  \param nr_nuclei number of atomic centers
 *  \param charges atomic charges
 *  \param coordinates coordinates of the atomic centers
 *  \param grid_size the number of cavity points
 *  \param grid the cavity points
 *  \return the nuclear MEP
 *  \warning Caller deallocats returned array
 */
double * nuclear_mep(int nr_nuclei, double charges[nr_nuclei],
    double coordinates[3*nr_nuclei], size_t grid_size, double grid[grid_size]);

/*! \brief Compares calculated and reference values within a threshold
 *  \author Roberto Di Remigio
 *  \date 2015
 *  \param calculated calculated value
 *  \param reference  reference value
 *  \param threshold  comparison threshold
 *  \return true if the calculated and reference values are compatible, false otherwise
 */
bool check_unsigned_error(double calculated, double reference, double threshold);

int main()
{

  output = fopen("C_host.log", "w+");
  if (!pcmsolver_is_compatible_library())
  {
    fprintf(stderr, "%s\n", "PCMSolver library not compatible");
    exit(EXIT_FAILURE);
  }

  fprintf(output, "%s\n", "Starting a PCMSolver calculation");
  double charges[NR_NUCLEI] = {6.0, 1.0, 1.0, 6.0, 1.0, 1.0};
  double coordinates[3 * NR_NUCLEI] = { 0.0,  0.000000,  1.257892,
                    0.0,  1.745462,  2.342716,
                    0.0, -1.745462,  2.342716,
                    0.0,  0.000000, -1.257892,
                    0.0,  1.745462, -2.342716,
                    0.0, -1.745462, -2.342716
                  };
  // This means the molecular point group has three generators:
  // the Oxy, Oxz and Oyz planes
  int symmetry_info[4] = {3, 4, 2, 1};
  PCMInput_t host_input = pcmsolver_input();

  pcmsolver_context_t * pcm_context = pcmsolver_new(PCMSOLVER_READER_HOST,
      NR_NUCLEI, charges, coordinates, symmetry_info, host_input);

  pcmsolver_print(pcm_context);

  size_t grid_size = pcmsolver_get_cavity_size(pcm_context);
  size_t irr_grid_size = pcmsolver_get_irreducible_cavity_size(pcm_context);
  double * grid = (double *) calloc(3*grid_size, sizeof(double));
  pcmsolver_get_centers(pcm_context, grid);

  double * mep = nuclear_mep(NR_NUCLEI, charges, coordinates, grid_size, grid);
  const char * mep_lbl = {"NucMEP"};
  pcmsolver_set_surface_function(pcm_context, grid_size, mep, mep_lbl);
  const char * asc_lbl = {"NucASC"};
  // This is the Ag irreducible representation (totally symmetric)
  int irrep = 0;
  pcmsolver_compute_asc(pcm_context, mep_lbl, asc_lbl, irrep);
  double * asc_Ag = (double *) calloc(grid_size, sizeof(double));
  pcmsolver_get_surface_function(pcm_context, grid_size, asc_Ag, asc_lbl);

  double energy = pcmsolver_compute_polarization_energy(pcm_context, mep_lbl, asc_lbl);

  fprintf(output, "Polarization energy: %20.12f\n", energy);

  double * asc_neq_B3g = (double *) calloc(grid_size, sizeof(double));
  const char * asc_neq_B3g_lbl = {"OITASC"};
  // This is the B3g irreducible representation
  irrep = 3;
  pcmsolver_compute_response_asc(pcm_context, mep_lbl, asc_neq_B3g_lbl, irrep);
  pcmsolver_get_surface_function(pcm_context, grid_size, asc_neq_B3g, asc_neq_B3g_lbl);

  // Equilibrium ASC in B3g symmetry.
  // This is an internal check: the relevant segment of the vector
  // should be the same as the one calculated using pcmsolver_compute_response_asc
  double * asc_B3g = (double *) calloc(grid_size, sizeof(double));
  const char * asc_B3g_lbl = {"ASCB3g"};
  pcmsolver_compute_asc(pcm_context, mep_lbl, asc_B3g_lbl, irrep);
  pcmsolver_get_surface_function(pcm_context, grid_size, asc_B3g, asc_B3g_lbl);

  // Check that everything calculated is OK
  // Cavity size
  const size_t ref_size = 576;
  if (grid_size != ref_size) {
    fprintf(stderr, "%s\n", "Error in the cavity size, please file an issue on: https://github.com/PCMSolver/pcmsolver");
    exit(EXIT_FAILURE);
  } else {
    fprintf(output, "%s\n", "Test on cavity size: PASSED");
  }
  // Irreducible cavity size
  const size_t ref_irr_size = 72;
  if (irr_grid_size != ref_irr_size) {
    fprintf(stderr, "%s\n", "Error in the irreducible cavity size, please file an issue on: https://github.com/PCMSolver/pcmsolver");
    exit(EXIT_FAILURE);
  } else {
    fprintf(output, "%s\n", "Test on irreducible cavity size: PASSED");
  }
  // Polarization energy
  const double ref_energy = -0.437960027982;
  if (!check_unsigned_error(energy, ref_energy, 1.0e-7)) {
    fprintf(stderr, "%s\n", "Error in the polarization energy, please file an issue on: https://github.com/PCMSolver/pcmsolver");
    exit(EXIT_FAILURE);
  } else {
    fprintf(output, "%s\n", "Test on polarization energy: PASSED");
  }
  // Surface functions
  //test_surface_functions(grid_size, mep, asc_Ag, asc_B3g, asc_neq_B3g)

  pcmsolver_write_timings(pcm_context);

  pcmsolver_delete(pcm_context);

  fclose(output);

  free(grid);
  free(mep);
  free(asc_Ag);
  free(asc_B3g);
  free(asc_neq_B3g);

  return EXIT_SUCCESS;
}

PCMInput_t pcmsolver_input()
{
  PCMInput_t host_input;

  // These parameters would be set by the host input reading
  // Length and area parameters are all assumed to be in Angstrom,
  // the module will convert to Bohr internally
  strcpy(host_input.cavity_type, "gepol");
  host_input.patch_level  = 2;
  host_input.coarsity     = 0.5;
  host_input.area         = 0.2;
  host_input.min_distance = 0.1;
  host_input.der_order    = 4;
  host_input.scaling      = true;
  strcpy(host_input.radii_set, "bondi");
  strcpy(host_input.restart_name, "cavity.npz");
  host_input.min_radius   = 100.0;

  strcpy(host_input.solver_type, "iefpcm");
  strcpy(host_input.solvent, "water");
  strcpy(host_input.equation_type, "secondkind");
  host_input.correction    = 0.0;
  host_input.probe_radius  = 1.0;

  strcpy(host_input.inside_type, "vacuum");
  host_input.outside_epsilon    = 1.0;
  strcpy(host_input.outside_type, "uniformdielectric");

  return host_input;
}

double * nuclear_mep(int nr_nuclei, double charges[nr_nuclei],
    double coordinates[3*nr_nuclei], size_t grid_size, double grid[3*grid_size])
{
  double * mep = (double *) calloc(grid_size, sizeof(double));
  for (int i = 0; i < nr_nuclei; i++) {
    for (size_t j = 0; j < grid_size; j++) {
      // Column-major ordering. Offsets: col_idx * nr_rows + row_idx
      double dist = pow((coordinates[i*3]     - grid[j*3]), 2)
                  + pow((coordinates[i*3 + 1] - grid[j*3 + 1]), 2)
                  + pow((coordinates[i*3 + 2] - grid[j*3 + 2]), 2);
      mep[j] += charges[i] / sqrt(dist);
    }
  }
  return mep;
}

void host_writer(const char * message, size_t UNUSED(message_length))
{
  fprintf(output, "%s\n", message);
}

bool check_unsigned_error(double calculated, double reference, double threshold)
{
  double err = calculated - reference;
  return (fabs(err) <= fabs(calculated) * threshold);
}

