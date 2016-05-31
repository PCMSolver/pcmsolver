#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "pcmsolver.h"
#include "PCMInput.h"

#include "C_host-functions.h"

#define NR_NUCLEI 6

FILE * output;

void host_writer(const char * message, int UNUSED(message_length))
{
  fprintf(output, "%s\n", message);
}

int main()
{

  output = fopen("fail-C_host.log", "w+");
  if (!pcmsolver_is_compatible_library())
  {
    fprintf(stderr, "%s\n", "PCMSolver library not compatible");
    exit(EXIT_FAILURE);
  }

  fprintf(output, "%s\n", "Starting a PCMSolver calculation");
  // Use C2H4 in D2h symmetry
  double charges[NR_NUCLEI] = {42.0, 1.0, 1.0, 6.0, 1.0, 1.0};
  double coordinates[3 * NR_NUCLEI] = { 0.0,  0.000000,  1.257892,
                    0.0,  1.745462,  2.342716,
                    0.0, -1.745462,  2.342716,
                    0.0,  0.000000, -1.257892,
                    0.0,  1.745462, -2.342716,
                    0.0, -1.745462, -2.342716
                  };
  // This means the molecular point group has three generators:
  // the Oxy, Oxz and Oyz planes
  int symmetry_info[4] = {0, 0, 0, 0};
  struct PCMInput host_input = pcmsolver_input();

  pcmsolver_context_t * pcm_context = pcmsolver_new(PCMSOLVER_READER_HOST,
      NR_NUCLEI, charges, coordinates, symmetry_info, &host_input);

  pcmsolver_delete(pcm_context);

  fclose(output);

  return EXIT_SUCCESS;
}
