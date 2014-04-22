/* warning-disabler-start */

#if (defined(__GNUC__) || defined(__GNUG__)) && !(defined(__clang__) || defined(__INTEL_COMPILER))
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wall"
#pragma GCC diagnostic ignored "-Weffc++"
#pragma GCC diagnostic ignored "-Wextra"
#elif defined(__ICC) || defined(__INTEL_COMPILER)
#pragma warning push
#pragma warning disable "-Wall"
#elif defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wall"
#pragma clang diagnostic ignored "-Weffc++"
#pragma clang diagnostic ignored "-Wextra"
#pragma clang diagnostic ignored "-Wdocumentation"
#pragma clang diagnostic ignored "-Wdeprecated-register"
#pragma clang diagnostic ignored "-Wincompatible-pointer-types"
#pragma clang diagnostic ignored "-Wempty-body"
#endif

/* warning-disabler-end */


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "vector3.h"

int main(int argc, char **argv)
{

    unsigned int n = 0;         /* n*n Elemente pro Patch auf dem Level m */
    unsigned int k, l;          /* Laufindizes                            */
    unsigned int j1, j2, j3;    /* Punktindizes                           */
    unsigned int m, p;
    float x, y, z;              /* Punktkoordinaten                       */
    FILE *infile;               /* Ausgabe-File                           */
    FILE *outfile;              /* Ausgabe-File                           */

    printf("Infile:%s\n", argv[1]);
    printf("Outfile:%s\n", argv[2]);

    /* file with geometric data */
    infile = fopen(argv[1], "r");
    outfile = fopen(argv[2], "w");

    fscanf(infile, "%d\n%d\n", &m, &p);
    fprintf(outfile, "%d\n%d\n", m, p);
    n = 1 << m;

    for (k = 0; k < p * (n + 1) * (n + 1); k++) {
        fscanf(infile, "%d %d %d %g %g %g\n", &j1, &j3, &j2, &x, &y, &z);
        fprintf(outfile, "%d %d %d %.16f %.16f %.16f\n", j1, j2, j3, x, y, z);
    }

    fclose(infile);
    fclose(outfile);

}
/* warning-disabler-start */

#if (defined(__GNUC__) || defined(__GNUG__)) && !(defined(__clang__) || defined(__INTEL_COMPILER))
#pragma GCC diagnostic pop
#elif defined(__ICC) || defined(__INTEL_COMPILER)
#pragma warning pop
#elif defined(__clang__)
#pragma clang diagnostic pop
#endif

/* warning-disabler-end */

