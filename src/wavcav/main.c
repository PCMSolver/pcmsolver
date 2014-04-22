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
#endif

/* warning-disabler-end */

#include "test.h"

int main(int argc, char* argv[])
{
    char* infile = "cavity.inp";
    double probeRadius = 0.4;
    double coarsity = 0.4;
    int patchLevel = 2;
    int info;
    printf("nr of args %d\n", argc); 
    switch (argc) {
    case 5:
        coarsity = atof(argv[4]);
    case 4:
        probeRadius = atof(argv[3]);
    case 3:
        patchLevel = atoi(argv[2]);
    case 2:
        infile = argv[1];
    case 1:
        break;
    default :
        printf("Usage:\n"); 
        printf("cav.x input_file [patchLevel[probeRadius[coarsity]]]\n"); 
        exit(-1);
    }
    printf ("Calling cavity generator with the following input\n");
    printf ("Input file: %s Patch Level: %2d   Radius: %7.4f   Coarsity: %7.4f\n",
            infile, patchLevel, probeRadius, coarsity);
    info = waveletCavityDrv_(probeRadius, coarsity, patchLevel, infile);
    return info;
};

/* warning-disabler-start */

#if (defined(__GNUC__) || defined(__GNUG__)) && !(defined(__clang__) || defined(__INTEL_COMPILER))
#pragma GCC diagnostic pop
#elif defined(__ICC) || defined(__INTEL_COMPILER)
#pragma warning pop
#elif defined(__clang__)
#pragma clang diagnostic pop
#endif

/* warning-disabler-end */

