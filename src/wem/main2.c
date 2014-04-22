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

#include "dalton_wem.h"
#include "vector3.h"
#include <stdio.h>
#include <stdlib.h>


const unsigned int k_m = 12;
const vector3 x_m[12] = {
    {5.274, 1.999, -8.568},
    {6.627, 2.018, -8.209},
    {7.366, 0.829, -8.202},
    {6.752, -0.379, -8.554},
    {5.399, -0.398, -8.912},
    {4.660, 0.791, -8.919},
    {4.704, 2.916, -8.573},
    {7.101, 2.950, -7.938},
    {8.410, 0.844, -7.926},
    {7.322, -1.296, -8.548},
    {4.925, -1.330, -9.183},
    {3.616, 0.776, -9.196}
};
const double alpha_m[12] = { 6, 6, 6, 6, 6, 6, 1, 1, 1, 1, 1, 1 };

/*
const unsigned int k_m=2;
const vector3 x_m[2] = {
  {0.0, 0.0, 0.0},
  {0.0, 0.0, 0.1}
};
const double alpha_m[2]={1.0,1.0};
*/

double f_m(vector3 a)
{
    double c = 0;
    unsigned int i;
    for (i = 0; i < k_m; i++)
        c += alpha_m[i] / vector3_norm(vector3_make(a.x - x_m[i].x, a.y - x_m[i].y, a.z - x_m[i].z));
    return (c);
}

double total_charge_inside_m()
{
    unsigned int i;
    double tot = 0.0;
    for (i = 0; i < k_m; i++)
        tot += alpha_m[i];
    return (tot);
}

int main()
{

    int num_points, i, num_charges;
    double *xpotential_coordinates;
    double *ypotential_coordinates;
    double *zpotential_coordinates;
    double *potential_values;
    double *charges;
    double *xcharge_positions;
    double *ycharge_positions;
    double *zcharge_positions;
    double x, y, z, energy, energy1;

    dalton_wem_initialize_(&num_points);

    xpotential_coordinates = (double *) malloc(num_points * sizeof(double));
    ypotential_coordinates = (double *) malloc(num_points * sizeof(double));
    zpotential_coordinates = (double *) malloc(num_points * sizeof(double));
    potential_values = (double *) malloc(num_points * sizeof(double));

    dalton_wem_get_potential_coordinates_(xpotential_coordinates, ypotential_coordinates, zpotential_coordinates);

    //Calculate potential
    printf("Potential:\n");
    for (i = 0; i < num_points; i++) {
        x = xpotential_coordinates[i];
        y = ypotential_coordinates[i];
        z = zpotential_coordinates[i];
        //if (i < 5)
        printf("%f %f %f - %lf \n", x, y, z, potential_values[i]);
        potential_values[i] = f_m(vector3_make(x, y, z));
    }


    /*
       // read potential from file
       FILE *in=fopen("pot.dat","r");
       if(in==NULL) {
       printf("terrorerror\n");
       return(-1);
       }
       for(i=0;i<num_points;i++){
       fscanf(in,"%lf %lf %lf - %lf\n",&x,&y,&z,&(potential_values[i]));
       printf("%lf\n",potential_values[i]);
       }
       fclose(in);
     */

    xcharge_positions = (double *) malloc(num_points * sizeof(double));
    ycharge_positions = (double *) malloc(num_points * sizeof(double));
    zcharge_positions = (double *) malloc(num_points * sizeof(double));
    charges = (double *) malloc(num_points * sizeof(double));

    dalton_wem_get_charges_(&num_points, &num_charges, potential_values, charges, xcharge_positions, ycharge_positions, zcharge_positions);

    //Calculate energy
    printf("Charges:\n");
    energy1 = 0.0;
    for (i = 0; i < num_charges; i++) {
        x = xcharge_positions[i];
        y = ycharge_positions[i];
        z = zcharge_positions[i];
        //if (i < 5)
        printf("%f %f %f - %lf \n", x, y, z, charges[i]);
        //energy1 += f_m(vector3_make(x, y, z)) * charges[i];
        energy1 += potential_values[i] * charges[i];
    }
    energy1 = -0.5 * energy1;
    //factor from Helmut

    // Using full quadrature in energy calculations
    dalton_wem_energy_(potential_values, &energy);

    printf("Energy, Helmut: %lf    Ville: %lf\n", energy, energy1);
    printf("Energy divided by au2angs   : %lf\n", energy * 0.5291772);
    dalton_wem_finalize_();

    free(xpotential_coordinates);
    free(ypotential_coordinates);
    free(zpotential_coordinates);
    free(potential_values);
    free(charges);
    free(xcharge_positions);
    free(ycharge_positions);
    free(zcharge_positions);

    return 0;
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

