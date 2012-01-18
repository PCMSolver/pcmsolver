#ifndef DALTON_WEM
#define DALTON_WEM

/* 
 * Initializes the great WEM machine!
 *
 *   int (out)      - number of potential points
 *   double(out)[3] - a,b,dp wavelet parameters
 *   double(out)[2] - sparsities of the system matrix after and before (in %)
 *
 * return 0 - if all ok
 *
 */
int dalton_wem_initialize_(int *num_points, double *da, double *db, double *ddp,
			   double *apriori, double *aposteriori);


/*
 * Gets points for potential evaluation
 *
 *   data[number of potential points] (out) :
 *     [parameter domain][patch_level_y][patch_level_x][quadrature point]
 *   double * (out) - [x1,x2, ...]
 *   double * (out) - [y1,y2, ...]
 *   double * (out) - [z1,z2, ...]
 *   double * (out) - Area of the cavity
 *     memory should be provided by dalton
 *
 */
void dalton_wem_get_potential_coordinates_(double *xp, double *yp, double *zp, double *carea);

/*
 * as - areas of patches
 * xn - x-component of the normal
 * ...
 */
void dalton_wem_get_areas_normals_(double *as, double *xn, double *yn, double *zn);

/*
 * Calculates surface polarization charges
 *
 *   int * (in)       - maximum number of charge points
 *   int * (out)      - number of charges
 *   double * (in)  - potential values at pre-specified points
 *   double * (out) - resulting charges at the centers of subpatches
 *                    (allocate in DALTON, should be larger than #(potential))
 *
 */
void dalton_wem_get_charges_(int *max_charges, int *num_charges, double *potential, double *charges);

/*
 * Calculates surface polarization charges
 *
 *   int * (in)       - maximum number of charge points
 *   int * (out)      - number of charges
 *   double * (in)  - field values at pre-specified points, x-component
 *   double * (in)  - field values at pre-specified points, y-component
 *   double * (in)  - field values at pre-specified points, z-component
 *   double * (out) - resulting charges at the centers of subpatches
 *                    (allocate in DALTON, should be larger than #(potential))
 *
 */
void dalton_wem_get_charges2_(int *max_charges, int *num_charges, double *fieldx, 
			      double *fieldy, double *fieldz, double *charges);

/*
 * Calculates the total energy
 *
 * double *potential;   Potential points!
 * double *energy;      Results
 *
 */
void dalton_wem_energy_(double *potential, double *energy);

void dalton_wem_finalize_();


#endif
