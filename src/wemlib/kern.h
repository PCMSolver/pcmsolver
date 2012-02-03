#ifndef KERN
#define KERN
/**********
 * kern.h *
 **********/


extern const double epsilon;	/* dielectric constant of the solvent */
extern const double kappa;	/* constant related to ion screening  */

/* inverse dielectric tensor of the solvent */
extern const double epsilon11;
extern const double epsilon12;
extern const double epsilon13;
extern const double epsilon21;
extern const double epsilon22;
extern const double epsilon23;
extern const double epsilon31;
extern const double epsilon32;
extern const double epsilon33;


/*===========================*
 *  Einfachschichtpotential  *
 *===========================*/

double          SingleLayerInt(vector3 x, vector3 y);


double          SingleLayerExt(vector3 x, vector3 y);


double          SingleLayerAni(vector3 x, vector3 y);


/*==========================*
 *  Doppelschichtpotential  *
 *==========================*/

double          DoubleLayerInt(vector3 x, vector3 y, vector3 n_y);


double          DoubleLayerExt(vector3 x, vector3 y, vector3 n_y);


double          DoubleLayerAni(vector3 x, vector3 y, vector3 n_y);
#endif
