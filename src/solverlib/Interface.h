#ifndef INTERFACE_H_
#define INTERFACE_H_
/*

  Interface functions prototypes.

*/


/*

  Cavity related functions

*/

extern "C" void init_gepol_cavity_();

extern "C" void collect_nctot_(int * nuclei);

extern "C" void collect_atoms_(double * charges, double * centers, int * flag);

extern "C" void init_atoms_(VectorXd & charges, 
                            Matrix<double, 3, Dynamic> & sphereCenter);

extern "C" void init_wavelet_cavity_();

extern "C" void get_cavity_size_(int * nts);

extern "C" void get_total_surface_charge_(double * charge);

extern "C" void get_nuclear_surface_charge_(double * charge);

extern "C" void get_electronic_surface_charge_(double * charge);

extern "C" void get_tess_centers_(double * centers);

extern "C" void comp_pot_chg_pcm_(double *density, double *work, int *lwork);

extern "C" void comp_pol_ene_pcm_(double * energy, int & printlevel);

extern "C" void print_gepol_cavity_();


/*

  Solver related functions

 */

extern "C" void get_epsilon_static_(double * epsilon);

extern "C" void init_pcm_();

extern "C" void init_iefsolver_();

extern "C" void init_wemsolver_();

extern "C" void build_isotropic_matrix_();

extern "C" void build_anisotropic_matrix_();

//copying mechanism of the following routine needs to be revised
extern "C" void comp_charge_(double *potential_, double *charge_, int & type);

//      Subroutine PotExpVal(Density, Centers, Nts, Potential, Work, 
//     $                     LWork)
//extern "C" void ele_pot_pcm_(double *density, double* centers, int *nts, 
//			     double *potential, double *work, int *lwork);

extern "C" void ele_pot_pcm_(int * nr_points, double * tess_cent, double * el_pot,
                            double * density, double * work, int * lwork);

//extern "C" void nuc_pot_pcm_(double* centers, int *nts, double *potential);
extern "C" void nuc_pot_pcm_(int * nts, double * tess_cent, double * nuc_pot);

//      Subroutine Fock_PCMModule(Fock, Centers, Nts, Charges, Work, 
//     $                     LWork)
extern "C" void fock_pcm_(double *fock, double* centers, int *nts, 
			  double *charges, double *work, int *lwork);

#endif

