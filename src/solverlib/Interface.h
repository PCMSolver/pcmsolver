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

extern "C" void get_tess_centers_(double * centers);

extern "C" void comp_pot_chg_pcm_(double *density, double *work, int *lwork);

extern "C" void comp_pol_ene_pcm_(double * energy);

extern "C" void print_gepol_cavity_();

extern "C" void get_tess_cent_coord_(int * its, double * center);

/*

  Solver related functions

 */

extern "C" void get_epsilon_static_(double * epsilon);

extern "C" void init_pcm_();

extern "C" void init_iefsolver_();

extern "C" void init_pwcsolver_();

extern "C" void init_pwlsolver_();

extern "C" void build_isotropic_matrix_();

extern "C" void build_anisotropic_matrix_();

extern "C" void comp_charge_(double *potential_, double *charge_);

extern "C" void ele_pot_pcm_(int * nts, double * centers, 
                             double * potential, double * density, 
                             double * work, int * lwork);

extern "C" void nuc_pot_pcm_(int * nts, double * tess_cent, 
                             double * nuc_pot);

extern "C" void fock_pcm_(double * fock, double * centers, int * nts, 
			  double * charges, double * work, int * lwork);

extern "C" void init_spheres_implicit_(VectorXd & charges, 
                            Matrix<double, 3, Dynamic> & centers);

extern "C" void init_spheres_atoms_(VectorXd & charges, 
                                    Matrix<double, 3, Dynamic> & centers);

extern "C" void comp_pot_chg_pcm_(double *density, double *work, int *lwork);

extern "C" void comp_chg_pcm_(char* potString, char* chgString);

extern "C" void set_surface_function_(int * nts, double * values, char * name);

extern "C" void get_surface_function_(int * nts, double * values, char * name);

extern "C" void add_surface_function_(char * result, double * coeff, char * part);

extern "C" void print_surface_function_(char * name);

#endif

