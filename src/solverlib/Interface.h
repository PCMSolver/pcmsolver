#ifndef INTERFACE_H_
#define INTERFACE_H_
/*

  Interface functions prototypes.

*/


/*

	Functions visible to host program 

*/

extern "C" void init_pcm_();

extern "C" void comp_chg_pcm_(char* potString, char* chgString);

extern "C" void comp_pol_ene_pcm_(double * energy);

extern "C" void get_epsilon_static_(double * epsilon);

extern "C" void collect_nctot_(int * nuclei);

extern "C" void collect_atoms_(double * charges, double * centers, int * units);

extern "C" void get_cavity_size_(int * nts);

extern "C" void get_tess_centers_(double * centers);

extern "C" void get_tess_cent_coord_(int * its, double * center);

extern "C" void print_gepol_cavity_();

extern "C" void set_surface_function_(int * nts, double * values, char * name);

extern "C" void get_surface_function_(int * nts, double * values, char * name);

extern "C" void add_surface_function_(char * result, double * coeff, char * part);

extern "C" void print_surface_function_(char * name);

extern "C" bool surf_func_exists_(char * name);

extern "C" void clear_surf_func_(char * name);

extern "C" void append_surf_func_(char * name);

/*

	Functions not visible to host program	

 */

void init_gepol_cavity_();

void init_wavelet_cavity_();

void init_iefsolver_();

void init_pwcsolver_();

void init_pwlsolver_();

void build_isotropic_matrix_();

void build_anisotropic_matrix_();

void init_atoms_(VectorXd & charges, 
                 Matrix<double, 3, Dynamic> & sphereCenter);

void init_spheres_implicit_(VectorXd & charges, 
                 Matrix<double, 3, Dynamic> & centers);

void init_spheres_atoms_(VectorXd & charges, 
                                    Matrix<double, 3, Dynamic> & centers);


#endif

