#ifndef SPARSE2
#define SPARSE2
/**************
 *  sparse.h  *
 **************/


/*==================================================*
 *  Modul zur Rechnung mit duennbesetzten Matrizen  *
 *==================================================*/


/*====================*
 *  Typendeklaration  *
 *====================*/

typedef struct {
	unsigned int    m, n;
	unsigned int   *row_number, *max_row_number;
	unsigned int  **index;
	double        **value1;
	double        **value2;
} sparse2;



/*==========================================*
 *  Suchalgorithmus gemaess binary-search:  *
 *  Liefert den Index des Elements mit      *
 *  Index >= dem gesuchten Index j.         *
 *==========================================*/

unsigned int    search_sparse2(unsigned int *array, unsigned int rn, unsigned int j);
 /*
  * array = A->index[i] rn = Ende des Arrays, meist: A->row_number[i]
  */


/*=====================================*
 *  Initialisierung der sparse-Matrix  *
 *=====================================*/

void 
init_sparse2(sparse2 *A,
	     unsigned int m, unsigned int n, unsigned int n_max);
 /* fuer jede Zeile werden n_max Elemente belegt */


/*======================================*
 *  Initialisierung der Matrix-Pattern  *
 *======================================*/
void 
init_pattern(sparse2 *A,
	     unsigned int m, unsigned int n, unsigned int n_max);
 /* fuer jede Zeile werden n_max Elemente belegt */


/*===============================*
 *  Freigeben der sparse-Matrix  *
 *===============================*/

void            free_sparse2(sparse2 *A);


/*==================================*
 *  Fuege A(i,j) den Pattern hinzu  *
 *==================================*/

void            set_pattern(sparse2 *A, unsigned int i, unsigned int j);


/*=====================*
 *  A(i,j) := (z1,z2)  *
 *=====================*/

void            set_sparse2(sparse2 *A, unsigned int i, unsigned int j, double z1, double z2);


/*=====================*
 *  A(i,j) += (z1,z2)  *
 *=====================*/

void            add_sparse2(sparse2 *A, unsigned int i, unsigned int j, double z1, double z2);


/*=====================*
 *  (z1,z2) := A(i,j)  *
 *=====================*/

void            get_sparse2(sparse2 *A, unsigned int i, unsigned int j, double *z1, double *z2);


/*====================*
 *  nz = non_zero(A)  *
 *====================*/

unsigned int    non_zero(sparse2 *A);


/*=====================*
 *  Garbage-Collector  *
 *=====================*/

void            garbage_collect(sparse2 *A);


/*==================*
 *  Finish-Pattern  *
 *==================*/

void            finish_pattern(sparse2 *A);

void            print_sparse2(sparse2 *A);

#endif
