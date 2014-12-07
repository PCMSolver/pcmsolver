#ifndef SPARSE
#define SPARSE
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
    unsigned int m, n;
    unsigned int *row_number, *max_row_number;
    unsigned int **index;
    double **value;
} sparse;



/*==========================================*
 *  Suchalgorithmus gemaess binary-search:  *
 *  Liefert den Index des Elements mit      *
 *  Index >= dem gesuchten Index j.         *
 *==========================================*/

unsigned int search_sparse(unsigned int *array, unsigned int rn, unsigned int j);
 /*
  * array = A->index[i] rn = Ende des Arrays, meist: A->row_number[i]
  */


/*=====================================*
 *  Initialisierung der sparse-Matrix  *
 *=====================================*/

void init_sparse(sparse *A, unsigned int m, unsigned int n, unsigned int n_max);
 /* fuer jede Zeile werden n_max Elemente belegt */


/*===============================*
 *  Freigeben der sparse-Matrix  *
 *===============================*/

void free_sparse(sparse *A);


/*===============*
 *  A(i,j) := z  *
 *===============*/

void set_sparse(sparse *A, unsigned int i, unsigned int j, double z);


/*===============*
 *  A(i,j) += z  *
 *===============*/

void add_sparse(sparse *A, unsigned int i, unsigned int j, double z);


/*===============*
 *  z := A(i,j)  *
 *===============*/

double get_sparse(sparse *A, unsigned int i, unsigned int j);
#endif
