#ifndef MOLECULE_H_
#define MOLECULE_H_
/************
 *  molecule.h  *
 ************/


/* number of charges */
extern const unsigned int k;

/* position of charges: vectors (x,y,z) */
extern const vector3 x[];


/* charges: double values */
extern const double alpha[];


int read_molecule(char * filename);

double potmol(vector3 a);

vector3 field(vector3 a);

void free_molecule();
#endif
