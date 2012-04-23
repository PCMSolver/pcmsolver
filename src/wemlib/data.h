#ifndef DATA
#define DATA
/************
 *  data.h  *
 ************/


/* number of charges */
extern const unsigned int k;

/* position of charges: vectors (x,y,z) */
extern const vector3 x[];


/* charges: double values */
extern const double alpha[];


double f(vector3 a);


vector3 df(vector3 a);
#endif
