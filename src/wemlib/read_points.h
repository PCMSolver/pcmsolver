/*******************
 *  read_points.h  *
 *******************/

#ifdef __cplusplus
extern "C" {
#endif


/*==============================*
 *  Liest die Punkteliste ein.  *
 *==============================*/
 
void read_points(vector3 ****P, unsigned int *p, unsigned int *m);

int read_points1(vector3 ****P, unsigned int *p, unsigned int *m, char *fname);

void free_points(vector3 ****P, unsigned int p, unsigned int m);


#ifdef __cplusplus
}
#endif
