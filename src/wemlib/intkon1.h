/***************
 *  IntKon1.h  *
 ***************/


/*===========================================================*
 *  Fernfeld-Quadratur-Routine:			             *
 *  Die in der Funktion init_randwerte vorab berechneten     *
 *  Auswertepunkte und Gewichte der Gauss-Quadratur werden   *	
 *  in IntKon1 zum entsprechenden Integral zusammengefuegt.  *
 *===========================================================*/
 

typedef struct
{
vector3		*Chi, *n_Chi;
double		*det_dChi;
unsigned int	nop;
}  
randwerte;


void init_randwerte(randwerte **RW, unsigned int g_max);
/* Initialisiert die Randwerte */


void reset_randwerte(randwerte *RW, unsigned int g_max);
/* Resetted die Randwerte */


void free_randwerte(randwerte **RW, unsigned int g_max);
/* Gibt den Speicherplatz fuer die Randwerte frei */


void IntKon1(double *c, element *element1, element *element2, randwerte *RW, 
	cubature *Q1, cubature *Q2, vector3 ****P, unsigned int M, double (*SL)(), double (*DL)());
/* No-Problem-Quadrature-Routine */
