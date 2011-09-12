/*************
 *  Gamma.h  *
 *************/


/*===========================================*
 *  Dummy-Header fuer die Parametrisierung,  *
 *  im makefile ist dann das entsprechende   *
 *  c-File zu setzen.			     *
 *===========================================*/


typedef struct 
{
vector3		(*f)(vector2 a);		/* Parametrisierung     */
vector3		(*df_dx)(vector2 a);		/* Ableitung nach x     */
vector3		(*df_dy)(vector2 a);		/* Ableitung nach y     */
vector3		(*n_f)(vector2 a);		/* zugehoerende Normale */
}  parametrix;


unsigned int init_p(void);
/* Initialisierung: Liefert als Funktionsergebnis die 
   Anzahl p der Parametergebiete */

  
void init_Chi(parametrix **Chi);
/* allokiert den noetigen Speicherplatz fuer die 
   Parametrisierung und definiert Chi[0],...,Chi[p] */


void free_Chi(parametrix **Chi);
/* gibt den Speicherplatz fuer die Parametrisierung frei */
