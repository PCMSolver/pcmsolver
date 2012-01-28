#ifndef PHI
#define PHI
/***********
 *  Phi.h  *
 ***********/

 
/*============================================*
 *  Definiert die vier stueckweise linearen   *
 *  Basisfunktionen auf dem Einheitsquadrat.  *
 *  Es ist Phi0([0,0]) = 1,		      *
 *         Phi1([1,0]) = 1, 		      *
 * 	   Phi2([1,1]) = 1,		      *
 *         Phi3([0,1]) = 1.		      *
 *============================================*/


/* Ansatzfunktion 0 */
double Phi0(vector2 a);


/* Ansatzfunktion 1 */
double Phi1(vector2 a);


/* Ansatzfunktion 2 */
double Phi2(vector2 a);


/* Ansatzfunktion 3 */
double Phi3(vector2 a);


/* updated c_{i,j} um weight*phi_i(xi)*phi_j(eta) */ 
void Phi_times_Phi(double *c, double weight, vector2 xi, vector2 eta);


/* updated c_{i,j} um weight * < curl[phi_i(xi)],curl[phi_j(eta)] > */ 
void Curl_Phi_times_Curl_Phi(double *c, double weight, vector2 xi, vector2 eta,
	vector3 dChi_dx_s, vector3 dChi_dy_s, vector3 dChi_dx_t, vector3 dChi_dy_t);
#endif
