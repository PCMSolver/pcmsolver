/*****************
 *  Constants.c  *
 *****************/
 
 
/*================================================*
 *  Hier werden alle Paramter gesetzt, die im     *
 *  Wavelet-Galerkin-Verfahren benoetigt werden.  *
 *================================================*/
 
/* Operatorordnung */
const double		op = -1;

/* Kompression: a > 1, 0 < b < 1, d < dp < td-op */

const double		a = 2.5;
const double		b = 0.001;
const double		dp = 2.5;


/* Quadratur */
const unsigned int	g_max = 10;			/* maximaler Quadraturgrad             */
const unsigned int	min_quadrature_level = 2;	/* minimales Quadraturlevel            */
const double		q = 1;				/* Unterteilungskonstante q > 0.25     */
const double		scaling_factor = 0.7071;	/* Groesse des relativen Umkreisradius */

/* Genaugikeit der iterativen Loesung */
const double		eps = 1e-10;

/* Konstante fuer Feldverlaengerung */
const unsigned int	delta = 10;

/* Genaugikeit bei Punktevergleich */
const double		tol = 1e-10;
