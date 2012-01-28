/*****************
 *  Constants.h  *
 *****************/
 
 
/*================================================*
 *  Hier werden alle Paramter gesetzt, die im     *
 *  Wavelet-Galerkin-Verfahren benoetigt werden.  *
 *================================================*/
 
/* Operatorordnung */
extern const double		op;

/* Konstanten bezueglich der verwendeteten Wavelet-Basis */
extern const unsigned int	td;
extern const unsigned int	minLevel;

/* Kompression: a > 1, 0 < b < 1, d < dp < td-op */
extern const double		a;
extern const double		b;
extern const double		dp;

/* Quadratur */
extern const unsigned int	g_max;			/* maximaler Quadraturgrad             */
extern const unsigned int	min_quadrature_level;	/* minimales Quadraturlevel            */
extern const double		q;			/* Unterteilungskonstante q > 0.25     */
extern const double		scaling_factor;		/* Groesse des relativen Umkreisradius */

/* Genaugikeit der iterativen Loesung */
extern const double		eps;

/* Konstante fuer Feldverlaengerung */
extern const unsigned int	delta;

/* Genaugikeit bei Punktevergleich */
extern const double            	tol;
