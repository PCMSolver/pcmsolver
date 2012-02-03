/***************
 *  IntKon4.c  *
 ***************/


#include <math.h>
#include "intvector.h"
#include "vector2.h"
#include "vector3.h"
#include "basis.h"
#include "trafos.h"
#include "kern.h"
#include "cubature.h"
#include "interpolate.h"
#include "intkon4.h"


void 
IntKon4(c, element1, element2, ind_s, ind_t, Q, P, M, SingleLayer, DoubleLayer)
/* GEMEINSAME ECKE IM NULLPUNKT -> modifiziertes Skalarprodukt */
	double         *c;
	element        *element1, *element2;
	unsigned int    ind_s, ind_t;
	cubature       *Q;
	vector3     ****P;
	unsigned int    M;
	double          SingleLayer(), DoubleLayer();
{
	unsigned int    i, j;
	double          w, d1, d2, d3;
	vector2         s, t, xi, eta, a;
	vector3         x1, n_x1, x2, n_x2, y1, n_y1, y2, n_y2, z, n_z;
	double          h = 1. / (1 << element1->level);

	c[0] = c[1] = c[2] = 0;
	s = vector2_make(h * element1->index_s, h * element1->index_t);
	t = vector2_make(h * element2->index_s, h * element2->index_t);
	for (i = 0; i < Q->nop; i++) {
		xi = Q->xi[i];
		w = pow(xi.x, 3) * Q->w[i];
		xi.y *= xi.x;

		a = Kappa(s, Tau(xi.x, xi.y, ind_s), h);
		x1 = Chi(a, P[element1->patch], M);
		n_x1 = n_Chi(a, P[element1->patch], M);

		a = Kappa(s, Tau(xi.y, xi.x, ind_s), h);
		x2 = Chi(a, P[element1->patch], M);
		n_x2 = n_Chi(a, P[element1->patch], M);

		a = Kappa(t, Tau(xi.x, xi.y, ind_t), h);
		y1 = Chi(a, P[element2->patch], M);
		n_y1 = n_Chi(a, P[element2->patch], M);

		a = Kappa(t, Tau(xi.y, xi.x, ind_t), h);
		y2 = Chi(a, P[element2->patch], M);
		n_y2 = n_Chi(a, P[element2->patch], M);

		for (j = 0; j < Q->nop; j++) {
			eta = vector2_Smul(xi.x, Q->xi[j]);

			a = Kappa(t, Tau(eta.x, eta.y, ind_t), h);
			z = Chi(a, P[element2->patch], M);
			n_z = n_Chi(a, P[element2->patch], M);

			d1 = SingleLayer(x1, z) + SingleLayer(x2, z);
			d2 = DoubleLayer(x1, z, n_z) + DoubleLayer(x2, z, n_z);
			d3 = DoubleLayer(z, x1, n_x1) + DoubleLayer(z, x2, n_x2);

			a = Kappa(s, Tau(eta.x, eta.y, ind_s), h);
			z = Chi(a, P[element1->patch], M);
			n_z = n_Chi(a, P[element1->patch], M);

			d1 += SingleLayer(z, y1) + SingleLayer(z, y2);
			d2 += DoubleLayer(z, y1, n_y1) + DoubleLayer(z, y2, n_y2);
			d3 += DoubleLayer(y1, z, n_z) + DoubleLayer(y2, z, n_z);

			c[0] += w * Q->w[j] * d1;
			c[1] += w * Q->w[j] * d2;
			c[2] += w * Q->w[j] * d3;
		}
	}

	c[0] *= h * h;		/* L^2-normiert */
	c[1] *= h * h;
	c[2] *= h * h;
	return;
}
