/**
 * @file GaussSquare.hpp
 *
 * @brief defines the gauss quadrature formula for a 2D domain [0,1]x[0,1]
 */
#include <stdio.h>
#include <stdlib.h>
#include "Vector2.hpp"
#include "Quadrature.hpp"
#include "Cubature.hpp"
#include "GaussLegendre.hpp"
#include "GaussSquare.hpp"

// compute tensor product
void tensorRegelSquare(Cubature *Q, Quadrature *R) {
	unsigned int	k;

	Q->noP = R->nop*R->nop;
	Q->xi = (Vector2*) malloc(Q->noP*sizeof(Vector2));
	Q->weight  = (double*)  malloc(Q->noP*sizeof(double) );

	for (k=0; k<Q->noP; ++k) {
		Q->xi[k] = Vector2(R->xi[k/R->nop],R->xi[k%R->nop]);
		Q->weight[k]  = R->w[k/R->nop] * R->w[k%R->nop];
	}
	return;
}

// initialize the Gauss-Legendre polinomials and compute tensor product
void initGaussSquare(Cubature **Q,unsigned int g) {
	Quadrature      *R;
	unsigned int	k;

	// Fehler-Routine
	if (g > 15) {
		printf("g should be less equal 15\n");
		exit(0);
	}

	// Tensor-Produkte zusammenbasteln
	initGaussLegendre(&R,g);
	(*Q) = (Cubature*) malloc(g*sizeof(Cubature));
	for (k=0; k<g; k++)
		tensorRegelSquare(&(*Q)[k],&R[k]);
	freeGaussLegendre(&R,g);
	return;
}

// return memory to the system
void freeGaussSquare(Cubature **Q, unsigned int g) {
	unsigned int    k;

	for (k=0; k<g; ++k) {
		free((*Q)[k].xi);
		free((*Q)[k].weight);
	}
	free(*Q);
	return;
}
