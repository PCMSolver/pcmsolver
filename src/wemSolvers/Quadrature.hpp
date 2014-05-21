#ifndef QUADRATURE_HPP
#define QUADRATURE_HPP

/**
 * @file Quadrature.hpp
 */

/**
 * Cubature structure that contains the integration points and the weights on
 * an interval. 
 */

typedef struct {
	unsigned int nop; ///< number of integration points
	double *xi; ///< integration points
	double *w; ///< integration weights
} Quadrature;

#endif
