#ifndef CUBATURE_HPP
#define CUBATURE_HPP

/**
 * @file Cubature.hpp
 */

/**
 * Cubature structure that contains the integration points and the weights on
 * the reference domain
 */

#include "Vector2.hpp"

typedef struct {
	unsigned int noP; ///< number of integration points
	Vector2 *xi; ///< integration points
	double *weight; ///< integration weights 
} Cubature;

#endif
