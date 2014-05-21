#ifndef GAUSS_LEGENDRE_HPP
#define GAUSS_LEGENDRE_HPP

/**
 * @file GaussLegendre.hpp
 *
 * @brief points xi and weights weight of the gauss quadrature formula on
 * interval [0,1]
 */

void initGaussLegendre(Quadrature **Q, unsigned int g);
void freeGaussLegendre(Quadrature **Q, unsigned int g);

#endif
