#ifndef GAUSS_SQUARE_HPP
#define GAUSS_SQUARE_HPP
/**
 * @file GaussSquare.hpp
 *
 * @brief defines the gauss quadrature formula for a 2D domain [0,1]x[0,1]
 */
 
void initGaussSquare(Cubature **Q, unsigned int g);
void freeGaussSquare(Cubature **Q, unsigned int g);

#endif
