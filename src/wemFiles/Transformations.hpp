/**
 * @file Transformations.hpp
 *
 * @brief helper routines for quadrature
 */

#ifndef TRANSFORMATIONS_HPP
#define TRANSFORMATIONS_HPP

class Vector2;

// projection of unit square on element
Vector2 kappa(Vector2 a, Vector2 b, double h);


// calculate rotations for duffy trick
Vector2 tau(double b_x, double b_y, unsigned int CASE);
#endif
