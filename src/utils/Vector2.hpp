/**
 * @file Vector2.hpp
 *
 * @brief 2d vector class
 */

#ifndef VECTOR2_HPP
#define VECTOR2_HPP

#include <math.h>

/// type definition
class Vector2
{
public:
    double x, y;

    Vector2(double xi = 0.0, double yi = 0.0):x(xi),y(yi) {}

};

// vector addition
inline Vector2 vector2Add(Vector2 a, Vector2 b);

// vector subtration
inline Vector2 vector2Sub(Vector2 a, Vector2 b);

// multiplication by scalar
inline Vector2 vector2SMul(double a,Vector2 b);

// scalar product
inline double vector2Dot(Vector2 a,Vector2 b);

// euclidean norm
inline double vector2Norm(Vector2 a);
#endif
