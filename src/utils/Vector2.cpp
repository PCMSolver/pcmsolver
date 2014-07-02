/**
 * @file Vector2.cpp
 *
 * @brief 2d vector class
 */

#include <math.h>
#include "Vector2.hpp"

// vector addition
inline Vector2 vector2Add(Vector2 a, Vector2 b)
{
    return Vector2 (a.x+b.x, a.y+b.y);
}

// vector subtration
inline Vector2 vector2Sub(Vector2 a, Vector2 b)
{
    return Vector2(a.x-b.x, a.y-b.y);
}

// multiplication by scalar
inline Vector2 vector2SMul(double a,Vector2 b)
{
    return Vector2 (a*b.x, a*b.y);
}

// scalar product
inline double vector2Dot(Vector2 a,Vector2 b)
{
    return(a.x*b.x + a.y*b.y);
}

// euclidean norm
inline double vector2Norm(Vector2 a)
{
    return(sqrt(a.x*a.x + a.y*a.y));
}
