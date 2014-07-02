/**
 * @file Vector3.cpp
 *
 * @brief 3d vector class
 */

#include <math.h>
#include "Vector3.hpp"

// vector addition
inline Vector3 vector3Add(Vector3 a,Vector3 b)
{
    return Vector3(a.x+b.x, a.y+b.y, a.z+b.z);
}

// vector in-place addition
inline void vector3AddTo(Vector3* a,Vector3 b)
{
    a->x += b.x;
    a->y += b.y;
    a->z += b.z;
}

// vector subtration
inline Vector3 vector3Sub(Vector3 a,Vector3 b)
{
    return Vector3(a.x-b.x, a.y-b.y, a.z-b.z);
}

// cross product
inline Vector3 vector3Mul(Vector3 a,Vector3 b)
{
    return Vector3 (a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x);
}

// multiplication by scalar
inline Vector3 vector3SMul(double a,Vector3 b)
{
    return Vector3(a*b.x, a*b.y, a*b.z);
}

// scalar product
inline double vector3Dot(Vector3 a,Vector3 b)
{
    return(a.x*b.x + a.y*b.y + a.z*b.z);
}

// euclidean norm
inline double vector3Norm(Vector3 a)
{
    return(sqrt(a.x*a.x + a.y*a.y + a.z*a.z));
}

