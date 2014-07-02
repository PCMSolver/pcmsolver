/**
 * @file Vector3.hpp
 *
 * @brief 3d vector class
 */
#ifndef VECTOR3_HPP
#define VECTOR3_HPP

#include <math.h>

class Vector3
{
public:
    double x, y, z;

    Vector3 (double xi = 0.0, double yi = 0.0, double zi = 0.0) : x(xi), y(yi), z(zi) {}
};

// vector addition
inline Vector3 vector3Add(Vector3 a,Vector3 b);

// vector in-place addition
inline void vector3AddTo(Vector3* a,Vector3 b);

// vector subtration
inline Vector3 vector3Sub(Vector3 a,Vector3 b);

// cross product
inline Vector3 vector3Mul(Vector3 a,Vector3 b);

// multiplication by scalar
inline Vector3 vector3SMul(double a,Vector3 b);

// scalar product
inline double vector3Dot(Vector3 a,Vector3 b);

// euclidean norm
inline double vector3Norm(Vector3 a);

#endif
