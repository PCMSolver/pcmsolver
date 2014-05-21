/**
 * @file Transformations.cpp
 *
 * @brief helper routines for quadrature
 */

#include "Transformations.hpp"
#include "Vector2.hpp"

// projection of unit square on element
Vector2 kappa(Vector2 a,Vector2 b,double h) {
	return Vector2(a.x + h*b.x,a.y + h*b.y);
}


// calculate rotations for duffy trick
Vector2 tau(double b_x, double b_y, unsigned int CASE) {
	switch (CASE) {
		case 1:  return(Vector2(1-b_y,  b_x));	// pi/2   rotation 
		case 2:  return(Vector2(1-b_x,1-b_y));	// pi     rotation
		case 3:  return(Vector2(  b_y,1-b_x));	// 3*pi/2 rotation 
		default: return(Vector2(  b_x,  b_y));	// identity
	}
}
