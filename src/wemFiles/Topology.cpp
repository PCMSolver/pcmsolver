#include "Topology.hpp"
#include <math.h>

#include "Vector3.hpp"

// calculates union K(d,r) = K(d1,r1) union K(d2,r2)
void unify(Vector3 *d, double *r, Vector3 d1, double r1, Vector3 d2, double r2) {
	Vector3 z(d1.x-d2.x, d1.y-d2.y, d1.z-d2.z);
	double norm;
	
	norm = sqrt(z.x*z.x+z.y*z.y+z.z*z.z);
	
	// K(d2,r2) \subset K(d1,r1)
	if (norm+r2 <= r1){
		*d = d1;
		*r = r1;
	} else if (norm+r1 <= r2) { 
    // K(d1,r1) \subset K(d2,r2)
		*d = d2;
		*r = r2;
   } else { 
    // union is not a sphere
		(*d).x = 0.5*(d1.x+d2.x+(r1-r2)/norm*z.x);
		(*d).y = 0.5*(d1.y+d2.y+(r1-r2)/norm*z.y);
		(*d).z = 0.5*(d1.z+d2.z+(r1-r2)/norm*z.z);
		*r = 0.5*(r1+r2+norm);
	}
	return;
}

