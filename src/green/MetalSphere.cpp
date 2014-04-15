#include "MetalSphere.hpp"

#include <cmath>
#include <ostream>
#include <stdexcept>

#include "Config.hpp"

#include "EigenPimpl.hpp"
#include "TaylorPimpl.hpp"

#include "DerivativeTypes.hpp"
#include "GreensFunction.hpp"
#include "IGreensFunction.hpp"

extern "C" void greens_function(const double * epssol, const double * epsre, const double * epsim,
                     const double * sphRadius, const double * ps, const double * p1, const double * p2,
                     double * greenre, double * greenim);

double MetalSphere::operator()(double * source, double * probe) const 
{
    // Calculation of the value of the Green's Function
    double epsre, epsim;
    double greenre, greenim;
    double point1[3], point2[3], sphere[3];
    point1[0] = source[0];
    point1[1] = source[1];
    point1[2] = source[2];
    point2[0] =  probe[0];
    point2[1] =  probe[1];
    point2[2] =  probe[2];
    sphere[0] = sphPosition_(0);
    sphere[1] = sphPosition_(1);
    sphere[2] = sphPosition_(2);
    epsre = epsMetal_.real();
    epsim = epsMetal_.imag();
    // Call the Fortran subroutine
    greens_function(&epsSolvent_, &epsre, &epsim, &sphRadius_, sphere, point1, point2,
                &greenre, &greenim);

    return greenre;
}

double MetalSphere::derivative(const Eigen::Vector3d & direction,
                                        const Eigen::Vector3d & p1, const Eigen::Vector3d & p2) const
{
        return epsSolvent_ * (this->derivativeProbe(direction, p1, p2));
}

std::ostream & MetalSphere::printObject(std::ostream & os)
{
    os << "Green's function type: metal sphere" << std::endl;
    os << "Permittivity (real part)      = " << epsMetal_.real() << std::endl;
    os << "Permittivity (imaginary part) = " << epsMetal_.imag() << std::endl;
    os << "Sphere position               = " << sphPosition_ << std::endl;
    os << "Sphere radius                 = " << sphRadius_;
    return os;
}
