#ifndef SPHERE_HPP
#define SPHERE_HPP

#include <iosfwd>
#include <string>

#include "Config.hpp"

#include "EigenPimpl.hpp"

/*! \file Sphere.hpp
 *  \class Sphere
 *  \brief Class describing a sphere.
 *  \author Roberto Di Remigio
 *  \date 2011
 */

class Sphere
{
public:
    Sphere() {}
    Sphere(const Eigen::Vector3d & center, double radius,
           const std::string & colour = "Violet" )
        : center_(center), radius_(radius), colour_(colour) {}
    ~Sphere() {}
    double radius() const { return radius_; }
    void radius(double r) { radius_ = r; }
    const Eigen::Vector3d & center() const { return center_; }
    double center(int i) const { return center_(i); }
    void center(Eigen::Vector3d & coord) { center_ = coord; }
    const std::string & colour() const { return colour_; }
    void colour(std::string & col) { colour_ = col; }
    friend inline void swap(Sphere & left, Sphere & right);
    inline void swap(Sphere & other);
    /// Assignment operator.
    Sphere & operator=(Sphere other);
    friend std::ostream& operator<<(std::ostream & os, Sphere & sph) {
        return sph.printObject(os);
    }
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW /* See http://eigen.tuxfamily.org/dox/group__TopicStructHavingEigenMembers.html */
private:
    Eigen::Vector3d center_;
    double radius_;
    std::string colour_;
    std::ostream & printObject(std::ostream & os);
};

#endif // SPHERE_HPP
