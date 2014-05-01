#ifndef COLLOCATIONINTEGRATOR_HPP
#define COLLOCATIONINTEGRATOR_HPP

#include <iosfwd>

#include "Config.hpp"

#include "EigenPimpl.hpp"

#include "DiagonalIntegrator.hpp"

/*! \file CollocationIntegrator.hpp
 *  \class CollocationIntegrator
 *  \brief Implementation of diagonal elements of S and D using approximate collocation
 *  \author Roberto Di Remigio
 *  \date 2014
 *
 *  Calculates the diagonal elements of S as:
 *  \f[
 *  	S_{ii} = 1.07\sqrt{\frac{4\pi}{a_i}}
 *  \f]
 *  while the diagonal elements of D are:
 *  \f[
 *  	D_{ii} = -1.07\sqrt{\frac{\pi}{a_i}} \frac{1}{R_I}
 *  \f]
 */

class CollocationIntegrator : public DiagonalIntegrator
{
public:
    /*! S operator diagonal elements
     */
    virtual double operator()(const Eigen::Vector3d & sp,
                              const Eigen::Vector3d & pp) = 0;
    /*! D operator diagonal elements
     */
    virtual double operator()(const Eigen::Vector3d & direction,
                              const Eigen::Vector3d & sp, const Eigen::Vector3d & pp) = 0;
};

#endif // COLLOCATIONINTEGRATOR_HPP
