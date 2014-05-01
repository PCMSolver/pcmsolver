#ifndef PURISIMAINTEGRATOR_HPP
#define PURISIMAINTEGRATOR_HPP

#include <iosfwd>

#include "Config.hpp"

#include "EigenPimpl.hpp"

#include "DiagonalIntegrator.hpp"

/*! \file PurisimaIntegrator.hpp
 *  \class PurisimaIntegrator
 *  \brief Implementation of diagonal elements of S and D using Purisima's formula
 *  \author Roberto Di Remigio
 *  \date 2014
 *
 *  Calculates the diagonal elements of S as:
 *  \f[
 *  	S_{ii} = 1.07\sqrt{\frac{4\pi}{a_i}}
 *  \f]
 *  while the diagonal elements of D are:
 *  \f[
 *  	D_{ii} = -(2\pi + \sum_{j\neq i} D_{ij}a_j)\frac{1}{a_i}
 *  \f]
 */

class PurisimaIntegrator : public DiagonalIntegrator
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

#endif // PURISIMAINTEGRATOR_HPP
