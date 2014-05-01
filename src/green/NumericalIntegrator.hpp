#ifndef NUMERICALINTEGRATOR_HPP
#define NUMERICALINTEGRATOR_HPP

#include <iosfwd>

#include "Config.hpp"

#include "EigenPimpl.hpp"

#include <boost/numeric/quadrature/kronrodgauss.hpp>
#include <boost/numeric/quadrature/error_estimator.hpp>

#include "DiagonalIntegrator.hpp"

/*! \file NumericalIntegrator.hpp
 *  \class NumericalIntegrator
 *  \brief Implementation of diagonal elements of S and D using numerical integration
 *  \author Roberto Di Remigio
 *  \date 2014
 *
 *  Calculates the diagonal elements of S and D by collocation, using numerical
 *  integration.
 */

class NumericalIntegrator : public DiagonalIntegrator
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

#endif // NUMERICALINTEGRATOR_HPP
