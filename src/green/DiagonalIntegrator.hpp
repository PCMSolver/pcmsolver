#ifndef DIAGONALINTEGRATOR_HPP
#define DIAGONALINTEGRATOR_HPP

#include <iosfwd>

#include "Config.hpp"

#include "EigenPimpl.hpp"

/*! \file DiagonalIntegrator.hpp
 *  \class DiagonalIntegrator
 *  \brief Abstract Base Class for implementation of diagonal elements of S and D
 *  \author Roberto Di Remigio
 *  \date 2014
 *
 *  This class encapsulates the calculation of the diagonal elements of the S and D
 *  operators needed to set up the PCM equations.
 *  Based on the ideas of the Strategy Pattern.
 */

class DiagonalIntegrator
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

#endif // DIAGONALINTEGRATOR_HPP
