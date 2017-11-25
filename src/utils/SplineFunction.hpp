/*
 * PCMSolver, an API for the Polarizable Continuum Model
 * Copyright (C) 2017 Roberto Di Remigio, Luca Frediani and collaborators.
 *
 * This file is part of PCMSolver.
 *
 * PCMSolver is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * PCMSolver is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with PCMSolver.  If not, see <http://www.gnu.org/licenses/>.
 *
 * For information on the complete list of contributors to the
 * PCMSolver API, see: <http://pcmsolver.readthedocs.io/>
 */

#pragma once

#include <algorithm>

#include "Config.hpp"

#include <Eigen/Core>
#include <unsupported/Eigen/Splines>

/*! \file SplineFunction.hpp
 *  \class SplineFunction
 *  \brief Spline interpolation of a function
 *  \author Roberto Di Remigio
 *  \date 2015
 *
 *  Taken from StackOverflow http://stackoverflow.com/a/29825204/2528668
 */

class SplineFunction __final {
private:
  typedef Eigen::Spline<double, 1> CubicSpline;

public:
  /*! \brief Constructor from abscissa and function values
   *  \param[in] x vector with abscissa values
   *  \param[in] y vector with function values
   *
   *  Interpolation happens in the initialization of the spline_ data
   *  member. We use std::min to set the degree in order to ensure that no
   *  more than a cubic spline is used, but short vectors are still accepted.
   */
  SplineFunction(const Eigen::VectorXd & x, const Eigen::VectorXd & y)
      : xMin_(x.minCoeff()),
        xMax_(x.maxCoeff()),
        spline_(Eigen::SplineFitting<CubicSpline>::Interpolate(
            y.transpose(),
            std::min<int>(x.rows() - 1, 3),
            fitVector(x))) {}
  /*! \brief Evaluate spline at given point
   *  \param[in] x evaluation point
   */
  double operator()(double x) const { return spline_(fitScalar(x))(0); }

private:
  double xMin_;
  double xMax_;
  CubicSpline spline_;
  /*! \brief Scale a scalar value to [0, 1] interval
   *  \param[in] x value to be scaled
   *
   *  val is defined in [min, max] This function returns the variable
   *  as scaled to fit in the [0, 1] interval.
   */
  double fitScalar(double x) const { return (x - xMin_) / (xMax_ - xMin_); }
/*! \brief Scale a vector to [0, 1] interval
 *  \param[in] x_vec a vector
 *
 *  Given column vector, returns a *row* vector with the values scaled
 *  to fit a preset interval. The interval is embedded in the callable object.
 */
#ifdef HAS_CXX11_LAMBDA
  Eigen::RowVectorXd fitVector(const Eigen::VectorXd & x_vec) const {
    return x_vec.unaryExpr([this](double x) -> double { return fitScalar(x); })
        .transpose();
  }
#else  /* HAS_CXX11_LAMBDA */
  Eigen::RowVectorXd fitVector(const Eigen::VectorXd & x_vec) const {
    pcm::function<double(double)> fit =
        pcm::bind(&SplineFunction::fitScalar, this, pcm::_1);
    return x_vec.unaryExpr(fit).transpose();
  }
#endif /*HAS_CXX11_LAMBDA */
};
