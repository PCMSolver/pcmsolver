/* pcmsolver_copyright_start */
/*
 *     PCMSolver, an API for the Polarizable Continuum Model
 *     Copyright (C) 2013 Roberto Di Remigio, Luca Frediani and contributors
 *
 *     This file is part of PCMSolver.
 *
 *     PCMSolver is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU Lesser General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *
 *     PCMSolver is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU Lesser General Public License for more details.
 *
 *     You should have received a copy of the GNU Lesser General Public License
 *     along with PCMSolver.  If not, see <http://www.gnu.org/licenses/>.
 *
 *     For information on the complete list of contributors to the
 *     PCMSolver API, see: <http://pcmsolver.github.io/pcmsolver-doc>
 */
/* pcmsolver_copyright_end */

#ifndef SPLINEFUNCTION_HPP
#define SPLINEFUNCTION_HPP

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
   
/*! Scale abscissa value to [0, 1] */
inline double scale(double x) {
	return (x - xMin_) / (xMax_ - xMin_);
}

/*! Scale abscissa values to [0, 1] */
inline Eigen::RowVectorXd scale(const Eigen::VectorXd & x_vec) {
    return x_vec.unaryExpr([](double x) -> double { return scale(x); }).transpose();
}

class SplineFunction
{
private:
    typedef Eigen::Spline<double, 1, 3> SplineFit;
public:
    /*! \param[in] x vector with abscissa values 
     *  \param[in] y vector with function values
     */
    SplineFunction(const Eigen::VectorXd & x, const Eigen::VectorXd & y) 
	    : xMin_(x.minCoeff()), xMax_(x.maxCoeff()),
	    spline_(Eigen::SplineFitting<SplineFit>::Interpolate(y.transpose(), 
		      std::min<int>(x.rows() - 1, 3), scale(x))
	           )
	{}
    /*! \brief Evaluate spline at given point
     *  \param[in] x evaluation point
     */
    double operator()(double x) const {
	return spline_(scale(x))(0);
    }
private:
    double xMin_;
    double xMax_;
    SplineFit spline_;
};
#endif // SPLINEFUNCTION_HPP
