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

#define BOOST_TEST_MODULE NumericalQuadrature

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include <cmath>
#include <iomanip>
#include <limits>

#include "Config.hpp"

#include <Eigen/Dense>

#include <boost/numeric/quadrature/kronrodgauss.hpp>
#include <boost/numeric/quadrature/error_estimator.hpp>
#include <boost/lambda/lambda.hpp>

namespace quadrature = boost::numeric::quadrature;

/*! \class NumericalQuadrature
 *  \test \b NumericalQuadrature_x2sinx tests adaptive integration of x^2*sin(x) in [0, 5]
 */
BOOST_AUTO_TEST_CASE(x2sinx)
{
    // Define a functor to be passed to the integrator
    struct f {
        double operator()(double x) const { return (std::pow(x, 2) * std::sin(x)); }
    };
    double upper = 5.0;
    double lower = 0.0;
    double expected_result = -2.0 + 10.0 * std::sin(upper) - 23.0 * std::cos(upper);
    double actual_result;

    // declare the kernel
    quadrature::kronrod_gauss<15> integrator;

    // integrate x^2*sin(x) on [0,5]
    integrator(f(), lower, upper, actual_result);

    // display the result
    BOOST_TEST_MESSAGE("Integral(x^2*sin(x)) on [0, 5] is " << std::setprecision(
                           std::numeric_limits<long double>::digits10) << actual_result);
    BOOST_REQUIRE_CLOSE(expected_result, actual_result, 1.0e-12);
}
