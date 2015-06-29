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

#define BOOST_TEST_MODULE TanhOneLayer

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <tuple>

#include "Config.hpp"

#include <Eigen/Core>

#include "TanhDiffuse.hpp"

struct TanhOneLayerTest {
    TanhOneLayerTest() { SetUp(); }
    double eps1, eps2, center, width;
    double point;
    void SetUp() {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(0.0, 100.0);
        point = dis(gen);
        eps1 = 78.39;
        eps2 = 1.0;
        center = 50.0;
        width = 10.0;
    }
};

/*! \class TanhDiffuse
 *  \test \b TanhOneLayerTest tests the evaluation of the one layer hyperbolic tangent profile against the analytic value
 */
BOOST_FIXTURE_TEST_CASE(onelayer, TanhOneLayerTest)
{
    TanhDiffuse diffuse(eps1, eps2, width, center);
    double value = 0.0, deriv = 0.0;
    std::tie(value, deriv) = diffuse(point);

    double auto_v = diffuse.autovalue(point);
    BOOST_TEST_MESSAGE("value    = " << std::setprecision(std::numeric_limits<long double>::digits10) << value);
    BOOST_TEST_MESSAGE("auto_v = " << std::setprecision(std::numeric_limits<long double>::digits10) << auto_v);
    BOOST_CHECK_CLOSE(auto_v, value, 1.0e-12);

    double auto_d = diffuse.autoderivative(point);
    BOOST_TEST_MESSAGE("deriv    = " << std::setprecision(std::numeric_limits<long double>::digits10) << deriv);
    BOOST_TEST_MESSAGE("auto_d = " << std::setprecision(std::numeric_limits<long double>::digits10) << auto_d);
    BOOST_CHECK_CLOSE(auto_d, deriv, 1.0e-12);
}
