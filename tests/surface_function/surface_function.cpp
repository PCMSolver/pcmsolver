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

#define BOOST_TEST_MODULE SurfaceFunction

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include "Config.hpp"

#include <Eigen/Core>

#include "SurfaceFunction.hpp"

struct SurfaceFunctionTest {
    int nPoints;
    Eigen::VectorXd values1;
    Eigen::VectorXd values2;
    SurfaceFunction func1;
    SurfaceFunction func2;
    SurfaceFunctionTest() { SetUp(); }
    void SetUp() {
        nPoints = 1000;
        values1.resize(nPoints);
        values2.resize(nPoints);
        values1 = Eigen::VectorXd::Random(nPoints);
        values2 = Eigen::VectorXd::Random(nPoints);
        double * values1_ptr = values1.data();
        double * values2_ptr = values2.data();
        func1 = SurfaceFunction("TestFunction1", nPoints, values1_ptr);
        func2 = SurfaceFunction("TestFunction2", nPoints, values2_ptr);
    }
};

/*! \class SurfaceFunction
 *  \test \b SurfaceFunctionTest_addition tests addition of two SurfaceFunction
 */
BOOST_FIXTURE_TEST_CASE(addition, SurfaceFunctionTest)
{
    SurfaceFunction addition = func1 + func2;
    BOOST_REQUIRE_EQUAL(nPoints, addition.nPoints());
    Eigen::VectorXd result(nPoints);
    result = values1 + values2;
    for (int i = 0; i < nPoints; ++i) {
        BOOST_REQUIRE_CLOSE(result(i), addition.value(i), 1.0e-12);
    }
}

/*! \class SurfaceFunction
 *  \test \b SurfaceFunctionTest_subtraction tests subtraction of two SurfaceFunction
 */
BOOST_FIXTURE_TEST_CASE(subtraction, SurfaceFunctionTest)
{
    SurfaceFunction subtraction = func1 - func2;
    BOOST_REQUIRE_EQUAL(nPoints, subtraction.nPoints());
    Eigen::VectorXd result(nPoints);
    result = values1 - values2;
    for (int i = 0; i < nPoints; ++i) {
        BOOST_REQUIRE_CLOSE(result(i), subtraction.value(i), 1.0e-12);
    }
}

/*! \class SurfaceFunction
 *  \test \b SurfaceFunctionTest_multiply_by_scalar tests multiplication of a SurfaceFunction by a scalar
 */
BOOST_FIXTURE_TEST_CASE(multiply_by_scalar, SurfaceFunctionTest)
{
    SurfaceFunction scaled1 = 2.5 * func1;
    func2 *= 0.5;
    BOOST_REQUIRE_EQUAL(nPoints, scaled1.nPoints());
    BOOST_REQUIRE_EQUAL(nPoints, func2.nPoints());
    Eigen::VectorXd result1(nPoints), result2(nPoints);
    result1 = 2.5 * values1;
    result2 = 0.5 * values2;
    for (int i = 0; i < nPoints; ++i) {
        BOOST_REQUIRE_CLOSE(result1(i), scaled1.value(i), 1.0e-12);
        BOOST_REQUIRE_CLOSE(result2(i), func2.value(i), 1.0e-12);
    }
}

/*! \class SurfaceFunction
 *  \test \b SurfaceFunctionTest_multiply tests dot product of two SurfaceFunction
 */
BOOST_FIXTURE_TEST_CASE(multiply, SurfaceFunctionTest)
{
    double product = func1 * func2;
    double expected_product = values1.dot(values2);
    BOOST_REQUIRE_CLOSE(expected_product, product, 1.0e-12);
}
