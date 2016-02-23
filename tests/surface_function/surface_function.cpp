/* pcmsolver_copyright_start */
/*
 *     PCMSolver, an API for the Polarizable Continuum Model
 *     Copyright (C) 2013-2015 Roberto Di Remigio, Luca Frediani and contributors
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
 *     PCMSolver API, see: <http://pcmsolver.readthedocs.org/>
 */
/* pcmsolver_copyright_end */

#include "catch.hpp"

#include "Config.hpp"

#include <Eigen/Core>

#include "utils/SurfaceFunction.hpp"

TEST_CASE("SurfaceFunction", "[surface_function]")
{
    size_t nPoints = 1000;
    Eigen::VectorXd values1 = Eigen::VectorXd::Random(nPoints);
    Eigen::VectorXd values2 = Eigen::VectorXd::Random(nPoints);
    double * values1_ptr = values1.data();
    double * values2_ptr = values2.data();
    SurfaceFunction func1 = SurfaceFunction(nPoints, values1_ptr);
    SurfaceFunction func2 = SurfaceFunction(nPoints, values2_ptr);

    /*! \class SurfaceFunction
     *  \test \b SurfaceFunctionTest_addition tests addition of two SurfaceFunction
     */
    SECTION("SurfaceFunction addition")
    {
        SurfaceFunction addition = func1 + func2;
        REQUIRE( addition.nPoints() == nPoints );
        Eigen::VectorXd result = values1 + values2;
        for (size_t i = 0; i < nPoints; ++i) {
            REQUIRE( addition.value(i) == Approx(result(i)) );
        }
    }

    /*! \class SurfaceFunction
     *  \test \b SurfaceFunctionTest_subtraction tests subtraction of two SurfaceFunction
     */
    SECTION("SurfaceFunction subtraction")
    {
        SurfaceFunction subtraction = func1 - func2;
        REQUIRE( subtraction.nPoints() == nPoints );
        Eigen::VectorXd result = values1 - values2;
        for (size_t i = 0; i < nPoints; ++i) {
            REQUIRE( subtraction.value(i) == Approx(result(i)) );
        }
    }

    /*! \class SurfaceFunction
     *  \test \b SurfaceFunctionTest_multiply_by_scalar tests multiplication of a SurfaceFunction by a scalar
     */
    SECTION("SurfaceFunction multiply by scalar")
    {
        SurfaceFunction scaled1 = 2.5 * func1;
        func2 *= 0.5;
        REQUIRE( scaled1.nPoints() == nPoints );
        REQUIRE( func2.nPoints() == nPoints );
        Eigen::VectorXd result1 = 2.5 * values1;
        Eigen::VectorXd result2 = 0.5 * values2;
        for (size_t i = 0; i < nPoints; ++i) {
            REQUIRE( scaled1.value(i) == Approx(result1(i)) );
            REQUIRE( func2.value(i) == Approx(result2(i)) );
        }
    }

    /*! \class SurfaceFunction
     *  \test \b SurfaceFunctionTest_multiply tests dot product of two SurfaceFunction
     */
    SECTION("SurfaceFunction multiply")
    {
        double product = func1 * func2;
        double expected_product = values1.dot(values2);
        REQUIRE( product == Approx(expected_product) );
    }
}
