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

#define BOOST_TEST_MODULE GreensFunctionUniformDielectric

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include <cmath>
#include <iostream>

#include "Config.hpp"

#include <Eigen/Dense>

#include "DerivativeTypes.hpp"
#include "UniformDielectric.hpp"

struct UniformDielectricTest {
    Eigen::Array4d analyticEvaluate(double eps, const Eigen::Vector3d & spNormal,
                                    const Eigen::Vector3d & sp,
                                    const Eigen::Vector3d & ppNormal, const Eigen::Vector3d & pp) {
        Eigen::Array4d result = Eigen::Array4d::Zero();
        double distance = (sp - pp).norm();
        double distance_3 = std::pow(distance, 3);
        double distance_5 = std::pow(distance, 5);

        // Value of the function
        result(0) = 1.0 / (eps * distance);
        // Value of the directional derivative wrt probe
        result(1) = (sp - pp).dot(ppNormal) / (eps * distance_3);
        // Directional derivative wrt source
        result(2) = - (sp - pp).dot(spNormal) / (eps * distance_3);
        // Value of the Hessian
        result(3) = spNormal.dot(ppNormal) / (eps * distance_3) - 3 * ((
                        sp - pp).dot(spNormal))*((sp - pp).dot(
                                    ppNormal)) / (eps * distance_5);

        return result;
    }
    double epsilon;
    Eigen::Vector3d source, probe, sourceNormal, probeNormal;
    Eigen::Array4d result;
    UniformDielectricTest() { SetUp(); }
    void SetUp() {
        epsilon = 60.0;
        source = Eigen::Vector3d::Random();
        sourceNormal = source + Eigen::Vector3d::Random();
        sourceNormal.normalize();
        probe = Eigen::Vector3d::Random();
        probeNormal = probe + Eigen::Vector3d::Random();
        probeNormal.normalize();
        result = analyticEvaluate(epsilon, sourceNormal, source, probeNormal, probe);
    }
};

/*! \class UniformDielectric
 *  \test \b UniformDielectricTest_numerical tests the numerical evaluation of the UniformDielectric Green's function against analytical result
 */
BOOST_FIXTURE_TEST_CASE(numerical, UniformDielectricTest)
{
    UniformDielectric<double> gf(epsilon);
    double value = result(0);
    double gf_value = gf.function(source, probe);
    BOOST_REQUIRE_CLOSE(value, gf_value, 1.0e-12);

    double derProbe = result(1);
    double gf_derProbe = gf.derivativeProbe(probeNormal, source, probe);
    BOOST_REQUIRE_CLOSE(derProbe, gf_derProbe, 1.0e-06);

    double derSource = result(2);
    double gf_derSource = gf.derivativeSource(sourceNormal, source, probe);
    BOOST_REQUIRE_CLOSE(derSource, gf_derSource, 1.0e-06);
}

/*! \class UniformDielectric
 *  \test \b UniformDielectricTest_directional_AD tests the automatic evaluation (directional derivative only)
 *  of the UniformDielectric Green's function against analytical result
 */
BOOST_FIXTURE_TEST_CASE(directional_AD, UniformDielectricTest)
{
    Eigen::Array4d result = analyticEvaluate(epsilon, sourceNormal, source, probeNormal,
                            probe);

    UniformDielectric<AD_directional> gf(epsilon);
    double value = result(0);
    double gf_value = gf.function(source, probe);
    BOOST_REQUIRE_CLOSE(value, gf_value, 1.0e-12);

    double derProbe = result(1);
    double gf_derProbe = gf.derivativeProbe(probeNormal, source, probe);
    BOOST_REQUIRE_CLOSE(derProbe, gf_derProbe, 1.0e-12);

    double derSource = result(2);
    double gf_derSource = gf.derivativeSource(sourceNormal, source, probe);
    BOOST_REQUIRE_CLOSE(derSource, gf_derSource, 1.0e-12);
}

/*! \class UniformDielectric
 *  \test \b UniformDielectricTest_gradient_AD tests the automatic evaluation (full gradient)
 *  of the UniformDielectric Green's function against analytical result
 */
BOOST_FIXTURE_TEST_CASE(gradient_AD, UniformDielectricTest)
{
    Eigen::Array4d result = analyticEvaluate(epsilon, sourceNormal, source, probeNormal,
                            probe);

    UniformDielectric<AD_gradient> gf(epsilon);
    double value = result(0);
    double gf_value = gf.function(source, probe);
    BOOST_REQUIRE_CLOSE(value, gf_value, 1.0e-12);

    double derProbe = result(1);
    double gf_derProbe = gf.derivativeProbe(probeNormal, source, probe);
    BOOST_REQUIRE_CLOSE(derProbe, gf_derProbe, 1.0e-12);

    double derSource = result(2);
    double gf_derSource = gf.derivativeSource(sourceNormal, source, probe);
    BOOST_REQUIRE_CLOSE(derSource, gf_derSource, 1.0e-12);
}

/*! \class UniformDielectric
 *  \test \b UniformDielectricTest_hessian_AD tests the automatic evaluation (full hessian)
 *  of the UniformDielectric Green's function against analytical result
 */
BOOST_FIXTURE_TEST_CASE(hessian_AD, UniformDielectricTest)
{
    Eigen::Array4d result = analyticEvaluate(epsilon, sourceNormal, source, probeNormal,
                            probe);

    UniformDielectric<AD_hessian> gf(epsilon);
    double value = result(0);
    double gf_value = gf.function(source, probe);
    BOOST_REQUIRE_CLOSE(value, gf_value, 1.0e-12);

    double derProbe = result(1);
    double gf_derProbe = gf.derivativeProbe(probeNormal, source, probe);
    BOOST_REQUIRE_CLOSE(derProbe, gf_derProbe, 1.0e-12);

    double derSource = result(2);
    double gf_derSource = gf.derivativeSource(sourceNormal, source, probe);
    BOOST_REQUIRE_CLOSE(derSource, gf_derSource, 1.0e-12);

    /*	double hessian = result(4);
    	double gf_hessian = gf.hessian(sourceNormal, source, probeNormal, probe);
    	BOOST_REQUIRE_CLOSE(hessian, gf_hessian, 1.0e-12);*/
}
