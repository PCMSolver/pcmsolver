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
 *     PCMSolver API, see: <https://repo.ctcc.no/projects/pcmsolver>
 */
/* pcmsolver_copyright_end */

#define BOOST_TEST_MODULE GreensFunctionAnisotropicLiquid

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include <cmath>
#include <iostream>

#include "Config.hpp"

#include <Eigen/Dense>

#include "AnalyticEvaluate.hpp"
#include "DerivativeTypes.hpp"
#include "AnisotropicLiquid.hpp"

struct AnisotropicLiquidTest {
    Eigen::Vector3d epsilon;
    Eigen::Vector3d euler;
    Eigen::Vector3d source, probe, sourceNormal, probeNormal;
    Eigen::Array4d result;
    AnisotropicLiquidTest() { SetUp(); }
    void SetUp() {
    	epsilon << 2.0, 80.0, 15.0;
    	euler << 6.0, 40.0, 15.0;
        source = Eigen::Vector3d::Random();
        sourceNormal = source + Eigen::Vector3d::Random();
        sourceNormal.normalize();
        probe = Eigen::Vector3d::Random();
        probeNormal = probe + Eigen::Vector3d::Random();
        probeNormal.normalize();
        result = analyticAnisotropicLiquid(epsilon, euler, sourceNormal, source, probeNormal, probe);
    }
};

BOOST_FIXTURE_TEST_SUITE(Aniso, AnisotropicLiquidTest)

/*! \class AnisotropicLiquid
 *  \test \b AnisotropicLiquidTest_numerical tests the numerical evaluation of the AnisotropicLiquid Green's function against analytical result
 */
BOOST_FIXTURE_TEST_CASE(numerical, AnisotropicLiquidTest)
{
    AnisotropicLiquid<double> gf(epsilon, euler);
    double value = result(0);
    double gf_value = gf.function(source, probe);
    BOOST_REQUIRE_CLOSE(value, gf_value, 1.0e-12);

    double derProbe = result(1);
    double gf_derProbe = gf.derivativeProbe(probeNormal, source, probe);
    BOOST_REQUIRE_CLOSE(derProbe, gf_derProbe, 1.0e-05);

    double derSource = result(2);
    double gf_derSource = gf.derivativeSource(sourceNormal, source, probe);
    BOOST_REQUIRE_CLOSE(derSource, gf_derSource, 1.0e-05);
}

/*! \class AnisotropicLiquid
 *  \test \b AnisotropicLiquidTest_directional_AD tests the automatic evaluation (directional derivative only)
 *  of the AnisotropicLiquid Green's function against analytical result
 */
BOOST_FIXTURE_TEST_CASE(directional_AD, AnisotropicLiquidTest)
{
    AnisotropicLiquid<AD_directional> gf(epsilon, euler);
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

/*! \class AnisotropicLiquid
 *  \test \b AnisotropicLiquidTest_gradient_AD tests the automatic evaluation (full gradient)
 *  of the AnisotropicLiquid Green's function against analytical result
 */
BOOST_FIXTURE_TEST_CASE(gradient_AD, AnisotropicLiquidTest)
{
    AnisotropicLiquid<AD_gradient> gf(epsilon, euler);
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

/*! \class AnisotropicLiquid
 *  \test \b AnisotropicLiquidTest_hessian_AD tests the automatic evaluation (full hessian)
 *  of the AnisotropicLiquid Green's function against analytical result
 */
BOOST_FIXTURE_TEST_CASE(hessian_AD, AnisotropicLiquidTest)
{
    AnisotropicLiquid<AD_hessian> gf(epsilon, euler);
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

BOOST_AUTO_TEST_SUITE_END()

// Define a uniform dielectric as an anisotropic dielectric and test the evaluation of the
// Green's function and derivatives against the analytic result for the uniform dielectric
struct AnisotropicLiquidUniformTest {
    Eigen::Vector3d epsilon;
    Eigen::Vector3d euler;
    Eigen::Vector3d source, probe, sourceNormal, probeNormal;
    Eigen::Array4d result;
    AnisotropicLiquidUniformTest() { SetUp(); }
    void SetUp() {
    	epsilon << 80.0, 80.0, 80.0;
    	euler << 0.0, 0.0, 0.0;
        source = Eigen::Vector3d::Random();
        sourceNormal = source + Eigen::Vector3d::Random();
        sourceNormal.normalize();
        probe = Eigen::Vector3d::Random();
        probeNormal = probe + Eigen::Vector3d::Random();
        probeNormal.normalize();
	result = analyticUniformDielectric(epsilon(0), sourceNormal, source, probeNormal, probe);
    }
};

BOOST_FIXTURE_TEST_SUITE(AnisotropicVsUniform, AnisotropicLiquidUniformTest)

/*! \class AnisotropicLiquid
 *  \test \b AnisotropicLiquidUniformTest_numerical tests the numerical evaluation of the AnisotropicLiquid Green's function against analytical result for a uniform dielectric
 */
BOOST_FIXTURE_TEST_CASE(numerical, AnisotropicLiquidUniformTest)
{
    AnisotropicLiquid<double> gf(epsilon, euler);
    double value = result(0);
    double gf_value = gf.function(source, probe);
    BOOST_REQUIRE_CLOSE(value, gf_value, 1.0e-12);

    double derProbe = result(1);
    double gf_derProbe = gf.derivativeProbe(probeNormal, source, probe);
    BOOST_REQUIRE_CLOSE(derProbe, gf_derProbe, 1.0e-05);

    double derSource = result(2);
    double gf_derSource = gf.derivativeSource(sourceNormal, source, probe);
    BOOST_REQUIRE_CLOSE(derSource, gf_derSource, 1.0e-05);
}

/*! \class AnisotropicLiquid
 *  \test \b AnisotropicLiquidUniformTest_directional_AD tests the automatic evaluation (directional derivative only)
 *  of the AnisotropicLiquid Green's function against analytical result for a uniform dielectric
 */
BOOST_FIXTURE_TEST_CASE(directional_AD, AnisotropicLiquidUniformTest)
{
    AnisotropicLiquid<AD_directional> gf(epsilon, euler);
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

/*! \class AnisotropicLiquid
 *  \test \b AnisotropicLiquidUniformTest_gradient_AD tests the automatic evaluation (full gradient)
 *  of the AnisotropicLiquid Green's function against analytical result for a uniform dielectric
 */
BOOST_FIXTURE_TEST_CASE(gradient_AD, AnisotropicLiquidUniformTest)
{
    AnisotropicLiquid<AD_gradient> gf(epsilon, euler);
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

/*! \class AnisotropicLiquid
 *  \test \b AnisotropicLiquidUniformTest_hessian_AD tests the automatic evaluation (full hessian)
 *  of the AnisotropicLiquid Green's function against analytical result for a uniform dielectric
 */
BOOST_FIXTURE_TEST_CASE(hessian_AD, AnisotropicLiquidUniformTest)
{
    AnisotropicLiquid<AD_hessian> gf(epsilon, euler);
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

BOOST_AUTO_TEST_SUITE_END()
