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

#define BOOST_TEST_MODULE GreensFunctionSphericalSharp

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include <cmath>
#include <iostream>

#include "Config.hpp"

#include <Eigen/Dense>

#include "CollocationIntegrator.hpp"
#include "AnalyticEvaluate.hpp"
#include "SphericalSharp.hpp"
#include "DerivativeTypes.hpp"

struct SphericalSharpTest {
    double epsilonSolvent;
    double eps;
    double radius;
    Eigen::Vector3d origin;
    Eigen::Vector3d source, probe, sourceNormal, probeNormal;
    Eigen::Array4d result;
    SphericalSharpTest() { SetUp(); }
    void SetUp() {
        epsilonSolvent = 1.0;
        eps = 78.39;
        radius = 50.0;
        origin = Eigen::Vector3d::Zero();
        source = 100 * Eigen::Vector3d::Random();
        sourceNormal = source + Eigen::Vector3d::Random();
        sourceNormal.normalize();
        probe = 100 * Eigen::Vector3d::Random();
        probeNormal = probe + Eigen::Vector3d::Random();
        probeNormal.normalize();
        result = analyticSphericalSharp(eps, epsilonSolvent, radius, origin, sourceNormal, source, probeNormal, probe);
    }
};

BOOST_FIXTURE_TEST_SUITE(SharpDielectricSphere, SphericalSharpTest)

/*! \class SphericalSharp
 *  \test \b SphericalSharpTest_numerical tests the numerical evaluation of the SphericalSharp Green's function against analytical result
 */
BOOST_FIXTURE_TEST_CASE(numerical, SphericalSharpTest)
{
    SphericalSharp<Numerical, CollocationIntegrator> gf(eps, epsilonSolvent, radius, origin);
    double analytic = result(0);
    double gf_value = gf.kernelS(source, probe);
    BOOST_REQUIRE_CLOSE(analytic, gf_value, 1.0e-12);

    double analytic_derProbe = result(1);
    double gf_derProbe = gf.derivativeProbe(probeNormal, source, probe);
    BOOST_REQUIRE_CLOSE(analytic_derProbe, gf_derProbe, 1.0e-05);

    double analytic_derSource = result(2);
    double gf_derSource = gf.derivativeSource(sourceNormal, source, probe);
    BOOST_REQUIRE_CLOSE(analytic_derSource, gf_derSource, 1.0e-05);
}

/*! \class SphericalSharp
 *  \test \b SphericalSharpTest_directional_AD tests the automatic evaluation (directional derivative only)
 *  of the SphericalSharp Green's function against analytical result
 */
BOOST_FIXTURE_TEST_CASE(directional_AD, SphericalSharpTest)
{
    SphericalSharp<AD_directional, CollocationIntegrator> gf(eps, epsilonSolvent, radius, origin);
    double analytic = result(0);
    double gf_value = gf.kernelS(source, probe);
    BOOST_REQUIRE_CLOSE(analytic, gf_value, 1.0e-12);

    double analytic_derProbe = result(1);
    double gf_derProbe = gf.derivativeProbe(probeNormal, source, probe);
    BOOST_REQUIRE_CLOSE(analytic_derProbe, gf_derProbe, 1.0e-12);

    double analytic_derSource = result(2);
    double gf_derSource = gf.derivativeSource(sourceNormal, source, probe);
    BOOST_REQUIRE_CLOSE(analytic_derSource, gf_derSource, 1.0e-12);
}

/*! \class SphericalSharp
 *  \test \b SphericalSharpTest_gradient_AD tests the automatic evaluation (full gradient)
 *  of the SphericalSharp Green's function against analytical result
 */
BOOST_FIXTURE_TEST_CASE(gradient_AD, SphericalSharpTest)
{
    SphericalSharp<AD_gradient, CollocationIntegrator> gf(eps, epsilonSolvent, radius, origin);
    double analytic = result(0);
    double gf_value = gf.kernelS(source, probe);
    BOOST_REQUIRE_CLOSE(analytic, gf_value, 1.0e-12);

    double analytic_derProbe = result(1);
    double gf_derProbe = gf.derivativeProbe(probeNormal, source, probe);
    BOOST_REQUIRE_CLOSE(analytic_derProbe, gf_derProbe, 1.0e-12);

    double analytic_derSource = result(2);
    double gf_derSource = gf.derivativeSource(sourceNormal, source, probe);
    BOOST_REQUIRE_CLOSE(analytic_derSource, gf_derSource, 1.0e-12);
}

/*! \class SphericalSharp
 *  \test \b SphericalSharpTest_hessian_AD tests the automatic evaluation (full hessian)
 *  of the SphericalSharp Green's function against analytical result
 */
BOOST_FIXTURE_TEST_CASE(hessian_AD, SphericalSharpTest)
{
    SphericalSharp<AD_hessian, CollocationIntegrator> gf(eps, epsilonSolvent, radius, origin);
    double analytic = result(0);
    double gf_value = gf.kernelS(source, probe);
    BOOST_REQUIRE_CLOSE(analytic, gf_value, 1.0e-12);

    double analytic_derProbe = result(1);
    double gf_derProbe = gf.derivativeProbe(probeNormal, source, probe);
    BOOST_REQUIRE_CLOSE(analytic_derProbe, gf_derProbe, 1.0e-12);

    double analytic_derSource = result(2);
    double gf_derSource = gf.derivativeSource(sourceNormal, source, probe);
    BOOST_REQUIRE_CLOSE(analytic_derSource, gf_derSource, 1.0e-12);
}

BOOST_AUTO_TEST_SUITE_END()
