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

#define BOOST_TEST_MODULE GreensFunctionMetalNP

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include <cmath>
#include <iostream>

#include "Config.hpp"

#include <Eigen/Dense>

#include "CollocationIntegrator.hpp"
#include "AnalyticEvaluate.hpp"
#include "MetalNP.hpp"
#include "DerivativeTypes.hpp"

struct MetalNPTest {
    double epsilonSolvent;
    double epsRe;
    double epsIm;
    double radius;
    Eigen::Vector3d origin;
    Eigen::Vector3d source, probe, sourceNormal, probeNormal;
    Eigen::Array4d result;
    MetalNPTest() { SetUp(); }
    void SetUp() {
        epsilonSolvent = 80.0;
        epsRe = -9.906;
        epsIm = -4.329;
        radius = 50.0;
        origin = Eigen::Vector3d::Zero();
        source = Eigen::Vector3d::Random();
        sourceNormal = source + Eigen::Vector3d::Random();
        sourceNormal.normalize();
        probe = Eigen::Vector3d::Random();
        probeNormal = probe + Eigen::Vector3d::Random();
        probeNormal.normalize();
        result = analyticUniformDielectric(epsilonSolvent, sourceNormal, source, probeNormal, probe);
        //result = analyticMetalNP(epsilonSolvent, espRe, epsIm, sourceNormal, source, probeNormal, probe);
    }
};

BOOST_FIXTURE_TEST_SUITE(Nanoparticle, MetalNPTest)

/*! \class MetalNP
 *  \test \b MetalNPTest_numerical tests the numerical evaluation of the MetalNP Green's function against analytical result
 */
BOOST_FIXTURE_TEST_CASE(numerical, MetalNPTest)
{
    MetalNP<Numerical, CollocationIntegrator> gf(epsilonSolvent, epsRe, epsIm, origin, radius);
    double value = result(0);
    double gf_value = gf.kernelS(source, probe);
    BOOST_REQUIRE_CLOSE(value, gf_value, 1.0e-12);

    double derProbe = result(1);
    double gf_derProbe = gf.derivativeProbe(probeNormal, source, probe);
    BOOST_REQUIRE_CLOSE(derProbe, gf_derProbe, 1.0e-05);

    double derSource = result(2);
    double gf_derSource = gf.derivativeSource(sourceNormal, source, probe);
    BOOST_REQUIRE_CLOSE(derSource, gf_derSource, 1.0e-05);
}

/*! \class MetalNP
 *  \test \b MetalNPTest_directional_AD tests the automatic evaluation (directional derivative only)
 *  of the MetalNP Green's function against analytical result
 */
BOOST_FIXTURE_TEST_CASE(directional_AD, MetalNPTest)
{
    MetalNP<AD_directional, CollocationIntegrator> gf(epsilonSolvent, epsRe, epsIm, origin, radius);
    double value = result(0);
    double gf_value = gf.kernelS(source, probe);
    BOOST_REQUIRE_CLOSE(value, gf_value, 1.0e-12);

    double derProbe = result(1);
    double gf_derProbe = gf.derivativeProbe(probeNormal, source, probe);
    BOOST_REQUIRE_CLOSE(derProbe, gf_derProbe, 1.0e-12);

    double derSource = result(2);
    double gf_derSource = gf.derivativeSource(sourceNormal, source, probe);
    BOOST_REQUIRE_CLOSE(derSource, gf_derSource, 1.0e-12);
}

/*! \class MetalNP
 *  \test \b MetalNPTest_gradient_AD tests the automatic evaluation (full gradient)
 *  of the MetalNP Green's function against analytical result
 */
BOOST_FIXTURE_TEST_CASE(gradient_AD, MetalNPTest)
{
    MetalNP<AD_gradient, CollocationIntegrator> gf(epsilonSolvent, epsRe, epsIm, origin, radius);
    double value = result(0);
    double gf_value = gf.kernelS(source, probe);
    BOOST_REQUIRE_CLOSE(value, gf_value, 1.0e-12);

    double derProbe = result(1);
    double gf_derProbe = gf.derivativeProbe(probeNormal, source, probe);
    BOOST_REQUIRE_CLOSE(derProbe, gf_derProbe, 1.0e-12);

    double derSource = result(2);
    double gf_derSource = gf.derivativeSource(sourceNormal, source, probe);
    BOOST_REQUIRE_CLOSE(derSource, gf_derSource, 1.0e-12);
}

/*! \class MetalNP
 *  \test \b MetalNPTest_hessian_AD tests the automatic evaluation (full hessian)
 *  of the MetalNP Green's function against analytical result
 */
BOOST_FIXTURE_TEST_CASE(hessian_AD, MetalNPTest)
{
    MetalNP<AD_hessian, CollocationIntegrator> gf(epsilonSolvent, epsRe, epsIm, origin, radius);
    double value = result(0);
    double gf_value = gf.kernelS(source, probe);
    BOOST_REQUIRE_CLOSE(value, gf_value, 1.0e-12);

    double derProbe = result(1);
    double gf_derProbe = gf.derivativeProbe(probeNormal, source, probe);
    BOOST_REQUIRE_CLOSE(derProbe, gf_derProbe, 1.0e-12);

    double derSource = result(2);
    double gf_derSource = gf.derivativeSource(sourceNormal, source, probe);
    BOOST_REQUIRE_CLOSE(derSource, gf_derSource, 1.0e-12);
}

BOOST_AUTO_TEST_SUITE_END()
