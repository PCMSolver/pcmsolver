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

#define BOOST_TEST_MODULE GreensFunctionIonicLiquid

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include <cmath>
#include <iostream>

#include "Config.hpp"

#include <Eigen/Dense>

#include "AnalyticEvaluate.hpp"
#include "DerivativeTypes.hpp"
#include "IonicLiquid.hpp"
#include "CollocationIntegrator.hpp"

struct IonicLiquidTest {
    double epsilon;
    double kappa;
    Eigen::Vector3d source, probe, sourceNormal, probeNormal;
    Eigen::Array4d result;
    IonicLiquidTest() { SetUp(); }
    void SetUp() {
        epsilon = 60.0;
        kappa = 5.0;
        source = Eigen::Vector3d::Random();
        sourceNormal = source + Eigen::Vector3d::Random();
        sourceNormal.normalize();
        probe = Eigen::Vector3d::Random();
        probeNormal = probe + Eigen::Vector3d::Random();
        probeNormal.normalize();
        result = analyticIonicLiquid(epsilon, kappa, sourceNormal, source, probeNormal, probe);
    }
};

/*! \class IonicLiquid
 *  \test \b IonicLiquidTest_numerical tests the numerical evaluation of the IonicLiquid Green's function against analytical result
 */
BOOST_FIXTURE_TEST_CASE(numerical, IonicLiquidTest)
{
    IonicLiquid<Numerical, CollocationIntegrator<Numerical, Yukawa> > gf(epsilon, kappa);
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

/*! \class IonicLiquid
 *  \test \b IonicLiquidTest_directional_AD tests the automatic evaluation (directional derivative only)
 *  of the IonicLiquid Green's function against analytical result
 */
BOOST_FIXTURE_TEST_CASE(directional_AD, IonicLiquidTest)
{
    IonicLiquid<AD_directional, CollocationIntegrator<AD_directional, Yukawa> > gf(epsilon, kappa);
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

/*! \class IonicLiquid
 *  \test \b IonicLiquidTest_gradient_AD tests the automatic evaluation (full gradient)
 *  of the IonicLiquid Green's function against analytical result
 */
BOOST_FIXTURE_TEST_CASE(gradient_AD, IonicLiquidTest)
{
    IonicLiquid<AD_gradient, CollocationIntegrator<AD_gradient, Yukawa> > gf(epsilon, kappa);
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

/*! \class IonicLiquid
 *  \test \b IonicLiquidTest_hessian_AD tests the automatic evaluation (full hessian)
 *  of the IonicLiquid Green's function against analytical result
 */
BOOST_FIXTURE_TEST_CASE(hessian_AD, IonicLiquidTest)
{
    IonicLiquid<AD_hessian, CollocationIntegrator<AD_hessian, Yukawa> > gf(epsilon, kappa);
    double value = result(0);
    double gf_value = gf.kernelS(source, probe);
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
