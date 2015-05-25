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

#define BOOST_TEST_MODULE CollocationIntegrator

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include <cmath>
#include <iostream>
#include <iomanip>
#include <limits>

#include "Config.hpp"

#include <Eigen/Dense>

#include "cnpyPimpl.hpp"
#include "DerivativeTypes.hpp"
#include "Element.hpp"
#include "GePolCavity.hpp"
#include "CollocationIntegrator.hpp"
#include "PhysicalConstants.hpp"
#include "SphericalDiffuse.hpp"
#include "TestingMolecules.hpp"
#include "UniformDielectric.hpp"
#include "Vacuum.hpp"

struct CollocationIntegratorTest {
    double radius;
    double epsilon;
    Eigen::Vector3d source, probe, sourceNormal, probeNormal;
    GePolCavity cavity;
    CollocationIntegratorTest() { SetUp(); }
    void SetUp() {
	    epsilon = 80.0;
        source = Eigen::Vector3d::Random();
        sourceNormal = source + Eigen::Vector3d::Random();
        sourceNormal.normalize();
        probe = Eigen::Vector3d::Random();
        probeNormal = probe + Eigen::Vector3d::Random();
        probeNormal.normalize();

	    Molecule molec = dummy<0>(1.44 / convertBohrToAngstrom);
        double area = 10.0;
        cavity = GePolCavity(molec, area, 0.0, 100.0);
    }
};

/*! \class CollocationIntegrator
 *  \test \b CollocationIntegratorTest_vacuum tests the numerical evaluation of the vacuum diagonal elements of S and D
 */
BOOST_FIXTURE_TEST_CASE(vacuum, CollocationIntegratorTest)
{
    Vacuum<AD_directional, CollocationIntegrator<AD_directional, Uniform> > gf;

    BOOST_TEST_MESSAGE("Vacuum");
    Eigen::VectorXd S_results = Eigen::VectorXd::Zero(cavity.size());
    for (int i = 0; i < cavity.size(); ++i) {
	S_results(i) = gf.diagonalS(cavity.elements(i));
    }
    // In case you need to update the reference files...
    /*
    unsigned int dim = static_cast<unsigned int>(cavity.size());
    const unsigned int shape[] = {dim};
    cnpy::npy_save("vacuum_S_collocation.npy", S_results.data(), shape, 1, "w", false);
    */
    cnpy::NpyArray raw_S_ref = cnpy::npy_load("vacuum_S_collocation.npy");
    int dim_read = raw_S_ref.shape[0];
    Eigen::VectorXd S_reference = Eigen::VectorXd::Zero(dim_read);
    S_reference = getFromRawBuffer<double>(dim_read, 1, raw_S_ref.data);
    for (int i = 0; i < cavity.size(); ++i) {
    	BOOST_REQUIRE_CLOSE(S_reference(i), S_results(i), 1.0e-12);
    }
    BOOST_TEST_MESSAGE("S operator diagonal by collocation");
    for (int i = 0; i < cavity.size(); ++i) {
	    BOOST_TEST_MESSAGE("S_{" << i+1 << ", " << i+1 << "} = "
			    << std::setprecision(std::numeric_limits<long double>::digits10) << S_results(i));
    }

    Eigen::VectorXd D_results = Eigen::VectorXd::Zero(cavity.size());
    for (int i = 0; i < cavity.size(); ++i) {
	D_results(i) = gf.diagonalD(cavity.elements(i));
    }
    // In case you need to update the reference files...
    /*
    cnpy::npy_save("vacuum_D_collocation.npy", D_results.data(), shape, 1, "w", false);
    */
    cnpy::NpyArray raw_D_ref = cnpy::npy_load("vacuum_D_collocation.npy");
    dim_read = raw_D_ref.shape[0];
    Eigen::VectorXd D_reference = Eigen::VectorXd::Zero(dim_read);
    D_reference = getFromRawBuffer<double>(dim_read, 1, raw_D_ref.data);
    for (int i = 0; i < cavity.size(); ++i) {
    	BOOST_REQUIRE_CLOSE(D_reference(i), D_results(i), 1.0e-12);
    }
    BOOST_TEST_MESSAGE("D operator diagonal by collocation");
    for (int i = 0; i < cavity.size(); ++i) {
	    BOOST_TEST_MESSAGE("D_{" << i+1 << ", " << i+1 << "} = "
			    << std::setprecision(std::numeric_limits<long double>::digits10) << D_results(i));
    }
}

/*! \class CollocationIntegrator
 *  \test \b CollocationIntegratorTest_uniformdielectric tests the numerical evaluation of the uniform dielectric diagonal elements of S and D
 */
BOOST_FIXTURE_TEST_CASE(uniformdielectric, CollocationIntegratorTest)
{
    UniformDielectric<AD_directional, CollocationIntegrator<AD_directional, Uniform> > gf(epsilon);

    BOOST_TEST_MESSAGE("UniformDielectric");
    Eigen::VectorXd S_results = Eigen::VectorXd::Zero(cavity.size());
    for (int i = 0; i < cavity.size(); ++i) {
	S_results(i) = gf.diagonalS(cavity.elements(i));
    }
    // In case you need to update the reference files...
    /*
    unsigned int dim = static_cast<unsigned int>(cavity.size());
    const unsigned int shape[] = {dim};
    cnpy::npy_save("uniformdielectric_S_collocation.npy", S_results.data(), shape, 1, "w", false);
    */
    cnpy::NpyArray raw_S_ref = cnpy::npy_load("uniformdielectric_S_collocation.npy");
    int dim_read = raw_S_ref.shape[0];
    Eigen::VectorXd S_reference = Eigen::VectorXd::Zero(dim_read);
    S_reference = getFromRawBuffer<double>(dim_read, 1, raw_S_ref.data);
    for (int i = 0; i < cavity.size(); ++i) {
    	BOOST_REQUIRE_CLOSE(S_reference(i), S_results(i), 1.0e-12);
    }
    BOOST_TEST_MESSAGE("S operator diagonal by collocation");
    for (int i = 0; i < cavity.size(); ++i) {
	    BOOST_TEST_MESSAGE("S_{" << i+1 << ", " << i+1 << "} = "
			    << std::setprecision(std::numeric_limits<long double>::digits10) << S_results(i));
    }

    Eigen::VectorXd D_results = Eigen::VectorXd::Zero(cavity.size());
    for (int i = 0; i < cavity.size(); ++i) {
	D_results(i) = gf.diagonalD(cavity.elements(i));
    }
    // In case you need to update the reference files...
    /*
    cnpy::npy_save("uniformdielectric_D_collocation.npy", D_results.data(), shape, 1, "w", false);
    */
    cnpy::NpyArray raw_D_ref = cnpy::npy_load("uniformdielectric_D_collocation.npy");
    dim_read = raw_D_ref.shape[0];
    Eigen::VectorXd D_reference = Eigen::VectorXd::Zero(dim_read);
    D_reference = getFromRawBuffer<double>(dim_read, 1, raw_D_ref.data);
    for (int i = 0; i < cavity.size(); ++i) {
    	BOOST_REQUIRE_CLOSE(D_reference(i), D_results(i), 1.0e-12);
    }
    BOOST_TEST_MESSAGE("D operator diagonal by collocation");
    for (int i = 0; i < cavity.size(); ++i) {
	    BOOST_TEST_MESSAGE("D_{" << i+1 << ", " << i+1 << "} = "
			    << std::setprecision(std::numeric_limits<long double>::digits10) << D_results(i));
    }
}

BOOST_FIXTURE_TEST_CASE(tanhsphericaldiffuse, CollocationIntegratorTest)
{
    double width = 5.0;
    double sphereRadius = 100.0;
    SphericalDiffuse<CollocationIntegrator<Numerical, TanhDiffuse>, TanhDiffuse> gf(epsilon, epsilon, width, sphereRadius, Eigen::Vector3d::Zero());

    BOOST_TEST_MESSAGE("TanhSphericalDiffuse");
    Eigen::VectorXd S_results = Eigen::VectorXd::Zero(cavity.size());
    for (int i = 0; i < cavity.size(); ++i) {
	    S_results(i) = gf.diagonalS(cavity.elements(i));
    }
    // In case you need to update the reference files...
    /*
    unsigned int dim = static_cast<unsigned int>(cavity.size());
    const unsigned int shape[] = {dim};
    cnpy::npy_save("tanhsphericaldiffuse_S_collocation.npy", S_results.data(), shape, 1, "w", false);
    */
    cnpy::NpyArray raw_S_ref = cnpy::npy_load("tanhsphericaldiffuse_S_collocation.npy");
    int dim_read = raw_S_ref.shape[0];
    Eigen::VectorXd S_reference = Eigen::VectorXd::Zero(dim_read);
    S_reference = getFromRawBuffer<double>(dim_read, 1, raw_S_ref.data);
    for (int i = 0; i < cavity.size(); ++i) {
    	BOOST_CHECK_CLOSE(S_reference(i), S_results(i), 1.0e-08);
    }
    BOOST_TEST_MESSAGE("S operator diagonal by collocation");
    for (int i = 0; i < cavity.size(); ++i) {
	    BOOST_TEST_MESSAGE("S_{" << i+1 << ", " << i+1 << "} = "
			    << std::setprecision(std::numeric_limits<long double>::digits10) << S_results(i));
    }

    Eigen::VectorXd D_results = Eigen::VectorXd::Zero(cavity.size());
    for (int i = 0; i < cavity.size(); ++i) {
	    D_results(i) = gf.diagonalD(cavity.elements(i));
    }
    // In case you need to update the reference files...
    /*
    cnpy::npy_save("tanhsphericaldiffuse_D_collocation.npy", D_results.data(), shape, 1, "w", false);
    */
    cnpy::NpyArray raw_D_ref = cnpy::npy_load("tanhsphericaldiffuse_D_collocation.npy");
    dim_read = raw_D_ref.shape[0];
    Eigen::VectorXd D_reference = Eigen::VectorXd::Zero(dim_read);
    D_reference = getFromRawBuffer<double>(dim_read, 1, raw_D_ref.data);
    for (int i = 0; i < cavity.size(); ++i) {
    	BOOST_CHECK_CLOSE(D_reference(i), D_results(i), 1.0e-06);
    }
    BOOST_TEST_MESSAGE("D operator diagonal by collocation");
    for (int i = 0; i < cavity.size(); ++i) {
	    BOOST_TEST_MESSAGE("D_{" << i+1 << ", " << i+1 << "} = "
			    << std::setprecision(std::numeric_limits<long double>::digits10) << D_results(i));
    }
}
