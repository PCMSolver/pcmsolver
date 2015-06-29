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

#include <Eigen/Core>

#include "cnpyPimpl.hpp"
#include "DerivativeTypes.hpp"
#include "IntegratorTypes.hpp"
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
 *  \test \b CollocationIntegratorTest_vacuum tests the evaluation by collocation of the vacuum matrix representations of S and D
 */
BOOST_FIXTURE_TEST_CASE(vacuum, CollocationIntegratorTest)
{
    Vacuum<AD_directional, CollocationIntegrator> gf;

    BOOST_TEST_MESSAGE("Vacuum");
    Eigen::MatrixXd S_results = gf.singleLayer(cavity.elements());
    cnpy::NpyArray raw_S_ref = cnpy::npy_load("vacuum_S_collocation.npy");
    int dim_read = raw_S_ref.shape[0];
    Eigen::MatrixXd S_reference = Eigen::MatrixXd::Zero(dim_read, dim_read);
    S_reference = getFromRawBuffer<double>(dim_read, dim_read, raw_S_ref.data);
    for (int i = 0; i < cavity.size(); ++i) {
        for (int j = 0; j < cavity.size(); ++j) {
            BOOST_REQUIRE_CLOSE(S_reference(i, j), S_results(i, j), 1.0e-12);
        }
    }

    Eigen::MatrixXd D_results = gf.doubleLayer(cavity.elements());
    cnpy::NpyArray raw_D_ref = cnpy::npy_load("vacuum_D_collocation.npy");
    dim_read = raw_D_ref.shape[0];
    Eigen::MatrixXd D_reference = Eigen::MatrixXd::Zero(dim_read, dim_read);
    D_reference = getFromRawBuffer<double>(dim_read, dim_read, raw_D_ref.data);
    for (int i = 0; i < cavity.size(); ++i) {
        for (int j = 0; j < cavity.size(); ++j) {
            BOOST_REQUIRE_CLOSE(D_reference(i, j), D_results(i, j), 1.0e-12);
        }
    }
}

/*! \class CollocationIntegrator
 *  \test \b CollocationIntegratorTest_uniformdielectric tests the evaluation by collocation of the uniform dielectric matrix representations of S and D
 */
BOOST_FIXTURE_TEST_CASE(uniformdielectric, CollocationIntegratorTest)
{
    UniformDielectric<AD_directional, CollocationIntegrator> gf(epsilon);

    BOOST_TEST_MESSAGE("UniformDielectric");
    Eigen::MatrixXd S_results = gf.singleLayer(cavity.elements());
    cnpy::NpyArray raw_S_ref = cnpy::npy_load("uniformdielectric_S_collocation.npy");
    int dim_read = raw_S_ref.shape[0];
    Eigen::MatrixXd S_reference = Eigen::MatrixXd::Zero(dim_read, dim_read);
    S_reference = getFromRawBuffer<double>(dim_read, dim_read, raw_S_ref.data);
    for (int i = 0; i < cavity.size(); ++i) {
        for (int j = 0; j < cavity.size(); ++j) {
            BOOST_REQUIRE_CLOSE(S_reference(i, j), S_results(i, j), 1.0e-12);
        }
    }

    Eigen::MatrixXd D_results = gf.doubleLayer(cavity.elements());
    cnpy::NpyArray raw_D_ref = cnpy::npy_load("uniformdielectric_D_collocation.npy");
    dim_read = raw_D_ref.shape[0];
    Eigen::MatrixXd D_reference = Eigen::MatrixXd::Zero(dim_read, dim_read);
    D_reference = getFromRawBuffer<double>(dim_read, dim_read, raw_D_ref.data);
    for (int i = 0; i < cavity.size(); ++i) {
        for (int j = 0; j < cavity.size(); ++j) {
            BOOST_REQUIRE_CLOSE(D_reference(i, j), D_results(i, j), 1.0e-12);
        }
    }
}

/*! \class CollocationIntegrator
 *  \test \b CollocationIntegratorTest_tanhsphericaldiffuse tests the evaluation by collocation of the spherical diffuse matrix representations of S and D
 */
BOOST_FIXTURE_TEST_CASE(tanhsphericaldiffuse, CollocationIntegratorTest)
{
    double width = 5.0;
    double sphereRadius = 100.0;
    SphericalDiffuse<CollocationIntegrator, OneLayerTanh> gf(epsilon, epsilon, width, sphereRadius, Eigen::Vector3d::Zero(), 3);

    BOOST_TEST_MESSAGE("TanhSphericalDiffuse");
    Eigen::MatrixXd S_results = gf.singleLayer(cavity.elements());
    cnpy::NpyArray raw_S_ref = cnpy::npy_load("tanhsphericaldiffuse_S_collocation.npy");
    int dim_read = raw_S_ref.shape[0];
    Eigen::MatrixXd S_reference = Eigen::MatrixXd::Zero(dim_read, dim_read);
    S_reference = getFromRawBuffer<double>(dim_read, dim_read, raw_S_ref.data);
    for (int i = 0; i < cavity.size(); ++i) {
        for (int j = 0; j < cavity.size(); ++j) {
            BOOST_REQUIRE_CLOSE(S_reference(i, j), S_results(i, j), 1.0e-12);
        }
    }

    Eigen::MatrixXd D_results = gf.doubleLayer(cavity.elements());
    cnpy::NpyArray raw_D_ref = cnpy::npy_load("tanhsphericaldiffuse_D_collocation.npy");
    dim_read = raw_D_ref.shape[0];
    Eigen::MatrixXd D_reference = Eigen::MatrixXd::Zero(dim_read, dim_read);
    D_reference = getFromRawBuffer<double>(dim_read, dim_read, raw_D_ref.data);
    for (int i = 0; i < cavity.size(); ++i) {
        for (int j = 0; j < cavity.size(); ++j) {
            BOOST_REQUIRE_CLOSE(D_reference(i, j), D_results(i, j), 1.0e-12);
        }
    }
}
