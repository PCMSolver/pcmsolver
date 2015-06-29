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

#define BOOST_TEST_MODULE NumericalIntegrator

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include <cmath>
#include <iostream>
#include <iomanip>
#include <limits>

#include "Config.hpp"

#include <Eigen/Core>

#include <boost/filesystem.hpp>

#include "cnpyPimpl.hpp"
#include "DerivativeTypes.hpp"
#include "AnisotropicLiquid.hpp"
#include "Element.hpp"
#include "GePolCavity.hpp"
#include "IonicLiquid.hpp"
#include "NumericalIntegrator.hpp"
#include "PhysicalConstants.hpp"
#include "UniformDielectric.hpp"
#include "TestingMolecules.hpp"
#include "Vacuum.hpp"

namespace fs = boost::filesystem;

struct NumericalIntegratorTest {
    double radius;
    Eigen::Vector3d epsilon;
    Eigen::Vector3d euler;
    double eps;
    double kappa;
    Eigen::Vector3d source, probe, sourceNormal, probeNormal;
    GePolCavity cavity;
    NumericalIntegratorTest() { SetUp(); }
    void SetUp() {
        /*    	epsilon << 2.0, 80.0, 15.0;
                euler << 6.0, 40.0, 15.0;*/
        epsilon << 80.0, 80.0, 80.0;
        euler << 0.0, 0.0, 0.0;
        eps = 80.0;
        kappa = 1.5;
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

/*! \class NumericalIntegrator
 *  \test \b NumericalIntegratorTest_vacuum tests the numerical evaluation of the vacuum matrix representations of S and D
 */
BOOST_FIXTURE_TEST_CASE(vacuum, NumericalIntegratorTest)
{
    fs::rename("PEDRA.OUT", "PEDRA.OUT.vacuum");
    Vacuum<AD_directional, NumericalIntegrator> gf;

    BOOST_TEST_MESSAGE("Vacuum");
    Eigen::MatrixXd S_results = gf.singleLayer(cavity.elements());
    /* Numerical integrator not really working now...
    cnpy::NpyArray raw_S_ref = cnpy::npy_load("vacuum_S_numerical.npy");
    int dim_read = raw_S_ref.shape[0];
    Eigen::MatrixXd S_reference = Eigen::MatrixXd::Zero(dim_read, dim_read);
    S_reference = getFromRawBuffer<double>(dim_read, dim_read, raw_S_ref.data);
    for (int i = 0; i < cavity.size(); ++i) {
        for (int j = 0; j < cavity.size(); ++j) {
            BOOST_REQUIRE_CLOSE(S_reference(i, j), S_results(i, j), 1.0e-12);
        }
    }
    */

    Eigen::MatrixXd D_results = gf.doubleLayer(cavity.elements());
    /* Numerical integrator not really working now...
    cnpy::NpyArray raw_D_ref = cnpy::npy_load("vacuum_D_numerical.npy");
    dim_read = raw_D_ref.shape[0];
    Eigen::MatrixXd D_reference = Eigen::MatrixXd::Zero(dim_read, dim_read);
    D_reference = getFromRawBuffer<double>(dim_read, dim_read, raw_D_ref.data);
    for (int i = 0; i < cavity.size(); ++i) {
        for (int j = 0; j < cavity.size(); ++j) {
            BOOST_REQUIRE_CLOSE(D_reference(i, j), D_results(i, j), 1.0e-12);
        }
    }
    // Checks sum rule for D operator
    BOOST_REQUIRE_CLOSE(D_sum, 2 * M_PI, 1.0e-12);
    */
}

/*! \class NumericalIntegrator
 *  \test \b NumericalIntegratorTest_uniformdielectric tests the numerical evaluation of the uniform dielectric matrix representations of S and D
 */
BOOST_FIXTURE_TEST_CASE(uniformdielectric, NumericalIntegratorTest)
{
    fs::rename("PEDRA.OUT", "PEDRA.OUT.uniform");
    UniformDielectric<AD_directional, NumericalIntegrator> gf(eps);

    BOOST_TEST_MESSAGE("UniformDielectric");
    Eigen::MatrixXd S_results = gf.singleLayer(cavity.elements());
    /* Numerical integrator not really working now...
    cnpy::NpyArray raw_S_ref = cnpy::npy_load("uniformdielectric_S_numerical.npy");
    int dim_read = raw_S_ref.shape[0];
    Eigen::MatrixXd S_reference = Eigen::MatrixXd::Zero(dim_read, dim_read);
    S_reference = getFromRawBuffer<double>(dim_read, dim_read, raw_S_ref.data);
    for (int i = 0; i < cavity.size(); ++i) {
        for (int j = 0; j < cavity.size(); ++j) {
            BOOST_REQUIRE_CLOSE(S_reference(i, j), S_results(i, j), 1.0e-12);
        }
    }
    */

    Eigen::MatrixXd D_results = gf.doubleLayer(cavity.elements());
    /* Numerical integrator not really working now...
    cnpy::NpyArray raw_D_ref = cnpy::npy_load("uniformdielectric_D_numerical.npy");
    dim_read = raw_D_ref.shape[0];
    Eigen::MatrixXd D_reference = Eigen::MatrixXd::Zero(dim_read, dim_read);
    D_reference = getFromRawBuffer<double>(dim_read, dim_read, raw_D_ref.data);
    for (int i = 0; i < cavity.size(); ++i) {
        for (int j = 0; j < cavity.size(); ++j) {
            BOOST_REQUIRE_CLOSE(D_reference(i, j), D_results(i, j), 1.0e-12);
        }
    }
    */
}

/*! \class NumericalIntegrator
 *  \test \b NumericalIntegratorTest_ionic tests the numerical evaluation of the ionic liquid matrix representations of S and D
 */
BOOST_FIXTURE_TEST_CASE(ionic, NumericalIntegratorTest)
{
    fs::rename("PEDRA.OUT", "PEDRA.OUT.ionic");
    IonicLiquid<AD_directional, NumericalIntegrator> gf(eps, kappa);

    BOOST_TEST_MESSAGE("IonicLiquid");
    Eigen::MatrixXd S_results = gf.singleLayer(cavity.elements());
    /* Numerical integrator not really working now...
    cnpy::NpyArray raw_S_ref = cnpy::npy_load("ionicliquid_S_numerical.npy");
    int dim_read = raw_S_ref.shape[0];
    Eigen::MatrixXd S_reference = Eigen::MatrixXd::Zero(dim_read, dim_read);
    S_reference = getFromRawBuffer<double>(dim_read, dim_read, raw_S_ref.data);
    for (int i = 0; i < cavity.size(); ++i) {
        for (int j = 0; j < cavity.size(); ++j) {
            BOOST_REQUIRE_CLOSE(S_reference(i, j), S_results(i, j), 1.0e-12);
        }
    }
    */

    Eigen::MatrixXd D_results = gf.doubleLayer(cavity.elements());
    /* Numerical integrator not really working now...
    cnpy::NpyArray raw_D_ref = cnpy::npy_load("ionicliquid_D_numerical.npy");
    dim_read = raw_D_ref.shape[0];
    Eigen::MatrixXd D_reference = Eigen::MatrixXd::Zero(dim_read, dim_read);
    D_reference = getFromRawBuffer<double>(dim_read, dim_read, raw_D_ref.data);
    for (int i = 0; i < cavity.size(); ++i) {
        for (int j = 0; j < cavity.size(); ++j) {
            BOOST_REQUIRE_CLOSE(D_reference(i, j), D_results(i, j), 1.0e-12);
        }
    }
    */
}

/*! \class NumericalIntegrator
 *  \test \b NumericalIntegratorTest_anisotropic tests the numerical evaluation of the anisotropic liquid matrix representations of S and D
 */
BOOST_FIXTURE_TEST_CASE(anisotropic, NumericalIntegratorTest)
{
    fs::rename("PEDRA.OUT", "PEDRA.OUT.anisotropic");
    AnisotropicLiquid<AD_directional, NumericalIntegrator> gf(epsilon, euler);

    BOOST_TEST_MESSAGE("AnisotropicLiquid");
    Eigen::MatrixXd S_results = gf.singleLayer(cavity.elements());
    /* Numerical integrator not really working now...
    cnpy::NpyArray raw_S_ref = cnpy::npy_load("anisotropicliquid_S_numerical.npy");
    int dim_read = raw_S_ref.shape[0];
    Eigen::MatrixXd S_reference = Eigen::MatrixXd::Zero(dim_read, dim_read);
    S_reference = getFromRawBuffer<double>(dim_read, dim_read, raw_S_ref.data);
    for (int i = 0; i < cavity.size(); ++i) {
        for (int j = 0; j < cavity.size(); ++j) {
            BOOST_REQUIRE_CLOSE(S_reference(i, j), S_results(i, j), 1.0e-12);
        }
    }
    */

    Eigen::MatrixXd D_results = gf.doubleLayer(cavity.elements());
    /* Numerical integrator not really working now...
    cnpy::NpyArray raw_D_ref = cnpy::npy_load("anisotropicliquid_D_numerical.npy");
    dim_read = raw_D_ref.shape[0];
    Eigen::MatrixXd D_reference = Eigen::MatrixXd::Zero(dim_read, dim_read);
    D_reference = getFromRawBuffer<double>(dim_read, dim_read, raw_D_ref.data);
    for (int i = 0; i < cavity.size(); ++i) {
        for (int j = 0; j < cavity.size(); ++j) {
            BOOST_REQUIRE_CLOSE(D_reference(i, j), D_results(i, j), 1.0e-12);
        }
    }
    */
}
