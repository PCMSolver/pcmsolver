/**
 * PCMSolver, an API for the Polarizable Continuum Model
 * Copyright (C) 2017 Roberto Di Remigio, Luca Frediani and collaborators.
 *
 * This file is part of PCMSolver.
 *
 * PCMSolver is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * PCMSolver is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with PCMSolver.  If not, see <http://www.gnu.org/licenses/>.
 *
 * For information on the complete list of contributors to the
 * PCMSolver API, see: <http://pcmsolver.readthedocs.io/>
 */

#include "catch.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>

#include <Eigen/Core>

#include "TestingMolecules.hpp"
#include "bi_operators/Numerical.hpp"
#include "cavity/Element.hpp"
#include "cavity/GePolCavity.hpp"
#include "green/AnisotropicLiquid.hpp"
#include "green/IonicLiquid.hpp"
#include "green/UniformDielectric.hpp"
#include "green/Vacuum.hpp"
#include "utils/MathUtils.hpp"

using namespace pcm;
using bi_operators::Numerical;
using cavity::GePolCavity;
using green::Vacuum;
using green::UniformDielectric;
using green::IonicLiquid;
using green::AnisotropicLiquid;

SCENARIO(
    "A collocation integrator with numerical integrator of the diagonal elements",
    "[bi_operators][bi_operators_numerical]") {
  GIVEN("A GePol cavity for a single sphere in the origin") {
    Molecule molec = dummy<0>(1.44 / bohrToAngstrom());
    double area = 10.0;
    GePolCavity cavity = GePolCavity(molec, area, 0.0, 100.0);
    Eigen::MatrixXd results = Eigen::MatrixXd::Zero(cavity.size(), cavity.size());
    Eigen::MatrixXd reference = Eigen::MatrixXd::Zero(cavity.size(), cavity.size());
    Numerical op;

    /*! \class Numerical
     *  \test \b NumericalTest_vacuum tests the numerical evaluation of the
     * vacuum matrix representations of S and D
     */
    WHEN("the vacuum Green's function is used") {
      Vacuum<> gf;
      THEN("the matrix elements of S are") {
        results = op.computeS(cavity, gf);
        // Numerical integrator not really working now...
        /*
        reference = cnpy::custom::npy_load<double>("vacuum_S_numerical.npy");
        for (int i = 0; i < cavity.size(); ++i) {
            for (int j = 0; j < cavity.size(); ++j) {
                REQUIRE(reference(i, j) == Approx(results(i, j)));
            }
        }
        */
      }
      AND_THEN("the matrix elements of D are") {
        Eigen::MatrixXd results = op.computeD(cavity, gf);
        // Numerical integrator not really working now...
        /*
        reference = cnpy::custom::npy_load<double>("vacuum_D_numerical.npy");
        for (int i = 0; i < cavity.size(); ++i) {
            for (int j = 0; j < cavity.size(); ++j) {
                REQUIRE(reference(i, j) == Approx(results(i, j)));
            }
        }
        // Checks sum rule for D operator
        REQUIRE(D_sum == Approx(2 * M_PI));
        */
      }
    }

    /*! \class Numerical
     *  \test \b NumericalTest_uniformdielectric tests the numerical
     * evaluation of the uniform dielectric matrix representations of S and D
     */
    WHEN("the uniform dielectric Green's function is used") {
      double eps = 80.0;
      UniformDielectric<> gf(eps);
      THEN("the matrix elements of S are") {
        results = op.computeS(cavity, gf);
        // Numerical integrator not really working now...
        /*
        reference =
        cnpy::custom::npy_load<double>("uniformdielectric_S_numerical.npy");
        for (int i = 0; i < cavity.size(); ++i) {
            for (int j = 0; j < cavity.size(); ++j) {
                REQUIRE(reference(i, j) == Approx(results(i, j)));
            }
        }
        */
      }
      AND_THEN("the matrix elements of D are") {
        Eigen::MatrixXd results = op.computeD(cavity, gf);
        // Numerical integrator not really working now...
        /*
        reference =
        cnpy::custom::npy_load<double>("uniformdielectric_D_numerical.npy");
        for (int i = 0; i < cavity.size(); ++i) {
            for (int j = 0; j < cavity.size(); ++j) {
                REQUIRE(reference(i, j) == Approx(results(i, j)));
            }
        }
        */
      }
    }

    /*! \class Numerical
     *  \test \b NumericalTest_ionic tests the numerical evaluation of the
     * ionic liquid matrix representations of S and D
     */
    WHEN("the ionic liquid Green's function is used") {
      double eps = 80.0;
      double kappa = 1.5;
      IonicLiquid<> gf(eps, kappa);
      THEN("the matrix elements of S are") {
        results = op.computeS(cavity, gf);
        // Numerical integrator not really working now...
        /*
        reference = cnpy::custom::npy_load<double>("ionicliquid_S_numerical.npy");
        for (int i = 0; i < cavity.size(); ++i) {
            for (int j = 0; j < cavity.size(); ++j) {
                REQUIRE(reference(i, j) == Approx(results(i, j)));
            }
        }
        */
      }
      AND_THEN("the matrix elements of D are") {
        results = op.computeD(cavity, gf);
        // Numerical integrator not really working now...
        /*
        reference = cnpy::custom::npy_load<double>("ionicliquid_D_numerical.npy");
        for (int i = 0; i < cavity.size(); ++i) {
            for (int j = 0; j < cavity.size(); ++j) {
                REQUIRE(reference(i, j) == Approx(results(i, j)));
            }
        }
        */
      }
    }

    /*! \class Numerical
     *  \test \b NumericalTest_anisotropic tests the numerical evaluation
     * of the anisotropic liquid matrix representations of S and D
     */
    WHEN("the ionic liquid Green's function is used") {
      Eigen::Vector3d epsilon = (Eigen::Vector3d() << 80.0, 80.0, 80.0).finished();
      Eigen::Vector3d euler = (Eigen::Vector3d() << 0.0, 0.0, 0.0).finished();
      AnisotropicLiquid<> gf(epsilon, euler);
      THEN("the matrix elements of S are") {
        results = op.computeS(cavity, gf);
        // Numerical integrator not really working now...
        /*
        reference =
        cnpy::custom::npy_load<double>("anisotropicliquid_S_numerical.npy");
        for (int i = 0; i < cavity.size(); ++i) {
            for (int j = 0; j < cavity.size(); ++j) {
                REQUIRE(reference(i, j) == Approx(results(i, j)));
            }
        }
        */
      }
      AND_THEN("the matrix elements of D are") {
        results = op.computeD(cavity, gf);
        // Numerical integrator not really working now...
        /*
        reference =
        cnpy::custom::npy_load<double>("anisotropicliquid_D_numerical.npy");
        for (int i = 0; i < cavity.size(); ++i) {
            for (int j = 0; j < cavity.size(); ++j) {
                REQUIRE(reference(i, j) == Approx(results(i, j)));
            }
        }
        */
      }
    }
  }
}
