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
#include "bi_operators/Collocation.hpp"
#include "cavity/Element.hpp"
#include "cavity/GePolCavity.hpp"
#include "green/SphericalDiffuse.hpp"
#include "green/UniformDielectric.hpp"
#include "green/Vacuum.hpp"
#include "utils/MathUtils.hpp"

using namespace pcm;
using bi_operators::Collocation;
using cavity::GePolCavity;
using green::Vacuum;
using green::UniformDielectric;
using green::SphericalDiffuse;

SCENARIO("A collocation integrator with approximate diagonal elements",
         "[bi_operators][bi_operators_collocation]") {
  GIVEN("A GePol cavity for a single sphere in the origin") {
    double radius = 1.44;
    Molecule molec = dummy<0>(1.44 / bohrToAngstrom());
    double area = 10.0;
    GePolCavity cavity = GePolCavity(molec, area, 0.0, 100.0);
    Eigen::MatrixXd results = Eigen::MatrixXd::Zero(cavity.size(), cavity.size());
    Eigen::MatrixXd reference = Eigen::MatrixXd::Zero(cavity.size(), cavity.size());
    Collocation op;

    /*! \class Collocation
     *  \test \b CollocationTest_vacuum tests the evaluation by collocation
     * of the vacuum matrix representations of S and D
     */
    WHEN("the vacuum Green's function is used") {
      Vacuum<> gf;
      THEN("the matrix elements of S are") {
        results = op.computeS(cavity, gf);
        reference = cnpy::custom::npy_load<double>("vacuum_S_collocation.npy");
        for (int i = 0; i < cavity.size(); ++i) {
          for (int j = 0; j < cavity.size(); ++j) {
            REQUIRE(reference(i, j) == Approx(results(i, j)));
          }
        }
      }
      AND_THEN("the matrix elements of D are") {
        results = op.computeD(cavity, gf);
        reference = cnpy::custom::npy_load<double>("vacuum_D_collocation.npy");
        for (int i = 0; i < cavity.size(); ++i) {
          for (int j = 0; j < cavity.size(); ++j) {
            REQUIRE(reference(i, j) == Approx(results(i, j)));
          }
        }
      }
    }

    /*! \class Collocation
     *  \test \b CollocationTest_uniformdielectric tests the evaluation by
     * collocation of the uniform dielectric matrix representations of S and D
     */
    AND_WHEN("the uniform dielectric Green's function is used") {
      double epsilon = 80.0;
      UniformDielectric<> gf(epsilon);
      THEN("the matrix elements of S are") {
        results = op.computeS(cavity, gf);
        reference =
            cnpy::custom::npy_load<double>("uniformdielectric_S_collocation.npy");
        for (int i = 0; i < cavity.size(); ++i) {
          for (int j = 0; j < cavity.size(); ++j) {
            REQUIRE(reference(i, j) == Approx(results(i, j)));
          }
        }
      }
      AND_THEN("the matrix elements of D are") {
        results = op.computeD(cavity, gf);
        reference =
            cnpy::custom::npy_load<double>("uniformdielectric_D_collocation.npy");
        for (int i = 0; i < cavity.size(); ++i) {
          for (int j = 0; j < cavity.size(); ++j) {
            REQUIRE(reference(i, j) == Approx(results(i, j)));
          }
        }
      }
    }

    /*! \class Collocation
     *  \test \b CollocationTest_tanhsphericaldiffuse tests the evaluation
     * by collocation of the spherical diffuse matrix representations of S and D
     */
    AND_WHEN("the spherical diffuse with a tanh profile Green's function is used") {
      double epsilon = 80.0;
      double width = 5.0;
      double sphereRadius = 100.0;
      SphericalDiffuse<> gf(
          epsilon, epsilon, width, sphereRadius, Eigen::Vector3d::Zero(), 3);
      THEN("the matrix elements of S are") {
        results = op.computeS(cavity, gf);
        reference =
            cnpy::custom::npy_load<double>("tanhsphericaldiffuse_S_collocation.npy");
        for (int i = 0; i < cavity.size(); ++i) {
          for (int j = 0; j < cavity.size(); ++j) {
            REQUIRE(reference(i, j) == Approx(results(i, j)));
          }
        }
      }
      AND_THEN("the matrix elements of D are") {
        results = op.computeD(cavity, gf);
        reference =
            cnpy::custom::npy_load<double>("tanhsphericaldiffuse_D_collocation.npy");
        for (int i = 0; i < cavity.size(); ++i) {
          for (int j = 0; j < cavity.size(); ++j) {
            REQUIRE(reference(i, j) == Approx(results(i, j)));
          }
        }
      }
    }
  }
}
