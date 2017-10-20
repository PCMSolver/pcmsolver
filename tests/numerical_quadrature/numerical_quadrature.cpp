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
#include <limits>
#include <sstream>

#include <Eigen/Core>

#include "TestingMolecules.hpp"
#include "bi_operators/Numerical.hpp"
#include "cavity/Element.hpp"
#include "cavity/GePolCavity.hpp"
#include "utils/MathUtils.hpp"
#include "utils/cnpy.hpp"

using namespace pcm;
using bi_operators::integrateS;
using bi_operators::integrateD;
using cavity::GePolCavity;

double constant(const Eigen::Vector3d & /* s */, const Eigen::Vector3d & /* p */);
double one_over_r(double r,
                  const Eigen::Vector3d & /* s */,
                  const Eigen::Vector3d & /* p */);

SCENARIO("Numerical quadrature of functions", "[numerical_quadrature]") {
  GIVEN("A function on the unit sphere") {
    double radius = 1.55;
    double area = 0.4;
    Molecule point = dummy<0>(radius);
    GePolCavity cavity(point, area, 0.0, 100.0, "");
    Eigen::VectorXd results = Eigen::VectorXd::Zero(cavity.size());
    /*! \class NumericalQuadrature
     *  \test \b NumericalQuadrature_sphere tests numerical quadrature on a sphere
     * for integrand = 1.0
     */
    WHEN("The function is a constant") {
      THEN("the integrals are the finite element areas") {
        for (int i = 0; i < cavity.size(); ++i) {
          results(i) = integrateS<32, 16>(pcm::bind(constant, pcm::_1, pcm::_2),
                                          cavity.elements(i));
          double diff = results(i) - cavity.elementArea(i);
          if (std::abs(diff) > 1.0e-12) {
            WARN("Test versus area for single sphere");
            WARN("Tessera n. " << i + 1);
            WARN("diff = " << results(i) - cavity.elementArea(i));
          }
        }
        for (int i = 0; i < cavity.size(); ++i) {
          REQUIRE(results(i) == Approx(cavity.elementArea(i)));
        }
      }
    }

    /*! \class NumericalQuadrature
     *  \test \b NumericalQuadrature_sphere_1r tests numerical quadrature on a sphere
     * for integrand 1.0/r
     */
    WHEN("The function is 1.0/r") {
      THEN("the integrals are the finite element areas divided by the sphere "
           "radius") {
        for (int i = 0; i < cavity.size(); ++i) {
          results(i) = integrateS<32, 16>(
              pcm::bind(one_over_r, 1.55, pcm::_1, pcm::_2), cavity.elements(i));
          double diff = results(i) - (cavity.elementArea(i) / radius);
          if (std::abs(diff) > 1.0e-11) {
            WARN("Test versus area divided by radius for single sphere");
            WARN("Tessera n. " << i + 1);
            WARN("diff = " << results(i) - (cavity.elementArea(i) / radius));
          }
        }
        for (int i = 0; i < cavity.size(); ++i) {
          REQUIRE(results(i) == Approx(cavity.elementArea(i) / radius));
        }
      }
    }
  }

  GIVEN("A function on a molecular surface") {
    double area = 0.2;
    double probeRadius = 1.385;
    double minRadius = 0.2;
    Molecule molecule = H2();
    GePolCavity cavity(molecule, area, probeRadius, minRadius, "");
    Eigen::VectorXd results = Eigen::VectorXd::Zero(cavity.size());
    Eigen::VectorXd reference = Eigen::VectorXd::Zero(cavity.size());

    /*! \class NumericalQuadrature
     *  \test \b NumericalQuadrature_molecule tests numerical quadrature function on
     * H2 cavity
     */
    WHEN("The function is a constant") {
      THEN("the integrals are the finite element areas") {
        for (int i = 0; i < cavity.size(); ++i) {
          results(i) = integrateS<64, 16>(pcm::bind(constant, pcm::_1, pcm::_2),
                                          cavity.elements(i));
          double diff = results(i) - cavity.elementArea(i);
          if (std::abs(diff) > 1.0e-11) {
            WARN("Test versus area for H2 molecule");
            WARN("Tessera n. " << i + 1);
            WARN("diff = " << results(i) - cavity.elementArea(i));
          }
        }
        /*
        // In case you need to update the reference files...
        cnpy::custom::npy_save("molecule.npy", results);
        */
        reference = cnpy::custom::npy_load<double>("molecule.npy");

        for (int i = 0; i < cavity.size(); ++i) {
          REQUIRE(results(i) == Approx(reference(i)));
        }
      }
    }

    /*! \class NumericalQuadrature
     *  \test \b NumericalQuadrature_molecule_1r tests numerical quadrature function
     * on H2 cavity
     */
    WHEN("The function is 1.0/r") {
      THEN("the integrals are the finite element areas divided by the sphere "
           "radius") {
        for (int i = 0; i < cavity.size(); ++i) {
          results(i) = integrateS<64, 16>(
              pcm::bind(one_over_r, 1.20, pcm::_1, pcm::_2), cavity.elements(i));
          double diff =
              results(i) - (cavity.elementArea(i) / molecule.spheres(0).radius);
          if (std::abs(diff) > 1.0e-11) {
            WARN("Test versus area divided by radius for H2 molecule");
            WARN("Tessera n. " << i + 1);
            WARN("diff = " << results(i) - (cavity.elementArea(i) /
                                            molecule.spheres(0).radius));
          }
        }
        /*
        // In case you need to update the reference files...
        cnpy::custom::npy_save("molecule_1r.npy", results);
        */
        reference = cnpy::custom::npy_load<double>("molecule_1r.npy");

        for (int i = 0; i < cavity.size(); ++i) {
          REQUIRE(results(i) == Approx(reference(i)));
        }
      }
    }
  }
}

double constant(const Eigen::Vector3d & /* s */, const Eigen::Vector3d & /* p */) {
  return 1.0;
}
double one_over_r(double r,
                  const Eigen::Vector3d & /* s */,
                  const Eigen::Vector3d & /* p */) {
  return 1.0 / r;
}
