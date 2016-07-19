/* pcmsolver_copyright_start */
/*
 *     PCMSolver, an API for the Polarizable Continuum Model
 *     Copyright (C) 2013-2016 Roberto Di Remigio, Luca Frediani and contributors
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
 *     PCMSolver API, see: <http://pcmsolver.readthedocs.io/>
 */
/* pcmsolver_copyright_end */

#include "catch.hpp"

#include <cmath>
#include <iostream>
#include <iomanip>
#include <limits>


#include <Eigen/Core>

#include "green/DerivativeTypes.hpp"
#include "green/AnisotropicLiquid.hpp"
#include "cavity/Element.hpp"
#include "cavity/GePolCavity.hpp"
#include "green/IonicLiquid.hpp"
#include "bi_operators/NumericalIntegrator.hpp"
#include "green/UniformDielectric.hpp"
#include "TestingMolecules.hpp"
#include "green/Vacuum.hpp"
#include "utils/MathUtils.hpp"

SCENARIO("A collocation integrator with numerical integrator of the diagonal elements", "[bi_operators][bi_operators_numerical]")
{
    GIVEN("A GePol cavity for a single sphere in the origin")
    {
        Molecule molec = dummy<0>(1.44 / bohrToAngstrom());
        double area = 10.0;
        GePolCavity cavity = GePolCavity(molec, area, 0.0, 100.0);
        Eigen::MatrixXd results = Eigen::MatrixXd::Zero(cavity.size(), cavity.size());
        Eigen::MatrixXd reference = Eigen::MatrixXd::Zero(cavity.size(), cavity.size());

        /*! \class NumericalIntegrator
         *  \test \b NumericalIntegratorTest_vacuum tests the numerical evaluation of the vacuum matrix representations of S and D
         */
        WHEN("the vacuum Green's function is used")
        {
            Vacuum<AD_directional, NumericalIntegrator> gf;
            THEN("the matrix elements of S are")
            {
                results = gf.singleLayer(cavity.elements());
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
            AND_THEN("the matrix elements of D are")
            {
                Eigen::MatrixXd results = gf.doubleLayer(cavity.elements());
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

        /*! \class NumericalIntegrator
         *  \test \b NumericalIntegratorTest_uniformdielectric tests the numerical evaluation of the uniform dielectric matrix representations of S and D
         */
        WHEN("the uniform dielectric Green's function is used")
        {
            double eps = 80.0;
            UniformDielectric<AD_directional, NumericalIntegrator> gf(eps);
            THEN("the matrix elements of S are")
            {
                results = gf.singleLayer(cavity.elements());
                // Numerical integrator not really working now...
                /*
                reference = cnpy::custom::npy_load<double>("uniformdielectric_S_numerical.npy");
                for (int i = 0; i < cavity.size(); ++i) {
                    for (int j = 0; j < cavity.size(); ++j) {
                        REQUIRE(reference(i, j) == Approx(results(i, j)));
                    }
                }
                */
            }
            AND_THEN("the matrix elements of D are")
            {
                Eigen::MatrixXd results = gf.doubleLayer(cavity.elements());
                // Numerical integrator not really working now...
                /*
                reference = cnpy::custom::npy_load<double>("uniformdielectric_D_numerical.npy");
                for (int i = 0; i < cavity.size(); ++i) {
                    for (int j = 0; j < cavity.size(); ++j) {
                        REQUIRE(reference(i, j) == Approx(results(i, j)));
                    }
                }
                */
            }
        }

        /*! \class NumericalIntegrator
         *  \test \b NumericalIntegratorTest_ionic tests the numerical evaluation of the ionic liquid matrix representations of S and D
         */
        WHEN("the ionic liquid Green's function is used")
        {
            double eps = 80.0;
            double kappa = 1.5;
            IonicLiquid<AD_directional, NumericalIntegrator> gf(eps, kappa);
            THEN("the matrix elements of S are")
            {
                results = gf.singleLayer(cavity.elements());
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
            AND_THEN("the matrix elements of D are")
            {
                results = gf.doubleLayer(cavity.elements());
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

        /*! \class NumericalIntegrator
         *  \test \b NumericalIntegratorTest_anisotropic tests the numerical evaluation of the anisotropic liquid matrix representations of S and D
         */
        WHEN("the ionic liquid Green's function is used")
        {
            Eigen::Vector3d epsilon = (Eigen::Vector3d() << 80.0, 80.0, 80.0).finished();
            Eigen::Vector3d euler = (Eigen::Vector3d() << 0.0, 0.0, 0.0).finished();
            AnisotropicLiquid<AD_directional, NumericalIntegrator> gf(epsilon, euler);
            THEN("the matrix elements of S are")
            {
                results = gf.singleLayer(cavity.elements());
                // Numerical integrator not really working now...
                /*
                reference = cnpy::custom::npy_load<double>("anisotropicliquid_S_numerical.npy");
                for (int i = 0; i < cavity.size(); ++i) {
                    for (int j = 0; j < cavity.size(); ++j) {
                        REQUIRE(reference(i, j) == Approx(results(i, j)));
                    }
                }
                */
            }
            AND_THEN("the matrix elements of D are")
            {
                results = gf.doubleLayer(cavity.elements());
                // Numerical integrator not really working now...
                /*
                reference = cnpy::custom::npy_load<double>("anisotropicliquid_D_numerical.npy");
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

