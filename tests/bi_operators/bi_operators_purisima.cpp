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
#include "bi_operators/IntegratorTypes.hpp"
#include "cavity/Element.hpp"
#include "cavity/GePolCavity.hpp"
#include "bi_operators/PurisimaIntegrator.hpp"
#include "green/SphericalDiffuse.hpp"
#include "TestingMolecules.hpp"
#include "green/UniformDielectric.hpp"
#include "green/Vacuum.hpp"
#include "utils/MathUtils.hpp"

SCENARIO("A collocation integrator with diagonal elements according to Purisima for D", "[bi_operators][bi_operators_purisima]")
{
    GIVEN("A GePol cavity for a single sphere in the origin")
    {
        double radius = 1.44;
        Molecule molec = dummy<0>(1.44 / bohrToAngstrom());
        double area = 10.0;
        GePolCavity cavity = GePolCavity(molec, area, 0.0, 100.0);
        Eigen::MatrixXd results = Eigen::MatrixXd::Zero(cavity.size(), cavity.size());
        Eigen::MatrixXd reference = Eigen::MatrixXd::Zero(cavity.size(), cavity.size());

        /*! \class PurisimaIntegrator
         *  \test \b PurisimaIntegratorTest_vacuum tests the evaluation by collocation of the vacuum matrix representations of S and D
         */
        WHEN("the vacuum Green's function is used")
        {
            Vacuum<AD_directional, PurisimaIntegrator> gf;
            THEN("the matrix elements of S are")
            {
                results = gf.singleLayer(cavity.elements());
                reference = cnpy::custom::npy_load<double>("vacuum_S_collocation.npy");
                for (int i = 0; i < cavity.size(); ++i) {
                    for (int j = 0; j < cavity.size(); ++j) {
                        REQUIRE(reference(i, j) == Approx(results(i, j)));
                    }
                }
            }
            AND_THEN("the matrix elements of D are")
            {
                results = gf.doubleLayer(cavity.elements());
                reference = cnpy::custom::npy_load<double>("vacuum_D_purisima.npy");
                for (int i = 0; i < cavity.size(); ++i) {
                    for (int j = 0; j < cavity.size(); ++j) {
                        REQUIRE(reference(i, j) == Approx(results(i, j)));
                    }
                }
            }
        }

        /*! \class PurisimaIntegrator
         *  \test \b PurisimaIntegratorTest_uniformdielectric tests the evaluation by collocation of the uniform dielectric matrix representations of S and D
         */
        AND_WHEN("the uniform dielectric Green's function is used")
        {
            double epsilon = 80.0;
            UniformDielectric<AD_directional, PurisimaIntegrator> gf(epsilon);
            THEN("the matrix elements of S are")
            {
                results = gf.singleLayer(cavity.elements());
                reference = cnpy::custom::npy_load<double>("uniformdielectric_S_collocation.npy");
                for (int i = 0; i < cavity.size(); ++i) {
                    for (int j = 0; j < cavity.size(); ++j) {
                        REQUIRE(reference(i, j) == Approx(results(i, j)));
                    }
                }
            }
            AND_THEN("the matrix elements of D are")
            {
                results = gf.doubleLayer(cavity.elements());
                reference = cnpy::custom::npy_load<double>("uniformdielectric_D_purisima.npy");
                for (int i = 0; i < cavity.size(); ++i) {
                    for (int j = 0; j < cavity.size(); ++j) {
                        REQUIRE(reference(i, j) == Approx(results(i, j)));
                    }
                }
            }
        }
    }
}
