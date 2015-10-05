/* pcmsolver_copyright_start */
/*
 *     PCMSolver, an API for the Polarizable Continuum Model
 *     Copyright (C) 2013-2015 Roberto Di Remigio, Luca Frediani and contributors
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

#include "catch.hpp"

#include "Config.hpp"

#include <Eigen/Core>

#include "CPCMSolver.hpp"
#include "CollocationIntegrator.hpp"
#include "DerivativeTypes.hpp"
#include "GePolCavity.hpp"
#include "Molecule.hpp"
#include "Symmetry.hpp"
#include "TestingMolecules.hpp"
#include "UniformDielectric.hpp"
#include "Vacuum.hpp"

/*! \class CPCMSolver
 *  \test \b NH3GePol tests CPCMSolver using ammonia and a GePol cavity
 */
TEST_CASE("Test solver for the C-PCM with NH3 molecule and a GePol cavity", "[solver][cpcm][cpcm_gepol-NH3]")
{
    Molecule molecule = NH3();

    double area = 0.4;
    double probeRadius = 0.0;
    double minRadius = 100.0;
    GePolCavity cavity = GePolCavity(molecule, area, probeRadius, minRadius);
    cavity.saveCavity("nh3.npz");

    double permittivity = 78.39;
    Vacuum<AD_directional, CollocationIntegrator> gfInside = Vacuum<AD_directional, CollocationIntegrator>();
    UniformDielectric<AD_directional, CollocationIntegrator> gfOutside =
    UniformDielectric<AD_directional, CollocationIntegrator>(permittivity);
    bool symm = true;
    double correction = 0.8;
    CPCMSolver solver(symm, correction);
    solver.buildSystemMatrix(cavity, gfInside, gfOutside);

    double Ncharge = 7.0;
    double Hcharge = 1.0;
    size_t size = cavity.size();
    Eigen::VectorXd fake_mep = computeMEP(molecule, cavity.elements());
    // The total ASC for a conductor is -Q
    // for CPCM it will be -Q*[(epsilon-1)/epsilon + correction]
    Eigen::VectorXd fake_asc = Eigen::VectorXd::Zero(size);
    fake_asc = solver.computeCharge(fake_mep);
    double totalASC = - (Ncharge + 3.0 * Hcharge) * (permittivity - 1) /
                      (permittivity + correction);
    double totalFakeASC = fake_asc.sum();
    CAPTURE(totalASC - totalFakeASC);
    REQUIRE(totalASC == Approx(totalFakeASC).epsilon(1.0e-03));
}
