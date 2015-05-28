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

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <limits>

#include "Config.hpp"

#include <Eigen/Dense>

#include "cnpyPimpl.hpp"
#include "CollocationIntegrator.hpp"
#include "DerivativeTypes.hpp"
#include "Element.hpp"
#include "GePolCavity.hpp"
#include "PhysicalConstants.hpp"
#include "SphericalDiffuse.hpp"
#include "TestingMolecules.hpp"
#include "UniformDielectric.hpp"
#include "Vacuum.hpp"

void save_vacuum();
void save_uniform_dielectric();
void save_tanh_spherical_diffuse();

int main()
{
    save_vacuum();
    save_uniform_dielectric();
    save_tanh_spherical_diffuse();

    return EXIT_SUCCESS;
}

void save_vacuum() {
    Molecule molec = dummy<0>(1.44 / convertBohrToAngstrom);
    double area = 10.0;
    GePolCavity cavity(molec, area, 0.0, 100.0);

    unsigned int dim = static_cast<unsigned int>(cavity.size());
    const unsigned int shape[] = {dim};

    Vacuum<AD_directional, CollocationIntegrator<AD_directional, Uniform> > gf;

    Eigen::VectorXd S_results = Eigen::VectorXd::Zero(cavity.size());
    for (int i = 0; i < cavity.size(); ++i) {
        S_results(i) = gf.diagonalS(cavity.elements(i));
    }
    cnpy::npy_save("vacuum_S_collocation.npy", S_results.data(), shape, 1, "w", false);

    Eigen::VectorXd D_results = Eigen::VectorXd::Zero(cavity.size());
    for (int i = 0; i < cavity.size(); ++i) {
        D_results(i) = gf.diagonalD(cavity.elements(i));
    }
    cnpy::npy_save("vacuum_D_collocation.npy", D_results.data(), shape, 1, "w", false);
}

void save_uniform_dielectric() {
    double epsilon = 80.0;
    Molecule molec = dummy<0>(1.44 / convertBohrToAngstrom);
    double area = 10.0;
    GePolCavity cavity(molec, area, 0.0, 100.0);

    unsigned int dim = static_cast<unsigned int>(cavity.size());
    const unsigned int shape[] = {dim};

    UniformDielectric<AD_directional, CollocationIntegrator<AD_directional, Uniform> > gf(epsilon);

    Eigen::VectorXd S_results = Eigen::VectorXd::Zero(cavity.size());
    for (int i = 0; i < cavity.size(); ++i) {
        S_results(i) = gf.diagonalS(cavity.elements(i));
    }
    cnpy::npy_save("uniformdielectric_S_collocation.npy", S_results.data(), shape, 1, "w", false);

    Eigen::VectorXd D_results = Eigen::VectorXd::Zero(cavity.size());
    for (int i = 0; i < cavity.size(); ++i) {
        D_results(i) = gf.diagonalD(cavity.elements(i));
    }
    cnpy::npy_save("uniformdielectric_D_collocation.npy", D_results.data(), shape, 1, "w", false);
}

void save_tanh_spherical_diffuse() {
    double epsilon = 80.0;
    double width = 5.0;
    double sphereRadius = 100.0;
    Molecule molec = dummy<0>(1.44 / convertBohrToAngstrom);
    double area = 10.0;
    GePolCavity cavity(molec, area, 0.0, 100.0);

    unsigned int dim = static_cast<unsigned int>(cavity.size());
    const unsigned int shape[] = {dim};

    SphericalDiffuse<CollocationIntegrator<Numerical, TanhDiffuse>, TanhDiffuse> gf(epsilon, epsilon, width, sphereRadius, Eigen::Vector3d::Zero());

    Eigen::VectorXd S_results = Eigen::VectorXd::Zero(cavity.size());
    for (int i = 0; i < cavity.size(); ++i) {
        S_results(i) = gf.diagonalS(cavity.elements(i));
    }
    cnpy::npy_save("tanhsphericaldiffuse_S_collocation.npy", S_results.data(), shape, 1, "w", false);

    Eigen::VectorXd D_results = Eigen::VectorXd::Zero(cavity.size());
    for (int i = 0; i < cavity.size(); ++i) {
        D_results(i) = gf.diagonalD(cavity.elements(i));
    }
    cnpy::npy_save("tanhsphericaldiffuse_D_collocation.npy", D_results.data(), shape, 1, "w", false);
}
