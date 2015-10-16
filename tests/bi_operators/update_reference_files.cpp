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

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <limits>

#include "Config.hpp"

#include <Eigen/Core>

#include "cnpy.hpp"
#include "CollocationIntegrator.hpp"
#include "PurisimaIntegrator.hpp"
#include "NumericalIntegrator.hpp"
#include "DerivativeTypes.hpp"
#include "Element.hpp"
#include "GePolCavity.hpp"
#include "PhysicalConstants.hpp"
#include "SphericalDiffuse.hpp"
#include "TestingMolecules.hpp"
#include "AnisotropicLiquid.hpp"
#include "IonicLiquid.hpp"
#include "SphericalDiffuse.hpp"
#include "UniformDielectric.hpp"
#include "Vacuum.hpp"

void save_vacuum_collocation();
void save_uniform_dielectric_collocation();
void save_tanh_spherical_diffuse_collocation();

void save_vacuum_purisima();
void save_uniform_dielectric_purisima();

void save_vacuum_numerical();
void save_uniform_dielectric_numerical();
void save_ionic_liquid_numerical();
void save_anisotropic_liquid_numerical();
void save_tanh_spherical_diffuse_numerical();

int main()
{
    save_vacuum_collocation();
    save_uniform_dielectric_collocation();
    save_tanh_spherical_diffuse_collocation();

    save_vacuum_purisima();
    save_uniform_dielectric_purisima();

    return EXIT_SUCCESS;
}

void save_vacuum_collocation() {
    Molecule molec = dummy<0>(1.44 / convertBohrToAngstrom);
    double area = 10.0;
    GePolCavity cavity(molec, area, 0.0, 100.0);

    unsigned int dim = static_cast<unsigned int>(cavity.size());
    const unsigned int shape[] = {dim, dim};

    Vacuum<AD_directional, CollocationIntegrator> gf;

    Eigen::MatrixXd S_results = gf.singleLayer(cavity.elements());
    cnpy::npy_save("vacuum_S_collocation.npy", S_results.data(), shape, 2, "w", true);

    Eigen::MatrixXd D_results = gf.doubleLayer(cavity.elements());
    cnpy::npy_save("vacuum_D_collocation.npy", D_results.data(), shape, 2, "w", true);
}

void save_uniform_dielectric_collocation() {
    double epsilon = 80.0;
    Molecule molec = dummy<0>(1.44 / convertBohrToAngstrom);
    double area = 10.0;
    GePolCavity cavity(molec, area, 0.0, 100.0);

    unsigned int dim = static_cast<unsigned int>(cavity.size());
    const unsigned int shape[] = {dim, dim};

    UniformDielectric<AD_directional, CollocationIntegrator> gf(epsilon);

    Eigen::MatrixXd S_results = gf.singleLayer(cavity.elements());
    cnpy::npy_save("uniformdielectric_S_collocation.npy", S_results.data(), shape, 2, "w", true);

    Eigen::MatrixXd D_results = gf.doubleLayer(cavity.elements());
    cnpy::npy_save("uniformdielectric_D_collocation.npy", D_results.data(), shape, 2, "w", true);
}

void save_tanh_spherical_diffuse_collocation() {
    double epsilon = 80.0;
    double width = 5.0;
    double sphereRadius = 100.0;
    Molecule molec = dummy<0>(1.44 / convertBohrToAngstrom);
    double area = 10.0;
    GePolCavity cavity(molec, area, 0.0, 100.0);

    unsigned int dim = static_cast<unsigned int>(cavity.size());
    const unsigned int shape[] = {dim, dim};

    SphericalDiffuse<CollocationIntegrator, OneLayerTanh> gf(epsilon, epsilon, width, sphereRadius, Eigen::Vector3d::Zero(), 3);

    Eigen::MatrixXd S_results = gf.singleLayer(cavity.elements());
    cnpy::npy_save("tanhsphericaldiffuse_S_collocation.npy", S_results.data(), shape, 2, "w", true);

    Eigen::MatrixXd D_results = gf.doubleLayer(cavity.elements());
    cnpy::npy_save("tanhsphericaldiffuse_D_collocation.npy", D_results.data(), shape, 2, "w", true);
}

void save_vacuum_purisima() {
    Molecule molec = dummy<0>(1.44 / convertBohrToAngstrom);
    double area = 10.0;
    GePolCavity cavity(molec, area, 0.0, 100.0);

    unsigned int dim = static_cast<unsigned int>(cavity.size());
    const unsigned int shape[] = {dim, dim};

    Vacuum<AD_directional, PurisimaIntegrator> gf;

    Eigen::MatrixXd D_results = gf.doubleLayer(cavity.elements());
    cnpy::npy_save("vacuum_D_purisima.npy", D_results.data(), shape, 2, "w", true);
}

void save_uniform_dielectric_purisima() {
    double epsilon = 80.0;
    Molecule molec = dummy<0>(1.44 / convertBohrToAngstrom);
    double area = 10.0;
    GePolCavity cavity(molec, area, 0.0, 100.0);

    unsigned int dim = static_cast<unsigned int>(cavity.size());
    const unsigned int shape[] = {dim, dim};

    UniformDielectric<AD_directional, PurisimaIntegrator> gf(epsilon);

    Eigen::MatrixXd D_results = gf.doubleLayer(cavity.elements());
    cnpy::npy_save("uniformdielectric_D_purisima.npy", D_results.data(), shape, 2, "w", true);
}

void save_vacuum_numerical() {
    Molecule molec = dummy<0>(1.44 / convertBohrToAngstrom);
    double area = 10.0;
    GePolCavity cavity(molec, area, 0.0, 100.0);

    unsigned int dim = static_cast<unsigned int>(cavity.size());
    const unsigned int shape[] = {dim, dim};

    Vacuum<AD_directional, NumericalIntegrator> gf;

    Eigen::MatrixXd S_results = gf.singleLayer(cavity.elements());
    cnpy::npy_save("vacuum_S_numerical.npy", S_results.data(), shape, 2, "w", true);

    Eigen::MatrixXd D_results = gf.doubleLayer(cavity.elements());
    cnpy::npy_save("vacuum_D_numerical.npy", D_results.data(), shape, 2, "w", true);
}

void save_uniform_dielectric_numerical() {
    double epsilon = 80.0;
    Molecule molec = dummy<0>(1.44 / convertBohrToAngstrom);
    double area = 10.0;
    GePolCavity cavity(molec, area, 0.0, 100.0);

    unsigned int dim = static_cast<unsigned int>(cavity.size());
    const unsigned int shape[] = {dim, dim};

    UniformDielectric<AD_directional, NumericalIntegrator> gf(epsilon);

    Eigen::MatrixXd S_results = gf.singleLayer(cavity.elements());
    cnpy::npy_save("uniformdielectric_S_numerical.npy", S_results.data(), shape, 2, "w", true);

    Eigen::MatrixXd D_results = gf.doubleLayer(cavity.elements());
    cnpy::npy_save("uniformdielectric_D_numerical.npy", D_results.data(), shape, 2, "w", true);
}

void save_ionic_liquid_numerical() {
    double eps = 80.0;
    double kappa = 1.5;
    Molecule molec = dummy<0>(1.44 / convertBohrToAngstrom);
    double area = 10.0;
    GePolCavity cavity(molec, area, 0.0, 100.0);

    unsigned int dim = static_cast<unsigned int>(cavity.size());
    const unsigned int shape[] = {dim, dim};

    IonicLiquid<AD_directional, NumericalIntegrator> gf(eps, kappa);

    Eigen::MatrixXd S_results = gf.singleLayer(cavity.elements());
    cnpy::npy_save("ionicliquid_S_numerical.npy", S_results.data(), shape, 2, "w", true);

    Eigen::MatrixXd D_results = gf.doubleLayer(cavity.elements());
    cnpy::npy_save("ionicliquid_D_numerical.npy", D_results.data(), shape, 2, "w", true);
}

void save_anisotropic_liquid_numerical() {
    Eigen::Vector3d epsilon = Eigen::Vector3d::Zero();
    epsilon << 80.0, 80.0, 80.0;
    Eigen::Vector3d euler = Eigen::Vector3d::Zero();
    euler << 0.0, 0.0, 0.0;
    Molecule molec = dummy<0>(1.44 / convertBohrToAngstrom);
    double area = 10.0;
    GePolCavity cavity(molec, area, 0.0, 100.0);

    unsigned int dim = static_cast<unsigned int>(cavity.size());
    const unsigned int shape[] = {dim, dim};

    AnisotropicLiquid<AD_directional, NumericalIntegrator> gf(epsilon, euler);

    Eigen::MatrixXd S_results = gf.singleLayer(cavity.elements());
    cnpy::npy_save("anisotropicliquid_S_numerical.npy", S_results.data(), shape, 2, "w", true);

    Eigen::MatrixXd D_results = gf.doubleLayer(cavity.elements());
    cnpy::npy_save("anisotropicliquid_D_numerical.npy", D_results.data(), shape, 2, "w", true);
}

void save_tanh_spherical_diffuse_numerical() {
    double epsilon = 80.0;
    double width = 5.0;
    double sphereRadius = 100.0;
    Molecule molec = dummy<0>(1.44 / convertBohrToAngstrom);
    double area = 10.0;
    GePolCavity cavity(molec, area, 0.0, 100.0);

    unsigned int dim = static_cast<unsigned int>(cavity.size());
    const unsigned int shape[] = {dim, dim};

    SphericalDiffuse<NumericalIntegrator, OneLayerTanh> gf(epsilon, epsilon, width, sphereRadius, Eigen::Vector3d::Zero(), 3);

    Eigen::MatrixXd S_results = gf.singleLayer(cavity.elements());
    cnpy::npy_save("tanhsphericaldiffuse_S_numerical.npy", S_results.data(), shape, 2, "w", true);

    Eigen::MatrixXd D_results = gf.doubleLayer(cavity.elements());
    cnpy::npy_save("tanhsphericaldiffuse_D_numerical.npy", D_results.data(), shape, 2, "w", true);
}
