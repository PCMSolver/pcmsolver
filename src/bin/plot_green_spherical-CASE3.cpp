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

#include <iostream>
#include <ostream>
#include <fstream>

#include <Eigen/Core>

#include "CollocationIntegrator.hpp"
#include "DerivativeTypes.hpp"
#include "UniformDielectric.hpp"
#include "SphericalDiffuse.hpp"
#include "LoggerInterface.hpp"

int main()
{
    int nPoints = 10000;
    double epsInside = 80.0;
    double epsOutside = 2.0;
    Eigen::Vector3d sphereCenter;
    sphereCenter << 0.0, 0.0, 0.0;
    double sphereRadius = 100.0;
    double width = 10.0;
    int maxL = 30;

    Eigen::Vector3d source;
    source << 4.0, 0.0, 100.0;
    Eigen::Vector3d probe;
    probe << 0.0, 0.0, 0.0;

    double zMin =  10.0;
    double zMax = 300.0;
    double step = (zMax - zMin) / nPoints;

    SphericalDiffuse<CollocationIntegrator, OneLayerTanh> gf(epsInside, epsOutside, width, sphereRadius, sphereCenter, maxL);
    LOG(gf);
    std::ofstream out;
    out.open("gf_spherical_CASE1.dat");
    out << "#" << '\t' << "Distance" << '\t' << "gf_value" << '\t' << "image" << '\t' << "Coulomb"
        << '\t' << "derivativeProbe" << '\t' << "der_image" << '\t' << "der_Coulomb" << std::endl;
    out.precision(16);
    for (int i = 0; i < nPoints; ++i) {
        probe(2) = zMin + i*step;
        out << '\t' << probe(2)
            << '\t' << gf.kernelS(source, probe)
            << '\t' << gf.imagePotential(source, probe)
            << '\t' << gf.Coulomb(source, probe)
            << '\t' << gf.derivativeProbe(Eigen::Vector3d::UnitZ(), source, probe)
            << '\t' << gf.CoulombDerivative(Eigen::Vector3d::UnitZ(), source, probe)
            << '\t' << gf.imagePotentialDerivative(Eigen::Vector3d::UnitZ(), source, probe) << std::endl;
    }
    out.close();
    LOG_TIME;

    UniformDielectric<AD_directional, CollocationIntegrator> gf_inside(epsInside);
    out.open("gf_uniform_inside_CASE1.dat");
    out << "#" << '\t' << "Distance" << '\t' << "gf_value" << '\t' << "derivativeProbe" << std::endl;
    out.precision(16);
    for (int i = 0; i < nPoints; ++i) {
        probe(2) = zMin + i*step;
        out << '\t' << probe(2)
            << '\t' << gf_inside.kernelS(source, probe)
            << '\t' << gf_inside.derivativeProbe(Eigen::Vector3d::UnitZ(), source, probe) << std::endl;
    }
    out.close();

    UniformDielectric<AD_directional, CollocationIntegrator> gf_outside(epsOutside);
    out.open("gf_uniform_outside_CASE1.dat");
    out << "#" << '\t' << "Distance" << '\t' << "gf_value" << '\t' << "derivativeProbe" << std::endl;
    out.precision(16);
    for (int i = 0; i < nPoints; ++i) {
        probe(2) = zMin + i*step;
        out << '\t' << probe(2)
            << '\t' << gf_outside.kernelS(source, probe)
            << '\t' << gf_outside.derivativeProbe(Eigen::Vector3d::UnitZ(), source, probe) << std::endl;
    }
    out.close();
}
