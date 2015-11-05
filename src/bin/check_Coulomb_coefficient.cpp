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
 *     PCMSolver API, see: <http://pcmsolver.readthedocs.org/>
 */
/* pcmsolver_copyright_end */

#include <iostream>
#include <ostream>
#include <fstream>

#include <Eigen/Core>

#include "DerivativeTypes.hpp"
#include "UniformDielectric.hpp"
#include "TanhSphericalDiffuse.hpp"

int main()
{
    int nPoints = 10000;
    double epsInside = 1.0;
    double epsOutside = 80.0;
    Eigen::Vector3d sphereCenter;
    sphereCenter << 0.0, 0.0, 0.0;
    double sphereRadius = 100.0;
    double width = 10.0;

    Eigen::Vector3d source;
    source << 0.0, 0.0, 0.0;
    Eigen::Vector3d probe;
    double delta = 0.05;

    double zMin = 5.0;
    double zMax = 200.0;
    double step = (zMax - zMin) / nPoints;

    TanhSphericalDiffuse gf(epsInside, epsOutside, width, sphereRadius, sphereCenter);
    // Plot Green's function, image potential and separation coefficient
    std::ofstream out, out_err;
    out.open("checkCoulomb_spherical.log");
    out << "#" << '\t' << "Distance" << '\t' << "coefficient" << std::endl;
    out.precision(16);
    out_err.open("checkCoulomb_diff.log");
    out_err << "#" << '\t' << "Distance" << '\t' << "difference coefficient" << std::endl;
    out_err.precision(16);
    for (int i = 0; i < nPoints; ++i) {
        source(2) = zMin + i*step;
        probe = source;
        probe(2) +=delta;
        double s = source.norm();
        double p = probe.norm();
        double ave = (s + p) / 2;
        double eps_ave = 0.0, epsPrime_ave = 0.0;
        gf.epsilon(eps_ave, epsPrime_ave, ave);
        double coeff = gf.coefficientCoulomb(source, probe);
        double coulomb = gf.Coulomb(source, probe);
        double ave_coulomb = 1.0 / (eps_ave * (source-probe).norm());
        out << '\t' << source(2) << '\t' << coulomb << '\t' << ave_coulomb << std::endl;
        out_err << '\t' << source(2) << '\t' << coeff << '\t' << eps_ave << '\t' << coeff-eps_ave << std::endl;
    }
    out.close();
    out_err.close();

    UniformDielectric<AD_directional> gf_inside(epsInside);
    out.open("checkCoulomb_uniform_inside.log");
    out << "#" << '\t' << "Distance" << '\t' << "gf_value" << std::endl;
    out.precision(16);
    for (int i = 0; i < nPoints; ++i) {
        source(2) = zMin + i*step;
        probe = source;
        probe(2) +=delta;
        out << '\t' << source(2) << '\t' << gf_inside.function(source, probe) << std::endl;
    }
    out.close();

    UniformDielectric<AD_directional> gf_outside(epsOutside);
    out.open("checkCoulomb_uniform_outside.log");
    out << "#" << '\t' << "Distance" << '\t' << "gf_value" << std::endl;
    out.precision(16);
    for (int i = 0; i < nPoints; ++i) {
        source(2) = zMin + i*step;
        probe = source;
        probe(2) +=delta;
        out << '\t' << source(2) << '\t' << gf_outside.function(source, probe) << std::endl;
    }
    out.close();
}
