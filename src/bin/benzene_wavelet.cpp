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

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>

#include "Config.hpp"

#include <Eigen/Dense>

#include "DerivativeTypes.hpp"
#include "PWCSolver.hpp"
#include "UniformDielectric.hpp"
#include "Vacuum.hpp"
#include "WaveletCavity.hpp"
#include "PhysicalConstants.hpp"

void pwc_C6H6();

int main() {
    pwc_C6H6();
}

void pwc_C6H6()
{
    // Molecular geometry. These coordinates are in Angstrom!
    Eigen::Vector3d C1( 5.274,  1.999, -8.568);
    Eigen::Vector3d C2( 6.627,  2.018, -8.209);
    Eigen::Vector3d C3( 7.366,  0.829, -8.202);
    Eigen::Vector3d C4( 6.752, -0.379, -8.554);
    Eigen::Vector3d C5( 5.399, -0.398, -8.912);
    Eigen::Vector3d C6( 4.660,  0.791, -8.919);
    Eigen::Vector3d H1( 4.704,  2.916, -8.573);
    Eigen::Vector3d H2( 7.101,  2.950, -7.938);
    Eigen::Vector3d H3( 8.410,  0.844, -7.926);
    Eigen::Vector3d H4( 7.322, -1.296, -8.548);
    Eigen::Vector3d H5( 4.925, -1.330, -9.183);
    Eigen::Vector3d H6( 3.616,  0.776, -9.196);
    // Set up cavity, read it from Maharavo's file benzene2.dat
    WaveletCavity cavity("benzene2.dat");
    cavity.scaleCavity(1./convertBohrToAngstrom);

    double permittivity = 78.39;
    double Hcharge = 1.0;
    double Ccharge = 6.0;
    double totalASC = - (6 * Ccharge + 6 * Hcharge) * ( permittivity - 1) / permittivity; 

    Vacuum<AD_directional> * gfInside = new Vacuum<AD_directional>();
    UniformDielectric<AD_directional> * gfOutside = new
    UniformDielectric<AD_directional>(permittivity);
    int firstKind = 0;
#ifdef DEBUG
    FILE* debugFile = fopen("debug.out","w");
    fclose(debugFile);
#endif
    PWCSolver solver(gfInside, gfOutside, firstKind);
    solver.buildSystemMatrix(cavity);
    cavity.uploadPoints(solver.getQuadratureLevel(), solver.getT_());

    int size = cavity.size();
    Eigen::VectorXd fake_mep = Eigen::VectorXd::Zero(size);
    for (int i = 0; i < size; ++i) {
        Eigen::Vector3d center = cavity.elementCenter(i);
        double C1mep = Ccharge/(center - C1/convertBohrToAngstrom).norm();
        double C2mep = Ccharge/(center - C2/convertBohrToAngstrom).norm();
        double C3mep = Ccharge/(center - C3/convertBohrToAngstrom).norm();
        double C4mep = Ccharge/(center - C4/convertBohrToAngstrom).norm();
        double C5mep = Ccharge/(center - C5/convertBohrToAngstrom).norm();
        double C6mep = Ccharge/(center - C6/convertBohrToAngstrom).norm();

        double H1mep = Hcharge/(center - H1/convertBohrToAngstrom).norm();
        double H2mep = Hcharge/(center - H2/convertBohrToAngstrom).norm();
        double H3mep = Hcharge/(center - H3/convertBohrToAngstrom).norm();
        double H4mep = Hcharge/(center - H4/convertBohrToAngstrom).norm();
        double H5mep = Hcharge/(center - H5/convertBohrToAngstrom).norm();
        double H6mep = Hcharge/(center - H6/convertBohrToAngstrom).norm();
        fake_mep(i) = C1mep + C2mep + C3mep + C4mep + C5mep + C6mep +
          H1mep + H2mep + H3mep + H4mep + H5mep + H6mep;
    }
    // The total ASC for a dielectric is -Q*[(epsilon-1)/epsilon]
    Eigen::VectorXd fake_asc = Eigen::VectorXd::Zero(size);
    solver.compCharge(fake_mep, fake_asc);
    double totalFakeASC = fake_asc.sum();

    double energy = 0.5 * (fake_mep.dot(fake_asc));
    
    std::ofstream report;
    report.open("pwc_C6H6_report.out", std::ios::out);
    report << " Piecewise constant wavelet solver, C6H6 molecule " << std::endl;
    report << cavity << std::endl;
    report << "------------------------------------------------------------" << std::endl;
    report << "totalASC     = " << std::setprecision(20) << totalASC     << std::endl;
    report << "totalFakeASC = " << std::setprecision(20) << totalFakeASC << std::endl;
    report << "Delta        = " << std::setprecision(20) << totalASC - totalFakeASC << std::endl;
    report << "Energy       = " << std::setprecision(20) << energy << std::endl;
    report << "------------------------------------------------------------" << std::endl;
    report.close();
}
