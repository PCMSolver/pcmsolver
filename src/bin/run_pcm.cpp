/*
 * PCMSolver, an API for the Polarizable Continuum Model
 * Copyright (C) 2020 Roberto Di Remigio, Luca Frediani and contributors.
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

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "Config.hpp"

#include <Eigen/Core>

#include "interface/Input.hpp"
#include "interface/Meddle.hpp"
#include "utils/ChargeDistribution.hpp"
#include "utils/Molecule.hpp"

std::ofstream pcmsolver_out;

void host_writer(const char * message);

std::string remove_extension(const std::string & filename);

int main(int argc, char * argv[]) {
  if (argc > 2)
    PCMSOLVER_ERROR("Too many arguments supplied");

  // Prepare output filename
  std::string filename(remove_extension(argv[1]).erase(0, 1) + ".out");
  pcmsolver_out.open(filename.c_str());

  using namespace pcm;

  TIMER_ON("Input parsing");
  Input input(argv[1]);
  TIMER_OFF("Input parsing");
  Meddle context_(input, host_writer);

  context_.printCitation();
  context_.printInfo();

  PCMSolverIndex size = context_.getCavitySize();

  // Form vector with electrostatic potential
  // First compute the potential from the classical point multipoles distribution
  // then add the one from the molecule
  TIMER_ON("Computing MEP");
  // FIXME
  // 1. Currently hardcoded to the dipole-dipole interaction potential in vacuum
  // 2. Try to re-write this such that the input object is not needed!!!
  Eigen::VectorXd mep = Eigen::VectorXd::Zero(size);
  if (input.MEPfromChargeDistribution())
    mep += computeDipolarPotential(context_.getCenters(), input.multipoles());
  if (input.MEPfromMolecule())
    mep += computeMEP(context_.molecule(), context_.getCenters());
  TIMER_OFF("Computing MEP");
  context_.setSurfaceFunction(mep.size(), mep.data(), "MEP");
  // Compute apparent surface charge
  int irrep = 0;
  TIMER_ON("Computing ASC");
  context_.computeASC("MEP", "ASC", irrep);
  Eigen::VectorXd asc(size);
  context_.getSurfaceFunction(asc.size(), asc.data(), "ASC");
  context_.computeResponseASC("MEP", "RspASC", irrep);
  Eigen::VectorXd rsp_asc(size);
  context_.getSurfaceFunction(rsp_asc.size(), rsp_asc.data(), "RspASC");
  TIMER_OFF("Computing ASC");
  context_.printSurfaceFunction("ASC");
  // Compute energy and print it out
  pcmsolver_out << "Solvation energy = "
                << std::setprecision(std::numeric_limits<long double>::digits10)
                << context_.computePolarizationEnergy("MEP", "ASC") << std::endl;
  pcmsolver_out << "DONE!" << std::endl;

  pcmsolver_out.close();
  // Write timings out
  context_.writeTimings();
  TIMER_DONE("pcmsolver.timer.dat");

  return EXIT_SUCCESS;
}

void host_writer(const char * message) {
  pcmsolver_out << std::string(message) << std::endl;
}

std::string remove_extension(const std::string & filename) {
  size_t lastdot = filename.find_last_of(".");
  if (lastdot == std::string::npos)
    return filename;
  return filename.substr(0, lastdot);
}
