/*
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

#include "Solvent.hpp"

#include <map>
#include <ostream>
#include <string>

#include "Config.hpp"

namespace pcm {
using utils::Solvent;

SolventMap & solvents() {
  static SolventMap availableSolvents;
  // ------------------------------------------------------------
  // Solvent("name", static permittivity, optical permittivity, probe radius)
  availableSolvents["N-HEPTANE"] = Solvent("N-heptane", 1.920, 1.918, 3.125);
  availableSolvents["CYCLOHEXANE"] = Solvent("Cyclohexane", 2.023, 2.028, 2.815);
  availableSolvents["CARBON TETRACHLORIDE"] =
      Solvent("Carbon tetrachloride", 2.228, 2.129, 2.685);
  availableSolvents["BENZENE"] = Solvent("Benzene", 2.247, 2.244, 2.630);
  // No data available for the probe radius of 1,4-dioxane: use that of benzene
  availableSolvents["1,4-DIOXANE"] = Solvent("1,4-Dioxane", 2.250, 2.023, 2.630);
  availableSolvents["TOLUENE"] = Solvent("Toluene", 2.379, 2.232, 2.820);
  availableSolvents["CHLOROFORM"] = Solvent("Chloroform", 4.900, 2.085, 2.480);
  availableSolvents["CHLOROBENZENE"] = Solvent("Chlorobenzene", 5.621, 2.320, 2.805);
  availableSolvents["ANILINE"] = Solvent("Aniline", 6.890, 2.506, 2.800);
  availableSolvents["TETRAHYDROFURANE"] =
      Solvent("Tetrahydrofurane", 7.580, 1.971, 2.900);
  availableSolvents["METHYLENECHLORIDE"] =
      Solvent("Methylenechloride", 8.930, 2.020, 2.270);
  availableSolvents["1,2-DICHLOROETHANE"] =
      Solvent("1,2-Dichloroethane", 10.360, 2.085, 2.505);
  availableSolvents["ACETONE"] = Solvent("Acetone", 20.700, 1.841, 2.380);
  availableSolvents["ETHANOL"] = Solvent("Ethanol", 24.550, 1.847, 2.180);
  availableSolvents["METHANOL"] = Solvent("Methanol", 32.630, 1.758, 1.855);
  availableSolvents["ACETONITRILE"] = Solvent("Acetonitrile", 36.640, 1.806, 2.155);
  availableSolvents["NITROMETHANE"] = Solvent("Nitromethane", 38.200, 1.904, 2.155);
  availableSolvents["DIMETHYLSULFOXIDE"] =
      Solvent("Dimethylsulfoxide", 46.700, 2.179, 2.455);
  // No data available for the probe radius of propylene carbonate: use that of water
  availableSolvents["PROPYLENECARBONATE"] =
      Solvent("Propylene Carbonate", 64.960, 2.019, 1.385);
  availableSolvents["WATER"] = Solvent("Water", 78.390, 1.776, 1.385);
  availableSolvents["EXPLICIT"] = Solvent("Explicit", 0.000, 0.000, 0.000);
  // ------------------------------------------------------------
  return availableSolvents;
}

std::ostream & operator<<(std::ostream & os, Solvent & solvent) {
  os << "Solvent name:          " << solvent.name << std::endl;
  os << "Static  permittivity = " << solvent.epsStatic << std::endl;
  os << "Optical permittivity = " << solvent.epsDynamic << std::endl;
  os << "Solvent radius =       " << solvent.probeRadius << " Ang";
  return os;
}
} // namespace pcm
