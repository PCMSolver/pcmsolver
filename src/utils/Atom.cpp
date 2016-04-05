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

#include "Atom.hpp"

#include <vector>

#include "Config.hpp"

#include <Eigen/Core>

bool invalid(const Atom & atom)
{
  return  ((std::abs(atom.radius) <= 1.0e-14) ? true : false);
}

std::vector<Atom> & initBondi()
{

    static std::vector<Atom> Bondi(92);
    Eigen::Vector3d Origin = Eigen::Vector3d::Zero();

// ------------------------------------------------------------
//              Atom("Name",    "Symbol", Charge,     Mass,   Radius, Center, Scaling)
    Bondi[0]  = Atom("Hydrogen",     "H",   1.0,   1.0078250, 1.20, Origin, 1.2);
    Bondi[1]  = Atom("Helium",       "He",  2.0,   4.0026030, 1.40, Origin, 1.2);
    Bondi[2]  = Atom("Lithium",      "Li",  3.0,   7.0160050, 1.82, Origin, 1.2);
    Bondi[3]  = Atom("Beryllium",    "Be",  4.0,   9.0121830, 1.53, Origin, 1.2);
    Bondi[4]  = Atom("Boron",        "B",   5.0,  11.0093050, 1.92, Origin, 1.2);
    Bondi[5]  = Atom("Carbon",       "C",   6.0,  12.0000000, 1.70, Origin, 1.2);
    Bondi[6]  = Atom("Nitrogen",     "N",   7.0,  14.0030740, 1.55, Origin, 1.2);
    Bondi[7]  = Atom("Oxygen",       "O",   8.0,  15.9949150, 1.52, Origin, 1.2);
    Bondi[8]  = Atom("Fluorine",     "F",   9.0,  18.9984030, 1.47, Origin, 1.2);
    Bondi[9]  = Atom("Neon",         "Ne", 10.0,  19.9924390, 1.54, Origin, 1.2);
    Bondi[10] = Atom("Sodium",       "Na", 11.0,  22.9897700, 2.27, Origin, 1.2);
    Bondi[11] = Atom("Magnesium",    "Mg", 12.0,  23.9850450, 1.73, Origin, 1.2);
    Bondi[12] = Atom("Aluminium",    "Al", 13.0,  26.9815410, 1.84, Origin, 1.2);
    Bondi[13] = Atom("Silicon",      "Si", 14.0,  27.9769280, 2.10, Origin, 1.2);
    Bondi[14] = Atom("Phosphorus",   "P",  15.0,  30.9737630, 1.80, Origin, 1.2);
    Bondi[15] = Atom("Sulphur",      "S",  16.0,  31.9720720, 1.80, Origin, 1.2);
    Bondi[16] = Atom("Chlorine",     "Cl", 17.0,  34.9688530, 1.75, Origin, 1.2);
    Bondi[17] = Atom("Argon",        "Ar", 18.0,  39.9623830, 1.88, Origin, 1.2);
    Bondi[18] = Atom("Potassium",    "K",  19.0,  38.9637080, 2.75, Origin, 1.2);
    Bondi[19] = Atom("Calcium",      "Ca", 20.0,  39.9625910, 2.31, Origin, 1.2);
    Bondi[20] = Atom("Scandium",     "Sc", 21.0,  44.9559140, 0.00, Origin, 1.2);
    Bondi[21] = Atom("Titanium",     "Ti", 22.0,  47.9479470, 0.00, Origin, 1.2);
    Bondi[22] = Atom("Vanadium",     "V",  23.0,  50.9439630, 0.00, Origin, 1.2);
    Bondi[23] = Atom("Chromium",     "Cr", 24.0,  51.9405100, 0.00, Origin, 1.2);
    Bondi[24] = Atom("Manganese",    "Mn", 25.0,  54.9380460, 0.00, Origin, 1.2);
    Bondi[25] = Atom("Iron",         "Fe", 26.0,  55.9349390, 0.00, Origin, 1.2);
    Bondi[26] = Atom("Cobalt",       "Co", 27.0,  58.9331980, 0.00, Origin, 1.2);
    Bondi[27] = Atom("Nickel",       "Ni", 28.0,  57.9353470, 1.63, Origin, 1.2);
    Bondi[28] = Atom("Copper",       "Cu", 29.0,  62.9295990, 1.40, Origin, 1.2);
    Bondi[29] = Atom("Zinc",         "Zn", 30.0,  63.9291450, 1.39, Origin, 1.2);
    Bondi[30] = Atom("Gallium",      "Ga", 31.0,  68.9255810, 1.87, Origin, 1.2);
    Bondi[31] = Atom("Germanium",    "Ge", 32.0,  73.9211790, 2.11, Origin, 1.2);
    Bondi[32] = Atom("Arsenic",      "As", 33.0,  74.9215960, 1.85, Origin, 1.2);
    Bondi[33] = Atom("Selenium",     "Se", 34.0,  79.9165210, 1.90, Origin, 1.2);
    Bondi[34] = Atom("Bromine",      "Br", 35.0,  78.9183360, 1.85, Origin, 1.2);
    Bondi[35] = Atom("Krypton",      "Kr", 36.0,  83.9115060, 2.02, Origin, 1.2);
    Bondi[36] = Atom("Rubidium",     "Rb", 37.0,  84.9118000, 3.03, Origin, 1.2);
    Bondi[37] = Atom("Strontium",    "Sr", 38.0,  87.9056250, 2.49, Origin, 1.2);
    Bondi[38] = Atom("Yttrium",      "Y",  39.0,  88.9058560, 0.00, Origin, 1.2);
    Bondi[39] = Atom("Zirconium",    "Zr", 40.0,  89.9047080, 0.00, Origin, 1.2);
    Bondi[40] = Atom("Niobium",      "Nb", 41.0,  92.9063780, 0.00, Origin, 1.2);
    Bondi[41] = Atom("Molybdenum",   "Mo", 42.0,  97.9054050, 0.00, Origin, 1.2);
    Bondi[42] = Atom("Technetium",   "Tc", 43.0,  96.0000000, 0.00, Origin, 1.2);
    Bondi[43] = Atom("Ruthenium",    "Ru", 44.0, 101.9043480, 0.00, Origin, 1.2);
    Bondi[44] = Atom("Rhodium",      "Rh", 45.0, 102.9055030, 0.00, Origin, 1.2);
    Bondi[45] = Atom("Palladium",    "Pd", 46.0, 105.9034750, 1.63, Origin, 1.2);
    Bondi[46] = Atom("Silver",       "Ag", 47.0, 106.9050950, 1.72, Origin, 1.2);
    Bondi[47] = Atom("Cadmium",      "Cd", 48.0, 113.9033610, 1.58, Origin, 1.2);
    Bondi[48] = Atom("Indium",       "In", 49.0, 114.9038750, 1.93, Origin, 1.2);
    Bondi[49] = Atom("Tin",          "Sn", 50.0, 119.9021990, 2.17, Origin, 1.2);
    Bondi[50] = Atom("Antimony",     "Sb", 51.0, 120.9038240, 2.06, Origin, 1.2);
    Bondi[51] = Atom("Tellurium",    "Te", 52.0, 129.9062290, 2.06, Origin, 1.2);
    Bondi[52] = Atom("Iodine",       "I",  53.0, 126.9044770, 1.98, Origin, 1.2);
    Bondi[53] = Atom("Xenon",        "Xe", 54.0, 131.9041480, 2.16, Origin, 1.2);
    Bondi[54] = Atom("Caesium",      "Cs", 55.0, 132.9054290, 3.43, Origin, 1.2);
    Bondi[55] = Atom("Barium",       "Ba", 56.0, 137.9052320, 2.68, Origin, 1.2);
    Bondi[56] = Atom("Lanthanum",    "La", 57.0, 138.9063470, 0.00, Origin, 1.2);
    Bondi[57] = Atom("Cerium",       "Ce", 58.0, 139.9054330, 0.00, Origin, 1.2);
    Bondi[58] = Atom("Praseodimium", "Pr", 59.0, 140.9076470, 0.00, Origin, 1.2);
    Bondi[59] = Atom("Neodymium",    "Nd", 60.0, 141.9077190, 0.00, Origin, 1.2);
    Bondi[60] = Atom("Promethium",   "Pm", 61.0, 144.9127430, 0.00, Origin, 1.2);
    Bondi[61] = Atom("Samarium",     "Sm", 62.0, 151.9197280, 0.00, Origin, 1.2);
    Bondi[62] = Atom("Europium",     "Eu", 63.0, 152.9212250, 0.00, Origin, 1.2);
    Bondi[63] = Atom("Gadolinium",   "Gd", 64.0, 157.9240190, 0.00, Origin, 1.2);
    Bondi[64] = Atom("Terbium",      "Tb", 65.0, 158.9253420, 0.00, Origin, 1.2);
    Bondi[65] = Atom("Dysprosium",   "Dy", 66.0, 163.9291710, 0.00, Origin, 1.2);
    Bondi[66] = Atom("Holmium",      "Ho", 67.0, 164.9303190, 0.00, Origin, 1.2);
    Bondi[67] = Atom("Erbium",       "Er", 68.0, 165.9302900, 0.00, Origin, 1.2);
    Bondi[68] = Atom("Thulium",      "Tm", 69.0, 168.9342120, 0.00, Origin, 1.2);
    Bondi[69] = Atom("Ytterbium",    "Yb", 70.0, 173.9388590, 0.00, Origin, 1.2);
    Bondi[70] = Atom("Lutetium",     "Lu", 71.0, 174.9407700, 0.00, Origin, 1.2);
    Bondi[71] = Atom("Hafnium",      "Hf", 72.0, 179.9465457, 0.00, Origin, 1.2);
    Bondi[72] = Atom("Tantalum",     "Ta", 73.0, 180.9479920, 0.00, Origin, 1.2);
    Bondi[73] = Atom("Tungsten",     "W",  74.0, 183.9509280, 0.00, Origin, 1.2);
    Bondi[74] = Atom("Rhenium",      "Re", 75.0, 186.9557440, 0.00, Origin, 1.2);
    Bondi[75] = Atom("Osmium",       "Os", 76.0, 191.9614670, 0.00, Origin, 1.2);
    Bondi[76] = Atom("Iridium",      "Ir", 77.0, 192.9629170, 0.00, Origin, 1.2);
    Bondi[77] = Atom("Platinum",     "Pt", 78.0, 194.9647660, 1.75, Origin, 1.2);
    Bondi[78] = Atom("Gold",         "Au", 79.0, 196.9665430, 1.66, Origin, 1.2);
    Bondi[79] = Atom("Mercury",      "Hg", 80.0, 201.9706170, 1.55, Origin, 1.2);
    Bondi[80] = Atom("Thallium",     "Tl", 81.0, 204.9744010, 1.96, Origin, 1.2);
    Bondi[81] = Atom("Lead",         "Pb", 82.0, 207.9766270, 2.02, Origin, 1.2);
    Bondi[82] = Atom("Bismuth",      "Bi", 83.0, 208.9803740, 2.07, Origin, 1.2);
    Bondi[83] = Atom("Polonium",     "Po", 84.0, 208.9824040, 1.97, Origin, 1.2);
    Bondi[84] = Atom("Astatine",     "At", 85.0, 209.9871260, 2.02, Origin, 1.2);
    Bondi[85] = Atom("Radon",        "Rn", 86.0, 222.0175710, 2.20, Origin, 1.2);
    Bondi[86] = Atom("Francium",     "Fr", 87.0, 223.0197330, 3.48, Origin, 1.2);
    Bondi[87] = Atom("Radium",       "Ra", 88.0, 226.0254030, 2.83, Origin, 1.2);
    Bondi[88] = Atom("Actinium",     "Ac", 89.0, 227.0277500, 0.00, Origin, 1.2);
    Bondi[89] = Atom("Thorium",      "Th", 90.0, 232.0380508, 0.00, Origin, 1.2);
    Bondi[90] = Atom("Protactinium", "Pa", 91.0, 231.0358800, 0.00, Origin, 1.2);
    Bondi[91] = Atom("Uranium",      "U",  92.0, 238.0507847, 1.86, Origin, 1.2);

// ------------------------------------------------------------

    return Bondi;
}

std::vector<Atom> & initUFF()
{

    static std::vector<Atom> UFF(92);
    Eigen::Vector3d Origin = Eigen::Vector3d::Zero();

// ------------------------------------------------------------
//            Atom("Name",    "Symbol", Charge,     Mass,   Radius, Center, Scaling)
    UFF[0]  = Atom("Hydrogen",     "H",   1.0,   1.0078250, 1.4430, Origin, 1.2);
    UFF[1]  = Atom("Helium",       "He",  2.0,   4.0026030, 1.8100, Origin, 1.2);
    UFF[2]  = Atom("Lithium",      "Li",  3.0,   7.0160050, 1.2255, Origin, 1.2);
    UFF[3]  = Atom("Beryllium",    "Be",  4.0,   9.0121830, 1.3725, Origin, 1.2);
    UFF[4]  = Atom("Boron",        "B",   5.0,  11.0093050, 2.0415, Origin, 1.2);
    UFF[5]  = Atom("Carbon",       "C",   6.0,  12.0000000, 1.9255, Origin, 1.2);
    UFF[6]  = Atom("Nitrogen",     "N",   7.0,  14.0030740, 1.8300, Origin, 1.2);
    UFF[7]  = Atom("Oxygen",       "O",   8.0,  15.9949150, 1.7500, Origin, 1.2);
    UFF[8]  = Atom("Fluorine",     "F",   9.0,  18.9984030, 1.6820, Origin, 1.2);
    UFF[9]  = Atom("Neon",         "Ne", 10.0,  19.9924390, 1.6215, Origin, 1.2);
    UFF[10] = Atom("Sodium",       "Na", 11.0,  22.9897700, 1.4915, Origin, 1.2);
    UFF[11] = Atom("Magnesium",    "Mg", 12.0,  23.9850450, 1.5105, Origin, 1.2);
    UFF[12] = Atom("Aluminium",    "Al", 13.0,  26.9815410, 2.2495, Origin, 1.2);
    UFF[13] = Atom("Silicon",      "Si", 14.0,  27.9769280, 2.1475, Origin, 1.2);
    UFF[14] = Atom("Phosphorus",   "P",  15.0,  30.9737630, 2.0735, Origin, 1.2);
    UFF[15] = Atom("Sulphur",      "S",  16.0,  31.9720720, 2.0175, Origin, 1.2);
    UFF[16] = Atom("Chlorine",     "Cl", 17.0,  34.9688530, 1.9735, Origin, 1.2);
    UFF[17] = Atom("Argon",        "Ar", 18.0,  39.9623830, 1.9340, Origin, 1.2);
    UFF[18] = Atom("Potassium",    "K",  19.0,  38.9637080, 1.9060, Origin, 1.2);
    UFF[19] = Atom("Calcium",      "Ca", 20.0,  39.9625910, 1.6995, Origin, 1.2);
    UFF[20] = Atom("Scandium",     "Sc", 21.0,  44.9559140, 1.6475, Origin, 1.2);
    UFF[21] = Atom("Titanium",     "Ti", 22.0,  47.9479470, 1.5875, Origin, 1.2);
    UFF[22] = Atom("Vanadium",     "V",  23.0,  50.9439630, 1.5720, Origin, 1.2);
    UFF[23] = Atom("Chromium",     "Cr", 24.0,  51.9405100, 1.5115, Origin, 1.2);
    UFF[24] = Atom("Manganese",    "Mn", 25.0,  54.9380460, 1.4805, Origin, 1.2);
    UFF[25] = Atom("Iron",         "Fe", 26.0,  55.9349390, 1.4560, Origin, 1.2);
    UFF[26] = Atom("Cobalt",       "Co", 27.0,  58.9331980, 1.4360, Origin, 1.2);
    UFF[27] = Atom("Nickel",       "Ni", 28.0,  57.9353470, 1.4170, Origin, 1.2);
    UFF[28] = Atom("Copper",       "Cu", 29.0,  62.9295990, 1.7475, Origin, 1.2);
    UFF[29] = Atom("Zinc",         "Zn", 30.0,  63.9291450, 1.3815, Origin, 1.2);
    UFF[30] = Atom("Gallium",      "Ga", 31.0,  68.9255810, 2.1915, Origin, 1.2);
    UFF[31] = Atom("Germanium",    "Ge", 32.0,  73.9211790, 2.1400, Origin, 1.2);
    UFF[32] = Atom("Arsenic",      "As", 33.0,  74.9215960, 2.1150, Origin, 1.2);
    UFF[33] = Atom("Selenium",     "Se", 34.0,  79.9165210, 2.1025, Origin, 1.2);
    UFF[34] = Atom("Bromine",      "Br", 35.0,  78.9183360, 2.0945, Origin, 1.2);
    UFF[35] = Atom("Krypton",      "Kr", 36.0,  83.9115060, 2.0705, Origin, 1.2);
    UFF[36] = Atom("Rubidium",     "Rb", 37.0,  84.9118000, 2.0570, Origin, 1.2);
    UFF[37] = Atom("Strontium",    "Sr", 38.0,  87.9056250, 1.8205, Origin, 1.2);
    UFF[38] = Atom("Yttrium",      "Y",  39.0,  88.9058560, 1.6725, Origin, 1.2);
    UFF[39] = Atom("Zirconium",    "Zr", 40.0,  89.9047080, 1.5620, Origin, 1.2);
    UFF[40] = Atom("Niobium",      "Nb", 41.0,  92.9063780, 1.5825, Origin, 1.2);
    UFF[41] = Atom("Molybdenum",   "Mo", 42.0,  97.9054050, 1.5260, Origin, 1.2);
    UFF[42] = Atom("Technetium",   "Tc", 43.0,  96.0000000, 1.4990, Origin, 1.2);
    UFF[43] = Atom("Ruthenium",    "Ru", 44.0, 101.9043480, 0.0000, Origin, 1.2);
    UFF[44] = Atom("Rhodium",      "Rh", 45.0, 102.9055030, 0.0000, Origin, 1.2);
    UFF[45] = Atom("Palladium",    "Pd", 46.0, 105.9034750, 0.0000, Origin, 1.2);
    UFF[46] = Atom("Silver",       "Ag", 47.0, 106.9050950, 0.0000, Origin, 1.2);
    UFF[47] = Atom("Cadmium",      "Cd", 48.0, 113.9033610, 0.0000, Origin, 1.2);
    UFF[48] = Atom("Indium",       "In", 49.0, 114.9038750, 0.0000, Origin, 1.2);
    UFF[49] = Atom("Tin",          "Sn", 50.0, 119.9021990, 0.0000, Origin, 1.2);
    UFF[50] = Atom("Antimony",     "Sb", 51.0, 120.9038240, 0.0000, Origin, 1.2);
    UFF[51] = Atom("Tellurium",    "Te", 52.0, 129.9062290, 0.0000, Origin, 1.2);
    UFF[52] = Atom("Iodine",       "I",  53.0, 126.9044770, 2.2500, Origin, 1.2);
    UFF[53] = Atom("Xenon",        "Xe", 54.0, 131.9041480, 0.0000, Origin, 1.2);
    UFF[54] = Atom("Caesium",      "Cs", 55.0, 132.9054290, 0.0000, Origin, 1.2);
    UFF[55] = Atom("Barium",       "Ba", 56.0, 137.9052320, 0.0000, Origin, 1.2);
    UFF[56] = Atom("Lanthanum",    "La", 57.0, 138.9063470, 0.0000, Origin, 1.2);
    UFF[57] = Atom("Cerium",       "Ce", 58.0, 139.9054330, 0.0000, Origin, 1.2);
    UFF[58] = Atom("Praseodimium", "Pr", 59.0, 140.9076470, 0.0000, Origin, 1.2);
    UFF[59] = Atom("Neodymium",    "Nd", 60.0, 141.9077190, 0.0000, Origin, 1.2);
    UFF[60] = Atom("Promethium",   "Pm", 61.0, 144.9127430, 0.0000, Origin, 1.2);
    UFF[61] = Atom("Samarium",     "Sm", 62.0, 151.9197280, 0.0000, Origin, 1.2);
    UFF[62] = Atom("Europium",     "Eu", 63.0, 152.9212250, 0.0000, Origin, 1.2);
    UFF[63] = Atom("Gadolinium",   "Gd", 64.0, 157.9240190, 0.0000, Origin, 1.2);
    UFF[64] = Atom("Terbium",      "Tb", 65.0, 158.9253420, 0.0000, Origin, 1.2);
    UFF[65] = Atom("Dysprosium",   "Dy", 66.0, 163.9291710, 0.0000, Origin, 1.2);
    UFF[66] = Atom("Holmium",      "Ho", 67.0, 164.9303190, 0.0000, Origin, 1.2);
    UFF[67] = Atom("Erbium",       "Er", 68.0, 165.9302900, 0.0000, Origin, 1.2);
    UFF[68] = Atom("Thulium",      "Tm", 69.0, 168.9342120, 0.0000, Origin, 1.2);
    UFF[69] = Atom("Ytterbium",    "Yb", 70.0, 173.9388590, 0.0000, Origin, 1.2);
    UFF[70] = Atom("Lutetium",     "Lu", 71.0, 174.9407700, 0.0000, Origin, 1.2);
    UFF[71] = Atom("Hafnium",      "Hf", 72.0, 179.9465457, 0.0000, Origin, 1.2);
    UFF[72] = Atom("Tantalum",     "Ta", 73.0, 180.9479920, 0.0000, Origin, 1.2);
    UFF[73] = Atom("Tungsten",     "W",  74.0, 183.9509280, 0.0000, Origin, 1.2);
    UFF[74] = Atom("Rhenium",      "Re", 75.0, 186.9557440, 0.0000, Origin, 1.2);
    UFF[75] = Atom("Osmium",       "Os", 76.0, 191.9614670, 0.0000, Origin, 1.2);
    UFF[76] = Atom("Iridium",      "Ir", 77.0, 192.9629170, 0.0000, Origin, 1.2);
    UFF[77] = Atom("Platinum",     "Pt", 78.0, 194.9647660, 0.0000, Origin, 1.2);
    UFF[78] = Atom("Gold",         "Au", 79.0, 196.9665430, 0.0000, Origin, 1.2);
    UFF[79] = Atom("Mercury",      "Hg", 80.0, 201.9706170, 0.0000, Origin, 1.2);
    UFF[80] = Atom("Thallium",     "Tl", 81.0, 204.9744010, 0.0000, Origin, 1.2);
    UFF[81] = Atom("Lead",         "Pb", 82.0, 207.9766270, 0.0000, Origin, 1.2);
    UFF[82] = Atom("Bismuth",      "Bi", 83.0, 208.9803740, 0.0000, Origin, 1.2);
    UFF[83] = Atom("Polonium",     "Po", 84.0, 208.9824040, 0.0000, Origin, 1.2);
    UFF[84] = Atom("Astatine",     "At", 85.0, 209.9871260, 0.0000, Origin, 1.2);
    UFF[85] = Atom("Radon",        "Rn", 86.0, 222.0175710, 0.0000, Origin, 1.2);
    UFF[86] = Atom("Francium",     "Fr", 87.0, 223.0197330, 0.0000, Origin, 1.2);
    UFF[87] = Atom("Radium",       "Ra", 88.0, 226.0254030, 0.0000, Origin, 1.2);
    UFF[88] = Atom("Actinium",     "Ac", 89.0, 227.0277500, 0.0000, Origin, 1.2);
    UFF[89] = Atom("Thorium",      "Th", 90.0, 232.0380508, 0.0000, Origin, 1.2);
    UFF[90] = Atom("Protactinium", "Pa", 91.0, 231.0358800, 0.0000, Origin, 1.2);
    UFF[91] = Atom("Uranium",      "U",  92.0, 238.0507847, 0.0000, Origin, 1.2);
// ------------------------------------------------------------

    return UFF;
}
