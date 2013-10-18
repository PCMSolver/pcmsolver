#include "Atom.hpp"

#include <vector>

#include "Config.hpp"

#include <Eigen/Dense>

std::vector<Atom> & Atom::initBondi() {
	
	static std::vector<Atom> Bondi(92);
	Eigen::Vector3d Origin = Eigen::Vector3d::Zero();

// ------------------------------------------------------------
	Bondi[0] = Atom("Hydrogen", "H", 1.0, 1.20, Origin, 1.2);
	Bondi[1] = Atom("Helium", "He", 2.0, 1.40, Origin, 1.2);
	Bondi[2] = Atom("Lithium", "Li", 3.0, 1.82, Origin, 1.2);
	Bondi[3] = Atom("Beryllium", "Be", 4.0, 1.53, Origin, 1.2);
	Bondi[4] = Atom("Boron", "B", 5.0, 1.92, Origin, 1.2);
	Bondi[5] = Atom("Carbon", "C", 6.0, 1.70, Origin, 1.2);
	Bondi[6] = Atom("Nitrogen", "N", 7.0, 1.55, Origin, 1.2);
	Bondi[7] = Atom("Oxygen", "O", 8.0, 1.52, Origin, 1.2);
	Bondi[8] = Atom("Fluorine", "F", 9.0, 1.47, Origin, 1.2);
	Bondi[9] = Atom("Neon", "Ne", 10.0, 1.54, Origin, 1.2);
	Bondi[10] = Atom("Sodium", "Na", 11.0, 2.27, Origin, 1.2);
	Bondi[11] = Atom("Magnesium", "Mg", 12.0, 1.73, Origin, 1.2);
	Bondi[12] = Atom("Aluminium", "Al", 13.0, 1.84, Origin, 1.2);
	Bondi[13] = Atom("Silicon", "Si", 14.0, 2.10, Origin, 1.2);
	Bondi[14] = Atom("Phosphorus", "P", 15.0, 1.80, Origin, 1.2);
	Bondi[15] = Atom("Sulphur", "S", 16.0, 1.80, Origin, 1.2);
	Bondi[16] = Atom("Chlorine", "Cl", 17.0, 1.75, Origin, 1.2);
	Bondi[17] = Atom("Argon", "Ar", 18.0, 1.88, Origin, 1.2);
	Bondi[18] = Atom("Potassium", "K", 19.0, 2.75, Origin, 1.2);
	Bondi[19] = Atom("Calcium", "Ca", 20.0, 2.31, Origin, 1.2);
	Bondi[20] = Atom("Scandium", "Sc", 21.0, 0.0, Origin, 1.2);
	Bondi[21] = Atom("Titanium", "Ti", 22.0, 0.0, Origin, 1.2);
	Bondi[22] = Atom("Vanadium", "V",  23.0, 0.0, Origin, 1.2);
	Bondi[23] = Atom("Chromium", "Cr", 24.0, 0.0, Origin, 1.2);
	Bondi[24] = Atom("Manganese", "Mn", 25.0, 0.0, Origin, 1.2);
	Bondi[25] = Atom("Iron", "Fe", 26.0, 0.0, Origin, 1.2);
	Bondi[26] = Atom("Cobalt", "Co", 27.0, 0.0, Origin, 1.2);
	Bondi[27] = Atom("Nickel", "Ni", 28.0, 1.63, Origin, 1.2);
	Bondi[28] = Atom("Copper", "Cu", 29.0, 1.4, Origin, 1.2);
	Bondi[29] = Atom("Zinc", "Zn", 30.0, 1.39, Origin, 1.2);
	Bondi[30] = Atom("Gallium", "Ga", 31.0, 1.87, Origin, 1.2);
	Bondi[31] = Atom("Germanium", "Ge", 32.0, 2.11, Origin, 1.2);
	Bondi[32] = Atom("Arsenic", "As", 33.0, 1.85, Origin, 1.2);
	Bondi[33] = Atom("Selenium", "Se", 34.0, 1.90, Origin, 1.2);
	Bondi[34] = Atom("Bromine", "Br", 35.0, 1.85, Origin, 1.2);
	Bondi[35] = Atom("Krypton", "Kr", 36.0, 2.02, Origin, 1.2);
	Bondi[36] = Atom("Rubidium", "Rb", 37.0, 3.03, Origin, 1.2);
	Bondi[37] = Atom("Strontium", "Sr", 38.0, 2.49, Origin, 1.2);
	Bondi[38] = Atom("Yttrium", "Y", 39.0, 0.0, Origin, 1.2);
	Bondi[39] = Atom("Zirconium", "Zr", 40.0, 0.0, Origin, 1.2);
	Bondi[40] = Atom("Niobium", "Nb", 41.0, 0.0, Origin, 1.2);
	Bondi[41] = Atom("Molybdenum", "Mo", 42.0, 0.0, Origin, 1.2);
	Bondi[42] = Atom("Technetium", "Tc", 43.0, 0.0, Origin, 1.2);
	Bondi[43] = Atom("Ruthenium", "Ru", 44.0, 0.0, Origin, 1.2);
	Bondi[44] = Atom("Rhodium", "Rh", 45.0, 0.0, Origin, 1.2);
	Bondi[45] = Atom("Palladium", "Pd", 46.0, 1.63, Origin, 1.2);
	Bondi[46] = Atom("Silver", "Ag", 47.0, 1.72, Origin, 1.2);
	Bondi[47] = Atom("Cadmium", "Cd", 48.0, 1.58, Origin, 1.2);
	Bondi[48] = Atom("Indium", "In", 49.0, 1.93, Origin, 1.2);
	Bondi[49] = Atom("Tin", "Sn", 50.0, 2.17, Origin, 1.2);
	Bondi[50] = Atom("Antimony", "Sb", 51.0, 2.06, Origin, 1.2);
	Bondi[51] = Atom("Tellurium", "Te", 52.0, 2.06, Origin, 1.2);
	Bondi[52] = Atom("Iodine", "I",  53.0, 1.98, Origin, 1.2);
	Bondi[53] = Atom("Xenon", "Xe", 54.0, 2.16, Origin, 1.2);
	Bondi[54] = Atom("Cesium", "Cs", 55.0, 3.43, Origin, 1.2);
	Bondi[55] = Atom("Barium", "Ba", 56.0, 2.68, Origin, 1.2);
	Bondi[56] = Atom("Lanthanum", "La", 57.0, 0.0, Origin, 1.2);
	Bondi[57] = Atom("Cerium", "Ce", 58.0, 0.0, Origin, 1.2);
	Bondi[58] = Atom("Praseodimium", "Pr", 59.0, 0.0, Origin, 1.2);
	Bondi[59] = Atom("Neodymium", "Nd", 60.0, 0.0, Origin, 1.2);
	Bondi[60] = Atom("Promethium", "Pm", 61.0, 0.0, Origin, 1.2);
	Bondi[61] = Atom("Samarium", "Sm", 62.0, 0.0, Origin, 1.2);
	Bondi[62] = Atom("Europium", "Eu", 63.0, 0.0, Origin, 1.2);
	Bondi[63] = Atom("Gadolinium", "Gd", 64.0, 0.0, Origin, 1.2);
	Bondi[64] = Atom("Terbium", "Tb", 65.0, 0.0, Origin, 1.2);
	Bondi[65] = Atom("Dysprosium", "Dy", 66.0, 0.0, Origin, 1.2);
	Bondi[66] = Atom("Holmium", "Ho", 67.0, 0.0, Origin, 1.2);
	Bondi[67] = Atom("Erbium", "Er", 68.0, 0.0, Origin, 1.2);
	Bondi[68] = Atom("Thulium", "Tm", 69.0, 0.0, Origin, 1.2);
	Bondi[69] = Atom("Ytterbium", "Yb", 70.0, 0.0, Origin, 1.2);
	Bondi[70] = Atom("Lutetium", "Lu", 71.0, 0.0, Origin, 1.2);
	Bondi[71] = Atom("Hafnium", "Hf", 72.0, 0.0, Origin, 1.2);
	Bondi[72] = Atom("Tantalum", "Ta", 73.0, 0.0, Origin, 1.2);
	Bondi[73] = Atom("Tungsten", "W", 74.0, 0.0, Origin, 1.2);
	Bondi[74] = Atom("Rhenium", "Re", 75.0, 0.0, Origin, 1.2);
	Bondi[75] = Atom("Osmium", "Os", 76.0, 0.0, Origin, 1.2);
	Bondi[76] = Atom("Iridium", "Ir", 77.0, 0.0, Origin, 1.2);
	Bondi[77] = Atom("Platinum", "Pt", 78.0, 1.75, Origin, 1.2);
	Bondi[78] = Atom("Gold", "Au", 79.0, 1.66, Origin, 1.2);
	Bondi[79] = Atom("Mercury", "Hg",  80.0, 1.55, Origin, 1.2);
	Bondi[80] = Atom("Tallium", "Tl", 81.0, 1.96, Origin, 1.2);
        Bondi[81] = Atom("Lead", "Pb", 82.0, 2.02, Origin, 1.2);
	Bondi[82] = Atom("Bismuth", "Bi", 83.0, 2.07, Origin, 1.2);
	Bondi[83] = Atom("Polonium", "Po", 84.0, 1.97, Origin, 1.2);
	Bondi[84] = Atom("Astatine", "At", 85.0, 2.02, Origin, 1.2);
	Bondi[85] = Atom("Radon", "Rn", 86.0, 2.20, Origin, 1.2);
	Bondi[86] = Atom("Francium", "Fr",  87.0, 3.48, Origin, 1.2);
	Bondi[87] = Atom("Radium", "Ra", 88.0, 2.83, Origin, 1.2);
	Bondi[88] = Atom("Actinium", "Ac", 89.0, 0.00, Origin, 1.2);
	Bondi[89] = Atom("Thorium", "Th", 90.0, 0.00, Origin, 1.2);
	Bondi[90] = Atom("Protactinium", "Pa", 91.0, 0.00, Origin, 1.2);
	Bondi[91] = Atom("Uranium", "U", 92.0, 1.86, Origin, 1.2);

// ------------------------------------------------------------
	
	return Bondi;
}

std::vector<Atom> & Atom::initUFF() {
	
	static std::vector<Atom> UFF(89);
	Eigen::Vector3d Origin = Eigen::Vector3d::Zero();

// ------------------------------------------------------------
	UFF[0] = Atom("Hydrogen", "H", 1.0, 1.4430, Origin, 1.0);
	UFF[1] = Atom("Helium", "He", 2.0, 1.810, Origin, 1.0);
	UFF[2] = Atom("Lithium", "Li", 3.0, 1.2255, Origin, 1.0);
	UFF[3] = Atom("Beryllium", "Be", 4.0, 1.3725, Origin, 1.0);
	UFF[4] = Atom("Boron", "B", 5.0, 2.0415, Origin, 1.0);
	UFF[5] = Atom("Carbon", "C", 6.0, 1.9255, Origin, 1.0);
	UFF[6] = Atom("Nitrogen", "N", 7.0, 1.83, Origin, 1.0);
	UFF[7] = Atom("Oxygen", "O", 8.0, 1.75, Origin, 1.0);
	UFF[8] = Atom("Fluorine", "F", 9.0, 1.6820, Origin, 1.0);
	UFF[9] = Atom("Neon", "Ne", 10.0, 1.6215, Origin, 1.0);
	UFF[10] = Atom("Sodium", "Na", 11.0, 1.4915, Origin, 1.0);
	UFF[11] = Atom("Magnesium", "Mg", 12.0, 1.5105, Origin, 1.0);
	UFF[12] = Atom("Aluminium", "Al", 13.0, 2.2495, Origin, 1.0);
	UFF[13] = Atom("Silicon", "Si", 14.0, 2.1475, Origin, 1.0);
	UFF[14] = Atom("Phosphorus", "P", 15.0, 2.0735, Origin, 1.0);
	UFF[15] = Atom("Sulphur", "S", 16.0, 2.0175, Origin, 1.0);
	UFF[16] = Atom("Chlorine", "Cl", 17.0, 1.9735, Origin, 1.0);
	UFF[17] = Atom("Argon", "Ar", 18.0, 1.9340, Origin, 1.0);
	UFF[18] = Atom("Potassium", "K", 19.0, 1.9060, Origin, 1.0);
	UFF[19] = Atom("Calcium", "Ca", 20.0, 1.6995, Origin, 1.0);
	UFF[20] = Atom("Scandium", "Sc", 21.0, 1.6475, Origin, 1.0);
	UFF[21] = Atom("Titanium", "Ti", 22.0, 1.5875, Origin, 1.0);
	UFF[22] = Atom("Vanadium", "V",  23.0, 1.5720, Origin, 1.0);
	UFF[23] = Atom("Chromium", "Cr", 24.0, 1.5115, Origin, 1.0);
	UFF[24] = Atom("Manganese", "Mn", 25.0, 1.4805, Origin, 1.0);
	UFF[25] = Atom("Iron", "Fe", 26.0, 1.4560, Origin, 1.0);
	UFF[26] = Atom("Cobalt", "Co", 27.0, 1.4360, Origin, 1.0);
	UFF[27] = Atom("Nickel", "Ni", 28.0, 1.4170, Origin, 1.0);
	UFF[28] = Atom("Copper", "Cu", 29.0, 1.7475, Origin, 1.0);
	UFF[29] = Atom("Zinc", "Zn", 30.0, 1.3815, Origin, 1.0);
	UFF[30] = Atom("Gallium", "Ga", 31.0, 2.1915, Origin, 1.0);
	UFF[31] = Atom("Germanium", "Ge", 32.0, 2.1400, Origin, 1.0);
	UFF[32] = Atom("Arsenic", "As", 33.0, 2.1150, Origin, 1.0);
	UFF[33] = Atom("Selenium", "Se", 34.0, 2.1025, Origin, 1.0);
	UFF[34] = Atom("Bromine", "Br", 35.0, 2.0945, Origin, 1.0);
	UFF[35] = Atom("Krypton", "Kr", 36.0, 2.0705, Origin, 1.0);
	UFF[36] = Atom("Rubidium", "Rb", 37.0, 2.0570, Origin, 1.0);
	UFF[37] = Atom("Strontium", "Sr", 38.0, 1.8205, Origin, 1.0);
	UFF[38] = Atom("Yttrium", "Y", 39.0, 1.6725, Origin, 1.0);
	UFF[39] = Atom("Zirconium", "Zr", 40.0, 1.5620, Origin, 1.0);
	UFF[40] = Atom("Niobium", "Nb", 41.0, 1.5825, Origin, 1.0);
	UFF[41] = Atom("Molybdenum", "Mo", 42.0, 1.5260, Origin, 1.0);
	UFF[42] = Atom("Technetium", "Tc", 43.0, 1.4990, Origin, 1.0);
	UFF[43] = Atom("Ruthenium", "Ru", 44.0, 0.0, Origin, 1.0);
	UFF[44] = Atom("Rhodium", "Rh", 45.0, 0.0, Origin, 1.0);
	UFF[45] = Atom("Palladium", "Pd", 46.0, 0.0, Origin, 1.0);
	UFF[46] = Atom("Silver", "Ag", 47.0, 0.0, Origin, 1.0);
	UFF[47] = Atom("Cadmium", "Cd", 48.0, 0.0, Origin, 1.0);
	UFF[48] = Atom("Indium", "In", 49.0, 0.0, Origin, 1.0);
	UFF[49] = Atom("Tin", "Sn", 50.0, 0.0, Origin, 1.0);
	UFF[50] = Atom("Antimony", "Sb", 51.0, 0.0, Origin, 1.0);
	UFF[51] = Atom("Tellurium", "Te", 52.0, 0.0, Origin, 1.0);
	UFF[52] = Atom("Iodine", "I",  53.0, 2.2500, Origin, 1.0);
	UFF[53] = Atom("Xenon", "Xe", 54.0, 0.0, Origin, 1.0);
	UFF[54] = Atom("Cesium", "Cs", 55.0, 0.0, Origin, 1.0);
	UFF[55] = Atom("Barium", "Ba", 56.0, 0.0, Origin, 1.0);
	UFF[56] = Atom("Lanthanum", "La", 57.0, 0.0, Origin, 1.0);
	UFF[57] = Atom("Cerium", "Ce", 58.0, 0.0, Origin, 1.0);
	UFF[58] = Atom("Praseodimium", "Pr", 59.0, 0.0, Origin, 1.0);
	UFF[59] = Atom("Neodymium", "Nd", 60.0, 0.0, Origin, 1.0);
	UFF[60] = Atom("Promethium", "Pm", 61.0, 0.0, Origin, 1.0);
	UFF[61] = Atom("Samarium", "Sm", 62.0, 0.0, Origin, 1.0);
	UFF[62] = Atom("Europium", "Eu", 63.0, 0.0, Origin, 1.0);
	UFF[63] = Atom("Gadolinium", "Gd", 64.0, 0.0, Origin, 1.0);
	UFF[64] = Atom("Terbium", "Tb", 65.0, 0.0, Origin, 1.0);
	UFF[65] = Atom("Dysprosium", "Dy", 66.0, 0.0, Origin, 1.0);
	UFF[66] = Atom("Holmium", "Ho", 67.0, 0.0, Origin, 1.0);
	UFF[67] = Atom("Erbium", "Er", 68.0, 0.0, Origin, 1.0);
	UFF[68] = Atom("Thulium", "Tm", 69.0, 0.0, Origin, 1.0);
	UFF[69] = Atom("Ytterbium", "Yb", 70.0, 0.0, Origin, 1.0);
	UFF[70] = Atom("Lutetium", "Lu", 71.0, 0.0, Origin, 1.0);
	UFF[71] = Atom("Hafnium", "Hf", 72.0, 0.0, Origin, 1.0);
	UFF[72] = Atom("Tantalum", "Ta", 73.0, 0.0, Origin, 1.0);
	UFF[73] = Atom("Tungsten", "W", 74.0, 0.0, Origin, 1.0);
	UFF[74] = Atom("Rhenium", "Re", 75.0, 0.0, Origin, 1.0);
	UFF[75] = Atom("Osmium", "Os", 76.0, 0.0, Origin, 1.0);
	UFF[76] = Atom("Iridium", "Ir", 77.0, 0.0, Origin, 1.0);
	UFF[77] = Atom("Platinum", "Pt", 78.0, 0.0, Origin, 1.0);
	UFF[78] = Atom("Gold", "Au", 79.0, 0.0, Origin, 1.0);
	UFF[79] = Atom("Mercury", "Hg",  80.0, 0.0, Origin, 1.0);
	UFF[80] = Atom("Tallium", "Tl", 81.0, 0.0, Origin, 1.0);
        UFF[81] = Atom("Lead", "Pb", 82.0, 0.0, Origin, 1.0);
	UFF[82] = Atom("Bismuth", "Bi", 83.0, 0.0, Origin, 1.0);
	UFF[83] = Atom("Polonium", "Po", 84.0, 0.0, Origin, 1.0);
	UFF[84] = Atom("Astatine", "At", 85.0, 0.0, Origin, 1.0);
	UFF[85] = Atom("Radon", "Rn", 86.0, 0.0, Origin, 1.0);
	UFF[86] = Atom("Francium", "Fr",  87.0, 0.0, Origin, 1.0);
	UFF[87] = Atom("Radium", "Ra", 88.0, 0.0, Origin, 1.0);
// ------------------------------------------------------------
	
	return UFF;
}

