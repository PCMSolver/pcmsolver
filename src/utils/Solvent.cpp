#include "Solvent.hpp"

#include <map>
#include <ostream>
#include <string>

#include "Config.hpp"

Solvent::SolventMap & Solvent::initSolventMap() 
{
	static SolventMap availableSolvents;
  // ------------------------------------------------------------
        availableSolvents["Water"] = Solvent("Water", 78.39, 1.776, 1.385);
	availableSolvents["Methanol"] = Solvent("Methanol", 32.63, 1.758, 1.855);
	availableSolvents["Ethanol"] = Solvent("Ethanol", 24.55, 1.847, 2.18);
	availableSolvents["Chloroform"] = Solvent("Chloroform", 4.90, 2.085, 2.48);
	availableSolvents["Methylenechloride"] = Solvent("Methylenechloride", 8.93, 2.020, 2.27);
	availableSolvents["1,2-Dichloroethane"] = Solvent("1,2-Dichloroethane", 10.36, 2.085, 2.505);
	availableSolvents["Carbon tetrachloride"] = Solvent("Carbon tetrachloride", 2.228, 2.129, 2.685);
	availableSolvents["Benzene"] = Solvent("Benzene", 2.247, 2.244, 2.630);
	availableSolvents["Toluene"] = Solvent("Toluene", 2.379, 2.232, 2.82);
	availableSolvents["Chlorobenzene"] = Solvent("Chlorobenzene", 5.621, 2.320, 2.805);
	availableSolvents["Nitromethane"] = Solvent("Nitromethane", 38.20, 1.904, 2.155);
	availableSolvents["N-heptane"] = Solvent("N-heptane", 1.92, 1.918, 3.125);
	availableSolvents["Cyclohexane"] = Solvent("Cyclohexane", 2.023, 2.028, 2.815);
	availableSolvents["Aniline"] = Solvent("Aniline", 6.89, 2.506, 2.80);
	availableSolvents["Acetone"] = Solvent("Acetone", 20.7, 1.841, 2.38);
	availableSolvents["Tetrahydrofurane"] = Solvent("Tetrahydrofurane", 7.58, 1.971, 2.9);
	availableSolvents["Dimethylsulfoxide"] = Solvent("Dimethylsulfoxide", 46.7, 2.179, 2.455);
	availableSolvents["Acetonitrile"] = Solvent("Acetonitrile", 36.64, 1.806, 2.155);	
	availableSolvents["Explicit"] = Solvent("Explicit", 0.0, 0.0, 0.0);	
  // ------------------------------------------------------------
	return availableSolvents;
}

std::ostream & Solvent::printSolvent(std::ostream & os) 
{
	os << "Solvent name:          " << name_ << std::endl;
	os << "Static  permittivity = " << epsStatic_ << std::endl;
	os << "Optical permittivity = " << epsOptical_ << std::endl;
	os << "Solvent radius =       " << probeRadius_;
	return os;
}
