#include "Symmetry.hpp"

#include <stdexcept>

#include "Config.hpp"

// Disable obnoxious warnings from Eigen headers
#if defined (__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wall" 
#pragma GCC diagnostic ignored "-Weffc++" 
#pragma GCC diagnostic ignored "-Wextra"
#include <Eigen/Dense>
#pragma GCC diagnostic pop
#elif (__INTEL_COMPILER)
#pragma warning push
#pragma warning disable "-Wall"
#include <Eigen/Dense>
#pragma warning pop
#endif
/*
 * A. Indexing of symmetry operations and their mapping to a bitstring:
 *      zyx         Parity
 *   0  000    E      1.0
 *   1  001   Oyz    -1.0
 *   2  010   Oxz    -1.0
 *   3  011   C2z     1.0
 *   4  100   Oxy    -1.0
 *   5  101   C2y     1.0
 *   6  110   C2x     1.0
 *   7  111    i     -1.0
 * B. Indexing of irreps for the Abelian groups
 * C1:  A  <-> 0
 * Cs:  A' <-> 0; A'' <-> 1
 * C2:  A  <-> 0; B   <-> 1
 * Ci:  Ag <-> 0; Au  <-> 1
 * C2h: Ag <-> 0; Au  <-> 1; Bu  <-> 2; Bg  <-> 3
 * D2:  A  <-> 0; B3  <-> 1; B2  <-> 2; B1  <-> 3
 * C2v: A1 <-> 0; B1  <-> 1; B2  <-> 2; A2  <-> 3
 * D2h: Ag <-> 0; B3u <-> 1; B2u <-> 2; B1g <-> 3; B1u <-> 4; B2g <-> 5; B3g <-> 6; Au <-> 7
 *
 */

Symmetry buildGroup(int _pGroup)
{
	Symmetry group;
	switch (_pGroup)
	{
		case C1:
			group = Symmetry(C1, "C1", 0);
			break;
		case Cs:
			group = Symmetry(Cs, "Cs", 1);
			break;
		case C2:
			group = Symmetry(C2, "C2", 1);
			break;
		case Ci:
			group = Symmetry(Ci, "Ci", 1);
			break;
		case C2h:
			group = Symmetry(C2h, "C2h", 2);
			break;
		case D2:
			group = Symmetry(D2, "D2", 2);
			break;
		case C2v:
			group = Symmetry(C2v, "C2v", 2);
			break;
		case D2h:
			group = Symmetry(D2h, "D2h", 3);
			break;
		default:
			throw std::runtime_error("Point group not recognized.");
	}
	return group;
}

double Symmetry::parity(int i)
{
	static Eigen::VectorXd parity(8);
	parity << 1.0, -1.0, -1.0, 1.0, -1.0, 1.0, 1.0, -1.0;
	return parity(i);
}
