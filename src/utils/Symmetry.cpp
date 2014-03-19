#include "Symmetry.hpp"

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
 */

Symmetry buildGroup(int _nr_gen, int _gen1, int _gen2, int _gen3)
{
	int gen[3];
	gen[0] = _gen1;
	gen[1] = _gen2;
	gen[2] = _gen3;
	
	return Symmetry(_nr_gen, gen);
}

double Symmetry::parity(int i)
{
	static Eigen::VectorXd parity(8);
	parity << 1.0, -1.0, -1.0, 1.0, -1.0, 1.0, 1.0, -1.0;
	return parity(i);
}
