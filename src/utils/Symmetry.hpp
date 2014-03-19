#ifndef SYMMETRY_HPP
#define SYMMETRY_HPP

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

/*! \file Symmetry.hpp
 *  \class Symmetry
 *  \brief Contains very basic info about symmetry (only Abelian groups)
 *  \author Roberto Di Remigio
 *  \date 2014
 *
 *  Just a wrapper around a vector containing the generators of the group
 */

class Symmetry
{
	private:
		/*!
		 * Number of generators
		 */
		int nrGenerators_;
		/*!
		 * Generators
		 */
		int generators_[3];
		/*!
		 * Number of irreps
		 */
		int nrIrrep_;
	public:
		Symmetry() {}
		Symmetry(int nr_gen, int gen[3]) : nrGenerators_(nr_gen) 
			{
				// Transfer the passed generators array into generators_
				std::copy(gen, gen + nrGenerators_, generators_);
				// We can now initialize the number of irreps
				nrIrrep_ = int(std::pow(2, nrGenerators_));
			}
		Symmetry(const Symmetry & other) : 
			nrGenerators_(other.nrGenerators_), nrIrrep_(other.nrIrrep_)
			{
				std::copy(other.generators_, other.generators_ + nrGenerators_, generators_);
			}
		~Symmetry() {}
		int nrGenerators() const { return nrGenerators_; }
		int generators(int i) const { return generators_[i]; }
		int nrIrrep() const { return nrIrrep_; }
		static double parity(int i);
};

Symmetry buildGroup(int _nr_gen, int _gen1, int _gen2, int _gen3);

#endif // SYMMETRY_HPP
