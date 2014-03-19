#ifndef SYMMETRY_HPP
#define SYMMETRY_HPP

#include <cmath>
#include <string>
#include <vector>

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
		 * Integer representing the group
		 */
		int groupInteger_;
		/*!
		 * The Schoenflies symbol
		 */
		std::string groupName_;
		/*!
		 * Number of generators
		 */
		int nrGenerators_;
		/*!
		 * Generators
		 */
		int generators_[3];
	public:
		Symmetry() {}
		Symmetry(int nr_gen, int gen[3]) : 
			groupInteger_(0), groupName_(""), nrGenerators_(nr_gen) 
			{
				std::copy(gen, gen + nrGenerators_, generators_);
			};
		Symmetry(int pGroup, const std::string & name, int nr_gen, int gen[3]) : 
			groupInteger_(pGroup), groupName_(name), nrGenerators_(nr_gen)
			{
				std::copy(gen, gen + nrGenerators_, generators_);
			};
		Symmetry(const Symmetry & other) : 
			groupInteger_(other.groupInteger_), groupName_(other.groupName_), 
			nrGenerators_(other.nrGenerators_)
			{
				std::copy(other.generators_, other.generators_ + nrGenerators_, generators_);
			};
		~Symmetry() {}
		int groupInteger() const { return groupInteger_; }
		std::string groupName() const { return groupName_; }
		int nrGenerators() const { return nrGenerators_; }
		int generators(int i) const { return generators_[i]; }
		/*!
		 * Number of irreducible representations,
		 * nrIrrep_ = 2**nrGenerators_
		 * This is also the number of operations in the group.
		 */
		int nrIrrep() const { return int(std::pow(2, nrGenerators_)); }
		/*!
		 * Number of nontrivial symmetry operations,
		 * nontrivialOps_ = 2**nrGenerators_ - 1 
		 */
		int nontrivialOps() const { return (int(std::pow(2, nrGenerators_)) - 1); }
		static double parity(int i);
};

/*!
 * int-to-point_group mapping
 */
enum pGroup { C1, C2, Cs, Ci, D2, C2v, C2h, D2h };

Symmetry buildGroup(int _nr_gen, int _gen1, int _gen2, int _gen3);

#endif // SYMMETRY_HPP
