#ifndef ATOM_HPP
#define ATOM_HPP

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

#include "PhysicalConstants.hpp"

/*! \file Atom.hpp
 *  \class Atom
 *  \brief Class describing an atom.
 *  \author Roberto Di Remigio
 *  \date 2011
 *
 *  This class contains all the radii sets available in the module.
 *  They can be obtained by the client through static functions initializing
 *  the sets via the lazy evaluation idiom.
 */

class Atom 
{
	public:
		Atom(){}
                Atom(const std::string & _element, const std::string & _symbol, double _charge,                              
                      double _radius, const Eigen::Vector3d & _coord, double _scaling = 1.0, const std::string & _colour = "Violet")
			: atomElement(_element), atomSymbol(_symbol), atomCharge(_charge), atomRadius(_radius),
			atomCoord(_coord), atomRadiusScaling(_scaling), atomColour(_colour) {}
                Atom(const std::string & _element, const std::string & _symbol, double _charge, double _radius)
			: atomElement(_element), atomSymbol(_symbol), atomCharge(_charge), atomRadius(_radius)
		{
			Eigen::Vector3d Origin(0.0, 0.0, 0.0);
			std::string colour = "Violet";
  			atomCoord = Origin;
  			atomColour = colour;
			atomRadiusScaling = 1.0;
		}
                ~Atom(){}
		std::string getAtomElement() { return atomElement; }
                void setAtomElement(const std::string & _element) { atomElement = _element; }
		std::string getAtomSymbol() { return atomSymbol; }
                void setAtomSymbol(const std::string & _symbol) { atomSymbol = _symbol; }
		Eigen::Vector3d & getAtomCoord() { return atomCoord; }
                void setAtomCoord(Eigen::Vector3d & _coord){ atomCoord = _coord; }
                double getAtomCharge() { return atomCharge; }
                void setAtomCharge( double _charge ){ atomCharge = _charge; }
                double getAtomRadius() { return (atomRadius / convertBohrToAngstrom); }
                void setAtomRadius( double _radius ) { atomRadius = _radius; }
                double getAtomRadiusScaling() { return atomRadiusScaling; }
                void setAtomRadiusScaling(double _scaling) { atomRadiusScaling = _scaling; }
		std::string getAtomColour() { return atomColour; }
                void setAtomColour(const std::string & _colour){ atomColour = _colour; }
                /*! \brief Returns a reference to a vector<Atom> containing Bondi van der Waals radii.                                                                        	
	         * 
		 * The van der Waals radii are taken from:
	         * --- A. Bondi, J. Phys. Chem. 68, 441-451 (1964) 
                 * complemented with the ones reported in:
                 * --- M. Mantina, A. C. Chamberlin, R. Valero, C. J. Cramer, D. G. Truhlar,
                 *     J. Phys. Chem. A, 113, 5806-5812 (2009)
                 * We are here using Angstrom as in the papers.
                 * The getAtomRadius method will perform the conversion Angstrom to AU.
	         */
		static std::vector<Atom> & initBondi();
                /*! \brief Returns a reference to a vector<Atom> containing UFF radii.
	         * 
		 * The UFF set of radii is taken from:
	         * --- A. Rappé, C. J. Casewit, K. S. Colwell, W. A. Goddard, W. M. Skiff,
                 *     J. Am. Chem. Soc., 114, 10024-10035 (1992)
                 * We are here using Angstrom as in the paper. 
                 * The getAtomRadius method will perform the conversion Angstrom to AU.
	         */
                static std::vector<Atom> & initUFF();
 	private:
		std::string atomElement;
		std::string atomSymbol;
                double atomCharge;
                double atomRadius;
		Eigen::Vector3d atomCoord;
                double atomRadiusScaling;
		std::string atomColour;
};

#endif // ATOM_HPP
