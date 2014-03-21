#ifndef ATOM_HPP
#define ATOM_HPP

#include <string>
#include <vector>

#include "Config.hpp"

#include "EigenPimpl.hpp"

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
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW // See http://eigen.tuxfamily.org/dox/group__TopicStructHavingEigenMembers.html
    Atom() {}
    Atom(const std::string & element, const std::string & symbol, double charge,
         double radius, const Eigen::Vector3d & coord, double scaling = 1.0,
         const std::string & colour = "Violet")
        : atomElement_(element), atomSymbol_(symbol), atomCharge_(charge),
          atomRadius_(radius),
          atomCoord_(coord), atomRadiusScaling_(scaling), atomColour_(colour) {}
    Atom(const std::string & element, const std::string & symbol, double charge,
         double radius)
        : atomElement_(element), atomSymbol_(symbol), atomCharge_(charge),
          atomRadius_(radius) {
        Eigen::Vector3d Origin(0.0, 0.0, 0.0);
        std::string colour = "Violet";
        atomCoord_ = Origin;
        atomColour_ = colour;
        atomRadiusScaling_ = 1.0;
    }
    ~Atom() {}
    std::string atomElement() const { return atomElement_; }
    void atomElement(const std::string & element) { atomElement_ = element; }
    std::string atomSymbol() const { return atomSymbol_; }
    void atomSymbol(const std::string & symbol) { atomSymbol_ = symbol; }
    const Eigen::Vector3d & atomCoord() const { return atomCoord_; }
    void atomCoord(Eigen::Vector3d & coord) { atomCoord_ = coord; }
    double atomCharge() const { return atomCharge_; }
    void atomCharge(double charge) { atomCharge_ = charge; }
    double atomRadius() const { return (atomRadius_ / convertBohrToAngstrom); }
    void atomRadius(double radius) { atomRadius_ = radius; }
    double atomRadiusScaling() const { return atomRadiusScaling_; }
    void atomRadiusScaling(double scaling) { atomRadiusScaling_ = scaling; }
    std::string atomColour() const { return atomColour_; }
    void atomColour(const std::string & colour) { atomColour_ = colour; }
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
    std::string atomElement_;
    std::string atomSymbol_;
    double atomCharge_;
    double atomRadius_;
    Eigen::Vector3d atomCoord_;
    double atomRadiusScaling_;
    std::string atomColour_;
};

#endif // ATOM_HPP
