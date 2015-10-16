/* pcmsolver_copyright_start */
/*
 *     PCMSolver, an API for the Polarizable Continuum Model
 *     Copyright (C) 2013 Roberto Di Remigio, Luca Frediani and contributors
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
 *     PCMSolver API, see: <https://repo.ctcc.no/projects/pcmsolver>
 */
/* pcmsolver_copyright_end */

#include "TsLessCavity.hpp"

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

#include "Config.hpp"
#include "FCMangle.hpp"

#include <Eigen/Core>

#include "Exception.hpp"
#include "Sphere.hpp"
#include "Symmetry.hpp"
#include "TimerInterface.hpp"

/*! \brief Fortran interface function to TsLess cavity generation
 *  \param[in] maxts maximum number of tesserae allowed
 *  \param[in] maxsph maximum number of spheres allowed
 *  \param[in] maxvert maximum number of vertices allowed
 *  \param[out] nesfp number of spheres (original + added)
 *  \param[out] nts number of generated tesserae
 *  \param[out] ntsirr number of generated irreducible tesserae
 *  \param[out] addsph number of added spheres
 *  \param[out] xtscor x-coordinate of tesserae centers (dimension maxts)
 *  \param[out] ytscor y-coordinate of tesserae centers (dimension maxts)
 *  \param[out] ztscor z-coordinate of tesserae centers (dimension maxts)
 *  \param[out] ar area of the tessera (dimension maxts)
 *  \param[out] xsphcor x-coordinate of the sphere center the tessera belongs to (dimension maxts)
 *  \param[out] ysphcor y-coordinate of the sphere center the tessera belongs to (dimension maxts)
 *  \param[out] zsphcor z-coordinate of the sphere center the tessera belongs to (dimension maxts)
 *  \param[out] rsph radii of the sphere the tessera belongs to, i.e. its curvature (dimension maxts)
 *  \param[out] xe x-coordinate of the sphere center (dimension nSpheres_ + maxAddedSpheres)
 *  \param[out] ye y-coordinate of the sphere center (dimension nSpheres_ + maxAddedSpheres)
 *  \param[out] ze z-coordinate of the sphere center (dimension nSpheres_ + maxAddedSpheres)
 *  \param[out] rin radius of the spheres (dimension nSpheres_ + maxAddedSpheres)
 *  \param[in] masses atomic masses (for inertia tensor formation in TSLESS)
 *  \param[in] nr_gen number of symmetry generators
 *  \param[in] gen1 first generator
 *  \param[in] gen2 second generator
 *  \param[in] gen3 third generator
 *  \param[in] avgArea average tesserae area
 *  \param[in] dmin mininal distance between sampling points
 *  \param[in] nord maximum order of continuous derivative of weight function
 *  \param[in] ifun whether to use the normalized or unnormalized form of the weight function
 *  \param[in] rsolv solvent probe radius
 *  \param[in] work scratch space
 */
#define tsless_driver \
    FortranCInterface_MODULE_(tsless_cavity, tsless_driver, TSLESS_CAVITY, TSLESS_DRIVER)
extern "C" void tsless_driver(size_t * maxts, size_t * maxsph, size_t * maxvert,
        int * nesfp, int * nts, int * ntsirr, int * addsph,
        double * xtscor, double * ytscor, double * ztscor, double * ar,
        double * xsphcor, double * ysphcor, double * zsphcor, double * rsph,
        double * xe, double * ye, double * ze, double * rin, double * masses,
        int * nr_gen, int * gen1, int * gen2, int * gen3,
        double * avgArea, double * dmin, int * nord, int * ifun, double * rsolv,
        double * work);

void TsLessCavity::build(size_t maxts, size_t maxsph, size_t maxvert)
{
    // This is a wrapper for the generatecavity_cpp_ function defined in the Fortran code TsLess.
    // Here we allocate the necessary arrays to be passed to TsLess, in particular we allow
    // for the insertion of additional spheres as in the most general formulation of the
    // GePol algorithm.

    int lwork = maxts*maxsph;
    double * xtscor  = new double[maxts];
    double * ytscor  = new double[maxts];
    double * ztscor  = new double[maxts];
    double * ar      = new double[maxts];
    double * xsphcor = new double[maxts];
    double * ysphcor = new double[maxts];
    double * zsphcor = new double[maxts];
    double * rsph    = new double[maxts];
    double * work    = new double[lwork];

    // Clean-up possible heap-crap
    std::fill_n(xtscor, maxts, 0.0);
    std::fill_n(ytscor, maxts, 0.0);
    std::fill_n(ztscor, maxts, 0.0);
    std::fill_n(ar, maxts, 0.0);
    std::fill_n(xsphcor, maxts, 0.0);
    std::fill_n(ysphcor, maxts, 0.0);
    std::fill_n(zsphcor, maxts, 0.0);
    std::fill_n(rsph, maxts, 0.0);
    std::fill_n(work, lwork, 0.0);

    int nts = 0;
    int ntsirr = 0;

    // If there's an overflow in the number of spheres TsLess will die.
    // The maximum number of spheres in PEDRA is set to 200 (primitive+additional)
    // so the integer here declared is just to have enough space C++ side to pass everything back.
    int maxAddedSpheres = 200;

    // Allocate vectors of size equal to nSpheres_ + maxAddedSpheres where maxAddedSpheres is the
    // maximum number of spheres we allow the algorithm to add to our original set.
    // If this number is exceeded, then the algorithm crashes (should look into this...)
    // After the cavity is generated we will update ALL the class data members, both related
    // to spheres and finite elements so that the cavity is fully formed.

    Eigen::VectorXd xv = Eigen::VectorXd::Zero(nSpheres_ + maxAddedSpheres);
    Eigen::VectorXd yv = Eigen::VectorXd::Zero(nSpheres_ + maxAddedSpheres);
    Eigen::VectorXd zv = Eigen::VectorXd::Zero(nSpheres_ + maxAddedSpheres);
    Eigen::VectorXd radii_scratch = Eigen::VectorXd::Zero(nSpheres_ +
                                    maxAddedSpheres); // Not to be confused with the data member inherited from Cavity!!!

    for ( int i = 0; i < nSpheres_; ++i ) {
        for ( int j = 0; j < 3; ++j ) {
            xv(i) = sphereCenter_(0, i);
            yv(i) = sphereCenter_(1, i);
            zv(i) = sphereCenter_(2, i);
        }
        radii_scratch(i) = sphereRadius_(i);
    }

    double * xe = xv.data();
    double * ye = yv.data();
    double * ze = zv.data();

    double *rin = radii_scratch.data();
    double * mass = new double[molecule_.nAtoms()];
    for (size_t i = 0; i < molecule_.nAtoms(); ++i) {
	    mass[i] = molecule_.masses(i);
    }

    addedSpheres = 0;
    // Number of generators and generators of the point group
    int nr_gen = molecule_.pointGroup().nrGenerators();
    int gen1 = molecule_.pointGroup().generators(0);
    int gen2 = molecule_.pointGroup().generators(1);
    int gen3 = molecule_.pointGroup().generators(2);

    int weightFunction = 1;
    // Go TsLess, Go!
    TIMER_ON("TsLessCavity::tsless_driver");
    tsless_driver(&maxts, &maxsph, &maxvert, &nSpheres_, &nts, &ntsirr, &addedSpheres,
                  xtscor, ytscor, ztscor, ar,
                  xsphcor, ysphcor, zsphcor, rsph,
                  xe, ye, ze, rin, mass,
                  &nr_gen, &gen1, &gen2, &gen3,
                  &averageArea_, &minDistance_, &derOrder_, &weightFunction, &probeRadius_, work);
    TIMER_OFF("TsLessCavity::tsless_driver");

    // The "intensive" part of updating the spheres related class data members will be of course
    // executed iff addedSpheres != 0
    if (addedSpheres != 0) {
        // Save the number of original spheres
        int orig = nSpheres_;
        // Update the nSpheres
        nSpheres_ += addedSpheres;
        // Resize sphereRadius and sphereCenter...
        sphereRadius_.resize(nSpheres_);
        sphereCenter_.resize(Eigen::NoChange, nSpheres_);
        // Transfer radii and centers.
        // Eigen has no push_back function, so we need to traverse all the spheres...
        sphereRadius_ = radii_scratch.head(nSpheres_);
        for ( int i = 0; i < nSpheres_; ++i ) {
            sphereCenter_(0, i) = xv(i);
            sphereCenter_(1, i) = yv(i);
            sphereCenter_(2, i) = zv(i);
        }
        // Now grow the vector<Sphere> containing the list of spheres
        for ( int i = orig;  i < nSpheres_; ++i ) {
            spheres_.push_back(Sphere(sphereCenter_.col(i), sphereRadius_(i)));
        }
    }

    nElements_ = static_cast<int>(nts);
    nIrrElements_ = static_cast<int>(ntsirr);
    elementCenter_.resize(Eigen::NoChange, nElements_);
    elementSphereCenter_.resize(Eigen::NoChange, nElements_);
    elementNormal_.resize(Eigen::NoChange, nElements_);
    elementArea_.resize(nElements_);
    elementRadius_.resize(nElements_);
    for( int i = 0; i < nElements_; ++i ) {
        elementCenter_(0,i) = xtscor[i];
        elementCenter_(1,i) = ytscor[i];
        elementCenter_(2,i) = ztscor[i];
        elementArea_(i) = ar[i];
        elementSphereCenter_(0,i) = xsphcor[i];
        elementSphereCenter_(1,i) = ysphcor[i];
        elementSphereCenter_(2,i) = zsphcor[i];
        elementRadius_(i) = rsph[i];
    }

    elementNormal_ = elementCenter_ - elementSphereCenter_;
    for( int i = 0; i < nElements_; ++i) {
        elementNormal_.col(i) /= elementNormal_.col(i).norm();
    }

    // Fill elements_ vector
    for (int i = 0; i < nElements_; ++i) {
        bool irr = false;
        int nv = 1; // TsLess does not generate spherical polygons!!
        // TsLess puts the irreducible tesserae first (? Check with Cris!)
        if (i < nIrrElements_) irr = true;
        Sphere sph(elementSphereCenter_.col(i), elementRadius_(i));
        Eigen::Matrix3Xd vertices, arcs;
        vertices.resize(Eigen::NoChange, nv);
        arcs.resize(Eigen::NoChange, nv);
        // FIXME index of the sphere the element belongs to
        elements_.push_back(Element(nv, 0,
                    elementArea_(i),
                    elementCenter_.col(i),
                    elementNormal_.col(i),
                    irr, sph,
                    vertices, arcs));
    }

    delete[] xtscor;
    delete[] ytscor;
    delete[] ztscor;
    delete[] ar;
    delete[] xsphcor;
    delete[] ysphcor;
    delete[] zsphcor;
    delete[] rsph;
    delete[] work;
    delete[] mass;

    built = true;

}

std::ostream & TsLessCavity::printCavity(std::ostream & os)
{
    os << "Cavity type: TsLess" << std::endl;
    os << "Average point weight = " << averageArea_ << " AU^2" << std::endl;
    os << "Minimal distance between sampling points = " << minDistance_ << " AU" << std::endl;
    os << "Switch function is of class C^" << derOrder_ << std::endl;
    os << "Addition of extra spheres enabled" << std::endl;
    os << "Probe radius = " << probeRadius_ << " AU" << std::endl;
    os << "Number of spheres = " << nSpheres_ << " [initial = " << nSpheres_ -
       addedSpheres << "; added = " << addedSpheres << "]" << std::endl;
    os << "Number of finite elements = " << nElements_;
    /*
    for (int i = 0; i < nElements_; i++)
    {
       os << std::endl;
       os << i+1 << " ";
       os << elementCenter_(0,i) << " ";
       os << elementCenter_(1,i) << " ";
       os << elementCenter_(2,i) << " ";
       os << elementArea_(i) << " ";
    }
    */
    return os;
}
