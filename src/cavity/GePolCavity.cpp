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
 *     PCMSolver API, see: <http://pcmsolver.github.io/pcmsolver-doc>
 */
/* pcmsolver_copyright_end */

#include "GePolCavity.hpp"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "Config.hpp"
#include "FCMangle.hpp"

#include <Eigen/Core>
#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>

#include "Sphere.hpp"
#include "Symmetry.hpp"

/*! \brief Interface to the Fortran PEDRA code
 *  \param[in] maxts maximum number of tesserae allowed
 *  \param[in] maxsph maximum number of spheres allowed
 *  \param[in] maxvert maximum number of vertices allowed
 *  \param[out] xtscor x-coordinate of tesserae centers (dimension maxts)
 *  \param[out] ytscor y-coordinate of tesserae centers (dimension maxts)
 *  \param[out] ztscor z-coordinate of tesserae centers (dimension maxts)
 *  \param[out] ar area of the tessera (dimension maxts)
 *  \param[out] xsphcor x-coordinate of the sphere center the tessera belongs to (dimension maxts)
 *  \param[out] ysphcor y-coordinate of the sphere center the tessera belongs to (dimension maxts)
 *  \param[out] zsphcor z-coordinate of the sphere center the tessera belongs to (dimension maxts)
 *  \param[out] rsph radii of the sphere the tessera belongs to, i.e. its curvature (dimension maxts)
 *  \param[out] nts number of generated tesserae
 *  \param[out] ntsirr number of generated irreducible tesserae
 *  \param[out] nesfp number of spheres (original + added)
 *  \param[out] addsph number of added spheres
 *  \param[out] xe x-coordinate of the sphere center (dimension nSpheres_ + maxAddedSpheres)
 *  \param[out] ye y-coordinate of the sphere center (dimension nSpheres_ + maxAddedSpheres)
 *  \param[out] ze z-coordinate of the sphere center (dimension nSpheres_ + maxAddedSpheres)
 *  \param[out] rin radius of the spheres (dimension nSpheres_ + maxAddedSpheres)
 *  \param[in] masses atomic masses (for inertia tensor formation in PEDRA)
 *  \param[in] avgArea average tesserae area
 *  \param[in] rsolv solvent probe radius
 *  \param[in] ret minimal radius for an added sphere
 *  \param[in] nr_gen number of symmetry generators
 *  \param[in] gen1 first generator
 *  \param[in] gen2 second generator
 *  \param[in] gen3 third generator
 *  \param[out] nvert number of vertices per tessera
 *  \param[out] vert coordinates of tesserae vertices
 *  \param[out] centr centers of arcs defining the edges of the tesserae
 */
#define pedra_driver\
    FortranCInterface_GLOBAL_(pedra_driver, PEDRA_DRIVER)
extern "C" void pedra_driver(size_t * maxts, size_t * maxsph, size_t * maxvert,
        double * xtscor, double * ytscor, double * ztscor, double * ar,
        double * xsphcor, double * ysphcor, double * zsphcor, double * rsph,
        int * nts, int * ntsirr, int * nesfp, int * addsph,
        double * xe, double * ye, double * ze, double * rin, double * masses,
        double * avgArea, double * rsolv, double * ret,
        int * nr_gen, int * gen1, int * gen2, int * gen3,
        int * nvert, double * vert, double * centr,
        int * isphe, char * pedra, int * len_f_pedra);


void GePolCavity::build(const std::string & suffix, size_t maxts, size_t maxsp, size_t maxvert)
{

    // This is a wrapper for the pedra_driver function defined in the Fortran code PEDRA.
    // Here we allocate the necessary arrays to be passed to PEDRA, in particular we allow
    // for the insertion of additional spheres as in the most general formulation of the
    // GePol algorithm.

    double * xtscor  = new double[maxts];
    double * ytscor  = new double[maxts];
    double * ztscor  = new double[maxts];
    double * ar      = new double[maxts];
    double * xsphcor = new double[maxts];
    double * ysphcor = new double[maxts];
    double * zsphcor = new double[maxts];
    double * rsph    = new double[maxts];
    int    * nvert   = new int[maxts];
    double * vert    = new double[30 * maxts];
    double * centr   = new double[30 * maxts];
    int    * isphe   = new int[maxts];

    // Clean-up possible heap-crap
    std::fill_n(xtscor, maxts, 0.0);
    std::fill_n(ytscor, maxts, 0.0);
    std::fill_n(ztscor, maxts, 0.0);
    std::fill_n(ar, maxts, 0.0);
    std::fill_n(xsphcor, maxts, 0.0);
    std::fill_n(ysphcor, maxts, 0.0);
    std::fill_n(zsphcor, maxts, 0.0);
    std::fill_n(rsph, maxts, 0.0);
    std::fill_n(nvert, maxts, 0);
    std::fill_n(vert, 30*maxts, 0.0);
    std::fill_n(centr, 30*maxts, 0.0);
    std::fill_n(isphe, maxts, 0);

    int nts = 0;
    int ntsirr = 0;

    // If there's an overflow in the number of spheres PEDRA will die.
    // The maximum number of spheres in PEDRA is set to 200 (primitive+additional)
    // so the integer here declared is just to have enough space C++ side to pass everything back.
    int maxAddedSpheres = 200;

    // Allocate vectors of size equal to nSpheres + maxAddedSpheres where maxAddedSpheres is the
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

    double * rin = radii_scratch.data();
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

    std::stringstream pedra;
    pedra << "PEDRA.OUT_" << suffix << "_" << ::getpid();
    char * f_pedra = const_cast<char *>(pedra.str().c_str());
    int len_f_pedra = std::strlen(f_pedra);
    // Go PEDRA, Go!
    TIMER_ON("GePolCavity::generatecavity_cpp");
    pedra_driver(&maxts, &maxsp, &maxvert,
                       xtscor, ytscor, ztscor, ar, xsphcor, ysphcor, zsphcor, rsph,
                       &nts, &ntsirr, &nSpheres_, &addedSpheres,
                       xe, ye, ze, rin, mass,
		       &averageArea, &probeRadius, &minimalRadius,
                       &nr_gen, &gen1, &gen2, &gen3,
                nvert, vert, centr, isphe, f_pedra, &len_f_pedra);
    TIMER_OFF("GePolCavity::generatecavity_cpp");

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

    // Now take care of updating the rest of the cavity info.
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
    // Check that no points are overlapping exactly
    // Do not perform float comparisons column by column.
    // Instead form differences between columns and evaluate if they differ
    // from zero by more than a fixed threshold.
    // The indices of the equal elements are gathered in a std::pair and saved into a std::vector
    double threshold = 1.0e-12;
    std::vector< std::pair<int, int> > equal_elements;
    for(int i = 0; i < nElements_; ++i) {
        for (int j = i + 1; j < nElements_; ++j) {
            Eigen::Vector3d difference = elementCenter_.col(i) - elementCenter_.col(j);
            if ( difference.isZero(threshold) ) {
                equal_elements.push_back(std::make_pair(i, j));
            }
        }
    }
    if (equal_elements.size() != 0) {
        // Not sure that printing the list of pairs is actually of any help...
        std::string list_of_pairs;
        for ( size_t i = 0; i < equal_elements.size(); ++i) {
            list_of_pairs += "(" + boost::lexical_cast<std::string>(equal_elements[i].first)
                             + ", " + boost::lexical_cast<std::string>(equal_elements[i].second) + ")\n";
        }
        // Prepare the error message:
        std::string message = boost::lexical_cast<std::string>(equal_elements.size()) +
                              " cavity finite element centers overlap exactly!\n" + list_of_pairs;
        PCMSOLVER_ERROR(message);
    }
    // Calculate normal vectors
    elementNormal_ = elementCenter_ - elementSphereCenter_;
    for( int i = 0; i < nElements_; ++i) {
        elementNormal_.col(i) /= elementNormal_.col(i).norm();
    }

    // Fill elements_ vector
    for (int i = 0; i < nElements_; ++i) {
        int i_off = i + 1;
        bool irr = false;
        // PEDRA puts the irreducible tesserae first
        if (i < nIrrElements_) irr = true;
        Sphere sph(elementSphereCenter_.col(i), elementRadius_(i));
        int nv = nvert[i];
        int isph = isphe[i]; // Back to C++ indexing starting from 0
        Eigen::Matrix3Xd vertices, arcs;
        vertices.resize(Eigen::NoChange, nv);
        arcs.resize(Eigen::NoChange, nv);
        // Populate vertices and arcs
        for (int j = 0; j < nv; ++j) {
            int j_off = (j + 1) * nElements_ - 1;
            for (int k = 0; k < 3; ++k) {
                int k_off = (k + 1) * nElements_ * nv;
                int offset = i_off + j_off + k_off;
                vertices(k, j) = vert[offset];
                arcs(k, j) = centr[offset];
            }
        }
        elements_.push_back(Element(nv, isph,
                    elementArea_(i),
                    elementCenter_.col(i),
                    elementNormal_.col(i),
                    irr, sph,
                    vertices, arcs));
    }

    // Clean-up
    delete[] xtscor;
    delete[] ytscor;
    delete[] ztscor;
    delete[] ar;
    delete[] xsphcor;
    delete[] ysphcor;
    delete[] zsphcor;
    delete[] rsph;
    delete[] nvert;
    delete[] vert;
    delete[] centr;
    delete[] mass;
    delete[] isphe;

    built = true;

    writeOFF(suffix);
}

void GePolCavity::writeOFF(const std::string & suffix)
{
    std::stringstream off;
    off << "cavity.off_" << suffix << "_" << ::getpid();

    std::ofstream fout;
    fout.open(off.str().c_str());

    int numv = 0;
    for (int i = 0; i < nElements_; ++i) {
        numv += elements_[i].nVertices();
    }
    fout << "COFF" << std::endl;
    fout << numv << " " << nElements_ << " " << numv << std::endl;

    int k = 0;
    double c1, c2, c3;
    Eigen::MatrixXi ivts = Eigen::MatrixXi::Zero(nElements_, 10);
    for (int i = 0; i < nElements_; ++i) {
        if (i == 0) fout << boost::format("# Sphere number %i\n") % elements_[i].iSphere();
        c1 = 1.0;
        c2 = 1.0;
        c3 = 1.0;
        for (int j = 0; j < elements_[i].nVertices(); ++j) {
            ivts(i, j) = k;
            k = k + 1;
            fout << boost::format("%20.14f\t%20.14f\t%20.14f\t%1.3f\t%1.3f\t%1.3f\t%1.3f   # Tess %i\n")
            % elements_[i].vertices()(0, j)
            % elements_[i].vertices()(1, j)
            % elements_[i].vertices()(2, j)
            % c1
            % c2
            % c3
            % 0.75
            % (i+1);
        }
    }
    for (int i = 0; i < nElements_; ++i) {
        fout << boost::format("%i ") % elements_[i].nVertices();
            for (int j = 0; j < elements_[i].nVertices(); ++j) {
                fout << boost::format("%i ") % ivts(i, j);
            }
        fout << std::endl;
    }

    fout.close();
}

std::ostream & GePolCavity::printCavity(std::ostream & os)
{
    os << "Cavity type: GePol" << std::endl;
    os << "Average area = " << averageArea << " AU^2" << std::endl;
    os << "Probe radius = " << probeRadius << " AU"   << std::endl;
    if (addedSpheres != 0) {
        os << "Addition of extra spheres enabled" << std::endl;
    }
    os << "Number of spheres = " << nSpheres_ << " [initial = " << nSpheres_ -
       addedSpheres << "; added = " << addedSpheres << "]" << std::endl;
    os << "Number of finite elements = " << nElements_;
    if (molecule_.pointGroup().nrGenerators() != 0) {
        os << "\nNumber of irreducible finite elements = " << nIrrElements_;
    }
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
