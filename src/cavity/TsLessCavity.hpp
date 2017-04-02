/**
 * PCMSolver, an API for the Polarizable Continuum Model
 * Copyright (C) 2017 Roberto Di Remigio, Luca Frediani and collaborators.
 *
 * This file is part of PCMSolver.
 *
 * PCMSolver is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * PCMSolver is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with PCMSolver.  If not, see <http://www.gnu.org/licenses/>.
 *
 * For information on the complete list of contributors to the
 * PCMSolver API, see: <http://pcmsolver.readthedocs.io/>
 */

#ifndef TSLESSCAVITY_HPP
#define TSLESSCAVITY_HPP

#include <iosfwd>
#include <string>
#include <vector>

#include "Config.hpp"
#include "FCMangle.hpp"

namespace pcm {
struct CavityData;
} // namespace pcm

#include "ICavity.hpp"

/*! \file TsLessCavity.hpp
 *  \class TsLessCavity
 *  \brief A class for TsLess cavity.
 *  \author Roberto Di Remigio
 *  \date 2013
 *
 *  This class is a wrapper for the Fortran routines generating
 *  a tessellationless grid on the cavity surface. The original
 *  algoritm is described in \cite Pomelli2004
 */

namespace pcm {
namespace cavity {
class TsLessCavity __final : public ICavity {
public:
  TsLessCavity() {}
  TsLessCavity(const Molecule & molec,
               double a,
               double pr,
               double minR,
               double minD,
               int der);
  TsLessCavity(const Sphere & sph,
               double a,
               double pr,
               double minR,
               double minD,
               int der);
  TsLessCavity(const std::vector<Sphere> & sph,
               double a,
               double pr,
               double minR,
               double minD,
               int der);
  virtual ~TsLessCavity() {}
  friend std::ostream & operator<<(std::ostream & os, TsLessCavity & cavity) {
    return cavity.printCavity(os);
  }

private:
  double averageArea_;
  double probeRadius_;
  double minimalRadius_;
  double minDistance_;
  int derOrder_;
  int addedSpheres;
  virtual std::ostream & printCavity(std::ostream & os) __override;
  virtual void makeCavity() __override { build(10000, 200, 25000); }
  /*! \brief Driver for TsLess Fortran module.
   *  \param[in]   maxts maximum number of tesserae
   *  \param[in]   maxsp maximum number of spheres (original + added)
   *  \param[in] maxvert maximum number of vertices
   */
  void build(int maxts, int maxsp, int maxvert);
};

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
extern "C" void tsless_driver(int * maxts, int * maxsph, int * maxvert,
        int * nesfp, int * nts, int * ntsirr, int * addsph,
        double * xtscor, double * ytscor, double * ztscor, double * ar,
        double * xsphcor, double * ysphcor, double * zsphcor, double * rsph,
        double * xe, double * ye, double * ze, double * rin, double * masses,
        int * nr_gen, int * gen1, int * gen2, int * gen3,
        double * avgArea, double * dmin, int * nord, int * ifun, double * rsolv,
        double * work);

ICavity * createTsLessCavity(const CavityData & data);
} // namespace cavity
} // namespace pcm

#endif // TSLESSCAVITY_HPP
