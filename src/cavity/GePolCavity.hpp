/*
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

#pragma once

#include <iosfwd>
#include <string>
#include <vector>

#include "Config.hpp"

/*! \file GePolCavity.hpp*/

namespace pcm {
struct CavityData;
} // namespace pcm

#include "ICavity.hpp"

namespace pcm {
namespace cavity {
/*! \class GePolCavity
 *  \brief A class for GePol cavity.
 *  \author Krzysztof Mozgawa, Roberto Di Remigio
 *  \date 2011, 2016
 *
 *  This class is an interface to the Fortran code PEDRA for the generation
 *  of cavities according to the GePol algorithm.
 */
class GePolCavity __final : public ICavity {
public:
  GePolCavity() {}
  GePolCavity(const Molecule & molec,
              double a,
              double pr,
              double minR,
              const std::string & suffix = "");
  GePolCavity(const Sphere & sph,
              double a,
              double pr,
              double minR,
              const std::string & suffix = "");
  GePolCavity(const std::vector<Sphere> & sph,
              double a,
              double pr,
              double minR,
              const std::string & suffix = "");
  virtual ~GePolCavity() {}
  friend std::ostream & operator<<(std::ostream & os, GePolCavity & cavity) {
    return cavity.printCavity(os);
  }

private:
  double averageArea;
  double probeRadius;
  double minimalRadius;
  int addedSpheres;
  virtual std::ostream & printCavity(std::ostream & os) __override;
  virtual void makeCavity() __override {
    build(std::string("PEDRA.OUT"), 10000, 200, 25000);
  }
  /*! \brief Driver for PEDRA Fortran module.
   *  \param[in]  suffix for the cavity.off and PEDRA.OUT files, the PID will also be
   * added
   *  \param[in]   maxts maximum number of tesserae
   *  \param[in]   maxsp maximum number of spheres (original + added)
   *  \param[in] maxvert maximum number of vertices
   */
  void build(const std::string & suffix, int maxts, int maxsp, int maxvert);
  /*! \brief Writes the cavity.off file for visualizing the cavity
   *  \param[in]  suffix for the cavity.off
   *  The full name of the visualization file will be cavity.off_suffix_PID
   */
  void writeOFF(const std::string & suffix);
};

/*! \fn extern "C" void generatecavity_cpp(int * maxts, int * maxsph, int * maxvert,
 *                                 double * xtscor, double * ytscor, double * ztscor,
 * double * ar,
 *                                 double * xsphcor, double * ysphcor, double *
 * zsphcor, double * rsph,
 *                                 int * nts, int * ntsirr, int * nesfp, int *
 * addsph,
 *                                 double * xe, double * ye, double * ze, double *
 * rin,
 *                                 double * avgArea, double * rsolv, double * ret,
 *                                 int * nr_gen, int * gen1, int * gen2, int * gen3,
 *                                 int * nvert, double * vert, double * centr)
 *  \param[in] maxts maximum number of tesserae allowed
 *  \param[in] maxsph maximum number of spheres allowed
 *  \param[in] maxvert maximum number of vertices allowed
 *  \param[out] xtscor x-coordinate of tesserae centers (dimension maxts)
 *  \param[out] ytscor y-coordinate of tesserae centers (dimension maxts)
 *  \param[out] ztscor z-coordinate of tesserae centers (dimension maxts)
 *  \param[out] ar area of the tessera (dimension maxts)
 *  \param[out] xsphcor x-coordinate of the sphere center the tessera belongs to
 * (dimension maxts)
 *  \param[out] ysphcor y-coordinate of the sphere center the tessera belongs to
 * (dimension maxts)
 *  \param[out] zsphcor z-coordinate of the sphere center the tessera belongs to
 * (dimension maxts)
 *  \param[out] rsph radii of the sphere the tessera belongs to, i.e. its curvature
 * (dimension maxts)
 *  \param[out] nts number of generated tesserae
 *  \param[out] ntsirr number of generated irreducible tesserae
 *  \param[out] nesfp number of spheres (original + added)
 *  \param[out] addsph number of added spheres
 *  \param[out] xe x-coordinate of the sphere center (dimension nSpheres_ +
 * maxAddedSpheres)
 *  \param[out] ye y-coordinate of the sphere center (dimension nSpheres_ +
 * maxAddedSpheres)
 *  \param[out] ze z-coordinate of the sphere center (dimension nSpheres_ +
 * maxAddedSpheres)
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
extern "C" void generatecavity_cpp(int * maxts,
                                   int * maxsph,
                                   int * maxvert,
                                   double * xtscor,
                                   double * ytscor,
                                   double * ztscor,
                                   double * ar,
                                   double * xsphcor,
                                   double * ysphcor,
                                   double * zsphcor,
                                   double * rsph,
                                   int * nts,
                                   int * ntsirr,
                                   int * nesfp,
                                   int * addsph,
                                   double * xe,
                                   double * ye,
                                   double * ze,
                                   double * rin,
                                   double * masses,
                                   double * avgArea,
                                   double * rsolv,
                                   double * ret,
                                   int * nr_gen,
                                   int * gen1,
                                   int * gen2,
                                   int * gen3,
                                   int * nvert,
                                   double * vert,
                                   double * centr,
                                   int * isphe,
                                   const char * pedra,
                                   int * len_f_pedra);

ICavity * createGePolCavity(const CavityData & data);
} // namespace cavity
} // namespace pcm
