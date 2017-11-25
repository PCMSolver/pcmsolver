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

#include "pcmsolver.h"

#include <map>
#include <string>

#include "Config.hpp"

#include <Eigen/Core>

namespace pcm {
class ICavity;
class IGreensFunction;
class ISolver;
} // namespace pcm
struct PCMInput;

#include "Input.hpp"
#include "utils/Molecule.hpp"
#include "utils/Symmetry.hpp"

/*! \file Meddle.hpp
 *  \author Roberto Di Remigio
 *  \date 2015
 */

/*! \namespace pcm */
namespace pcm {
unsigned int pcmsolver_get_version(void) attribute(const);

namespace detail {
Molecule initMolecule(const Input & inp,
                      const Symmetry & group,
                      int nuclei,
                      const Eigen::VectorXd & charges,
                      const Eigen::Matrix3Xd & centers);
void initSpheresAtoms(const Input &,
                      const Eigen::Matrix3Xd &,
                      std::vector<Sphere> &);
void print(const PCMInput &);
} // namespace detail

/*! \class Meddle
 *  \brief Contains functions exposing an interface to the module internals
 *  \author Roberto Di Remigio
 *  \date 2015-2017
 */
class Meddle __final {
public:
  /*! \brief CTOR from Input object
      *  \param[in] input an Input object
      *  \param[in] write the global HostWriter object
      *  \warning This CTOR is meant to be used with the standalone
      *  executable only
      */
  Meddle(const Input & input, const HostWriter & write);
  /*! \brief CTOR from own input reader
      *  \param[in] inputFileName name of the parsed, machine-readable input file
      *  \param[in] write the global HostWriter object
      *  \warning This CTOR is meant to be used with the standalone
      *  executable only
      */
  Meddle(const std::string & inputFileName, const HostWriter & write);
  /*! \brief Constructor
   *  \param[in] input_reading input processing strategy
   *  \param[in] nr_nuclei     number of atoms in the molecule
   *  \param[in] charges       atomic charges
   *  \param[in] coordinates   atomic coordinates
   *  \param[in] symmetry_info molecular point group information
   *  \param[in] host_input    input to the module, as read by the host
   *
   *  The molecular point group information is passed as an array
   *  of 4 integers: number of generators, first, second and third generator
   *  respectively. Generators map to integers as in table :ref:`symmetry-ops`
   */
  Meddle(pcmsolver_reader_t input_reading,
         int nr_nuclei,
         double charges[],
         double coordinates[],
         int symmetry_info[],
         const PCMInput & host_input,
         const HostWriter & write);
  ~Meddle();
  /*! \brief Getter for the molecule object */
  Molecule molecule() const attribute(pure);
  /*! \brief Getter for the number of finite elements composing the molecular cavity
   *  \return the size of the cavity
   */
  PCMSolverIndex getCavitySize() const attribute(pure);
  /*! \brief Getter for the number of irreducible finite elements composing the
   * molecular cavity
   *  \return the number of irreducible finite elements
   */
  PCMSolverIndex getIrreducibleCavitySize() const attribute(pure);
  /*! \brief Getter for the centers of the finite elements composing the molecular
   * cavity
   *  \param[out] centers array holding the coordinates of the finite elements
   * centers
   */
  void getCenters(double centers[]) const;
  /*! \brief Getter for the center of the i-th finite element
   *  \param[in] its index of the finite element
   *  \param[out] center array holding the coordinates of the finite element center
   */
  void getCenter(int its, double center[]) const;
  /*! \brief Getter for the centers of the finite elements composing the molecular
   * cavity
   *  \return a matrix with the finite elements centers (dimensions 3 x number of
   * finite elements)
   */
  Eigen::Matrix3Xd getCenters() const attribute(pure);
  /*! \brief Getter for the areas/weights of the finite elements
   *  \param[out] areas array holding the weights/areas of the finite elements
   */
  void getAreas(double areas[]) const;
  /*! \brief Computes ASC given a MEP and the desired irreducible representation
   *  \param[in] mep_name label of the MEP surface function
   *  \param[in] asc_name label of the ASC surface function
   *  \param[in] irrep index of the desired irreducible representation
   *  The module uses the surface function concept to handle potentials
   *  and charges. Given labels for each, this function retrieves the MEP
   *  and computes the corresponding ASC.
   */
  void computeASC(const char * mep_name, const char * asc_name, int irrep) const;
  /*! \brief Computes response ASC given a MEP and the desired irreducible
   * representation
   *  \param[in] mep_name label of the MEP surface function
   *  \param[in] asc_name label of the ASC surface function
   *  \param[in] irrep index of the desired irreducible representation
   *  If `Nonequilibrium = True` in the input, calculates a response
   *  ASC using the dynamic permittivity. Falls back to the solver with static
   * permittivity
   *  otherwise.
   */
  void computeResponseASC(const char * mep_name,
                          const char * asc_name,
                          int irrep) const;
  /*! \brief Computes the polarization energy
   *  \param[in] mep_name label of the MEP surface function
   *  \param[in] asc_name label of the ASC surface function
   *  \return the polarization energy
   *  This function calculates the dot product of the given MEP and ASC vectors.
   */
  double computePolarizationEnergy(const char * mep_name,
                                   const char * asc_name) const;
  /*! \brief Getter for the ASC dipole
         *  \param[in] asc_name label of the ASC surface function
         *  \param[out] dipole  the Cartesian components of the ASC dipole moment
         *  \return the ASC dipole, i.e. \sqrt{\sum_i \mu_i^2}
         */
  double getASCDipole(const char * asc_name, double dipole[]) const;
  /*! \brief Retrieves data wrapped in a given surface function
   *  \param[in] size the size of the surface function
   *  \param[in] values the values wrapped in the surface function
   *  \param[in] name label of the surface function
   */
  void getSurfaceFunction(PCMSolverIndex size,
                          double values[],
                          const char * name) const;
  /*! \brief Sets a surface function given data and label
   *  \param[in] size the size of the surface function
   *  \param[in] values the values to be wrapped in the surface function
   *  \param[in] name label of the surface function
   */
  void setSurfaceFunction(PCMSolverIndex size,
                          double values[],
                          const char * name) const;
  /*! \brief Prints surface function contents to host output
   *  \param[in] name label of the surface function
   */
  void printSurfaceFunction(const char * name) const;
  /*! \brief Dumps all currently saved surface functions to NumPy arrays in .npy
   * files
   */
  void saveSurfaceFunctions() const;
  /*! \brief Dumps a surface function to NumPy array in .npy file
   *  \param[in] name label of the surface function
   *
   *  \note The name parameter is the name of the NumPy array file
   *  **without** .npy extension
   */
  void saveSurfaceFunction(const char * name) const;
  /*! \brief Loads a surface function from a .npy file
   *  \param[in] name label of the surface function
   *
   *  \note The name parameter is the name of the NumPy array file
   *  **without** .npy extension
   */
  void loadSurfaceFunction(const char * name) const;
  /*! \brief Prints citation and set up information
   */
  void printInfo() const;
  /*! \brief Writes timing results for the API
   */
  void writeTimings() const;

private:
  typedef std::map<std::string, Eigen::VectorXd> SurfaceFunctionMap;
  typedef SurfaceFunctionMap::value_type SurfaceFunctionPair;
  typedef SurfaceFunctionMap::iterator SurfaceFunctionMapIter;
  typedef SurfaceFunctionMap::const_iterator SurfaceFunctionMapConstIter;

  struct Printer {
    Printer(const HostWriter & hw) : writer_(hw) {}
    HostWriter writer_;
    void operator()(const std::string & message) const { writer_(message.c_str()); }
    void operator()(const std::ostringstream & stream) const {
      writer_(stream.str().c_str());
    }
  };
  /*! Output redirect-or to host program output */
  Printer hostWriter_;
  /*! Input object */
  Input input_;
  /*! Cavity */
  ICavity * cavity_;
  /*! Solver with static permittivity */
  ISolver * K_0_;
  /*! Solver with dynamic permittivity */
  ISolver * K_d_;
  /*! PCMSolver set up information */
  mutable std::ostringstream infoStream_;
  /*! Whether K_d_ was initialized */
  bool hasDynamic_;
  /*! SurfaceFunction map */
  mutable SurfaceFunctionMap functions_;
  /*! Common implemenation for the CTOR-s */
  void CTORBody();
  /*! \brief Initialize input_
   *  \param[in] input_reading input processing strategy
   *  \param[in] nr_nuclei     number of atoms in the molecule
   *  \param[in] charges       atomic charges
   *  \param[in] coordinates   atomic coordinates
   *  \param[in] symmetry_info molecular point group information
   *  \param[in] host_input    input to the module, as read by the host
   */
  void initInput(pcmsolver_reader_t input_reading,
                 int nr_nuclei,
                 double charges[],
                 double coordinates[],
                 int symmetry_info[],
                 const PCMInput & host_input);
  /*! Initialize cavity_ */
  void initCavity();
  /*! Initialize static solver K_0_ */
  void initStaticSolver();
  /*! Initialize dynamic solver K_d_ */
  void initDynamicSolver();
  /*! Collect info on medium */
  void mediumInfo(IGreensFunction * gf_i, IGreensFunction * gf_o) const;
};
} // namespace pcm
