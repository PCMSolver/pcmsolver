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
 *     PCMSolver API, see: <http://pcmsolver.readthedocs.org/>
 */
/* pcmsolver_copyright_end */

#ifndef MEDDLE_HPP
#define MEDDLE_HPP

#include <string>

#include "Config.hpp"

#include <Eigen/Core>

#include <boost/container/flat_map.hpp>

class Cavity;
class IGreensFunction;
class Input;
struct PCMInput;
class PCMSolver;

#include "Input.hpp"
#include "utils/Symmetry.hpp"

/*! \file Meddle.hpp
 *  \author Roberto Di Remigio
 *  \date 2015
 */

/*! \namespace pcm */
namespace pcm {
    typedef boost::container::flat_map<std::string, Eigen::VectorXd> SurfaceFunctionMap;
    typedef SurfaceFunctionMap::value_type SurfaceFunctionPair;
    typedef SurfaceFunctionMap::const_iterator SurfaceFunctionMapConstIter;

    void printer(const std::string & message);
    void printer(const std::ostringstream & stream);
    void initMolecule(const Input & inp, const Symmetry & group,
            int nuclei, const Eigen::VectorXd & charges, const Eigen::Matrix3Xd & centers,
            Molecule & molecule);
    void initSpheresAtoms(const Input &, const Eigen::Matrix3Xd &, std::vector<Sphere> &);
    unsigned int pcmsolver_get_version(void) attribute(const);
    void print(const PCMInput &);

    /*! \class Meddle
     *  \brief Contains functions exposing an interface to the module internals
     */
    class Meddle __final
    {
        public:
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
            Meddle(pcmsolver_reader_t input_reading, int nr_nuclei, double charges[], double coordinates[], int symmetry_info[], const PCMInput & host_input);
            ~Meddle();
            /*! \brief Getter for the number of finite elements composing the molecular cavity
             *  \return the size of the cavity
             */
            PCMSolverIndex getCavitySize() const attribute(pure);
            /*! \brief Getter for the number of irreducible finite elements composing the molecular cavity
             *  \return the number of irreducible finite elements
             */
            PCMSolverIndex getIrreducibleCavitySize() const attribute(pure);
            /*! \brief Getter for the centers of the finite elements composing the molecular cavity
             *  \param[out] centers array holding the coordinates of the finite elements centers
             */
            void getCenters(double centers[]) const;
            /*! \brief Getter for the center of the i-th finite element
             *  \param[in] its index of the finite element
             *  \param[out] center array holding the coordinates of the finite element center
             */
            void getCenter(int its, double center[]) const;
            /*! \brief Getter for the areas/weights of the finite elements
             *  \param[in, out] context the PCM context object
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
            /*! \brief Computes response ASC given a MEP and the desired irreducible representation
             *  \param[in] mep_name label of the MEP surface function
             *  \param[in] asc_name label of the ASC surface function
             *  \param[in] irrep index of the desired irreducible representation
             *  If `Nonequilibrium = True` in the input, calculates a response
             *  ASC using the dynamic permittivity. Falls back to the solver with static permittivity
             *  otherwise.
             */
            void computeResponseASC(const char * mep_name, const char * asc_name, int irrep) const;
            /*! \brief Computes the polarization energy
             *  \param[in] mep_name label of the MEP surface function
             *  \param[in] asc_name label of the ASC surface function
             *  \return the polarization energy
             *  This function calculates the dot product of the given MEP and ASC vectors.
             */
            double computePolarizationEnergy(const char * mep_name, const char * asc_name) const;
            /*! \brief Retrieves data wrapped in a given surface function
             *  \param[in] size the size of the surface function
             *  \param[in] values the values wrapped in the surface function
             *  \param[in] name label of the surface function
             */
            void getSurfaceFunction(PCMSolverIndex size, double values[], const char * name) const;
            /*! \brief Sets a surface function given data and label
             *  \param[in] size the size of the surface function
             *  \param[in] values the values to be wrapped in the surface function
             *  \param[in] name label of the surface function
             */
            void setSurfaceFunction(PCMSolverIndex size, double values[], const char * name) const;
            /*! \brief Dumps all currently saved surface functions to NumPy arrays in .npy files
             */
            void saveSurfaceFunctions() const;
            /*! \brief Dumps a surface function to NumPy array in .npy file
             *  \param[in] name label of the surface function
             */
            void saveSurfaceFunction(const char * name) const;
            /*! \brief Loads a surface function from a .npy file
             *  \param[in] name label of the surface function
             */
            void loadSurfaceFunction(const char * name) const;
            /*! \brief Prints citation and set up information
             */
            void printInfo() const;
            /*! \brief Writes timing results for the API
             */
            void writeTimings() const;
        private:
            /*! Input object */
            Input input_;
            /*! Cavity */
            Cavity * cavity_;
            /*! Solver with static permittivity */
            PCMSolver * K_0_;
            /*! Solver with dynamic permittivity */
            PCMSolver * K_d_;
            /*! PCMSolver set up information */
            mutable std::ostringstream infoStream_;
            /*! Whether K_d_ was initialized */
            bool hasDynamic_;
            /*! SurfaceFunction map */
            mutable SurfaceFunctionMap functions_;
            /*! Initialize input_ */
            void initInput(pcmsolver_reader_t input_reading, int nr_nuclei, double charges[], double coordinates[], int symmetry_info[], const PCMInput & host_input);
            /*! Initialize cavity_ */
            void initCavity();
            /*! Initialize static solver K_0_ */
            void initStaticSolver();
            /*! Initialize dynamic solver K_d_ */
            void initDynamicSolver();
            /*! Collect info on medium */
            void mediumInfo(IGreensFunction * gf_i, IGreensFunction * gf_o) const;
    };
} /* end namespace pcm */

#endif /* MEDDLE_HPP */
