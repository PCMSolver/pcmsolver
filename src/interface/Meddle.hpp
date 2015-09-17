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
 *     PCMSolver API, see: <http://pcmsolver.github.io/pcmsolver-doc>
 */
/* pcmsolver_copyright_end */

#ifndef MEDDLE_HPP
#define MEDDLE_HPP

#include <string>

#include "Config.hpp"

class Cavity;
class Input;
class PCMSolver;

struct cavityInput;
struct greenInput;
struct solverInput;

#include "Input.hpp"
#include "SurfaceFunction.hpp"

namespace pcm {
    typedef function<int(void)> NrNucleiGetter;
    typedef function<void(double[], double[])> CoordinatesGetter;
    typedef function<void(const char *, size_t)> HostWriter;
    typedef function<void(int *, int *, int *, int *)> PointGroupSetter;
    typedef function<void(cavityInput &, solverInput &, greenInput &)> HostInput;
    typedef unordered_map<std::string, SurfaceFunction> SurfaceFunctionMap;

    void initMolecule(const Input & inp, const PointGroupSetter & set_group,
            int nuclei, const Eigen::VectorXd & charges, const Eigen::Matrix3Xd & centers,
            Molecule & molecule);
    void initSpheresAtoms(const Input &, const Eigen::Matrix3Xd &, std::vector<Sphere> &);
    unsigned int pcmsolver_get_version(void);

    class Meddle __final
    {
        public:
            Meddle(const NrNucleiGetter &, const CoordinatesGetter &, const HostWriter &,
                   const PointGroupSetter &);
            ~Meddle();
            size_t getCavitySize() const;
            size_t getIrreducibleCavitySize() const;
            void getCenters(double centers[]) const;
            void getCenter(int its, double center[]) const;
            void computeASC(const char * mep_name, const char * asc_name, int irrep) const;
            void computeResponseASC(const char * mep_name, const char * asc_name, int irrep) const;
            double computePolarizationEnergy(const char * mep_name, const char * asc_name) const;
            void getSurfaceFunction(size_t size, double values[], const char * name) const;
            void setSurfaceFunction(size_t size, double values[], const char * name) const;
            void saveSurfaceFunctions() const;
            void saveSurfaceFunction(const char * name) const;
            void loadSurfaceFunction(const char * name) const;
            void printInfo() const;
            void writeTimings() const;
        private:
            /*! Function to collect number of atoms in molecule */
            NrNucleiGetter nrNuclei_;
            /*! Function to collect atomic charges and coordinates */
            CoordinatesGetter chargesAndCoordinates_;
            /*! Function redirecting the output to the host */
            HostWriter hostWriter_;
            /*! Function setting the (Abelian) point group */
            PointGroupSetter pointGroup_;
            /*! Function reading input host-side */
            HostInput hostInputReader_;
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
            void initInput();
            /*! Initialize cavity_ */
            void initCavity();
            /*! Initialize solvers K_0_ and K_d_ */
            void initSolver();
            void printer(const std::string & message) const;
            void printer(const std::ostringstream & stream) const;
    };
} /* end namespace pcm */

#endif /* MEDDLE_HPP */
