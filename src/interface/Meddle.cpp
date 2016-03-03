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

#include "pcmsolver.h"
#include "PCMInput.h"
#include "Meddle.hpp"

#include <string>
#include <vector>

#include "Config.hpp"

#include <Eigen/Core>

#include <boost/foreach.hpp>

#include "cavity/Cavity.hpp"
#include "cavity/RegisterCavityToFactory.hpp"
#include "green/IGreensFunction.hpp"
#include "green/RegisterGreensFunctionToFactory.hpp"
#include "solver/PCMSolver.hpp"
#include "solver/RegisterSolverToFactory.hpp"
#include "utils/Atom.hpp"
#include "Citation.hpp"
#include "utils/cnpy.hpp"
#include "utils/PhysicalConstants.hpp"
#include "utils/Solvent.hpp"
#include "utils/Sphere.hpp"

#ifndef AS_TYPE
#define AS_TYPE(Type, Obj) reinterpret_cast<Type *>(Obj)
#endif

#ifndef AS_CTYPE
#define AS_CTYPE(Type, Obj) reinterpret_cast<const Type *>(Obj)
#endif

pcmsolver_context_t * pcmsolver_new(pcmsolver_reader_t input_reading, int
    nr_nuclei, double charges[], double coordinates[], int symmetry_info[],
    PCMInput * host_input)
{
  return AS_TYPE(pcmsolver_context_t, new pcm::Meddle(input_reading,
        nr_nuclei, charges, coordinates, symmetry_info, *host_input));
}

void pcmsolver_delete(pcmsolver_context_t * context)
{
  if (!context) return;
  delete AS_TYPE(pcm::Meddle, context);
}

bool pcmsolver_is_compatible_library(void)
{
  unsigned int major = (pcm::pcmsolver_get_version() >> 16);
  return (major == PROJECT_VERSION_MAJOR);
}

void pcmsolver_print(pcmsolver_context_t * context)
{
  AS_TYPE(pcm::Meddle, context)->printInfo();
}

size_t pcmsolver_get_cavity_size(pcmsolver_context_t * context)
{
  return (AS_TYPE(pcm::Meddle, context)->getCavitySize());
}

size_t pcmsolver_get_irreducible_cavity_size(pcmsolver_context_t * context)
{
  return (AS_TYPE(pcm::Meddle, context)->getIrreducibleCavitySize());
}

void pcmsolver_get_centers(pcmsolver_context_t * context, double centers[])
{
  AS_TYPE(pcm::Meddle, context)->getCenters(centers);
}

void pcmsolver_get_center(pcmsolver_context_t * context, int its, double center[])
{
  AS_TYPE(pcm::Meddle, context)->getCenter(its, center);
}

void pcmsolver_compute_asc(pcmsolver_context_t * context,
    const char * mep_name,
    const char * asc_name,
    int irrep)
{
  AS_TYPE(pcm::Meddle, context)->computeASC(mep_name, asc_name, irrep);
}

void pcmsolver_compute_response_asc(pcmsolver_context_t * context,
    const char * mep_name,
    const char * asc_name,
    int irrep)
{
  AS_TYPE(pcm::Meddle, context)->computeResponseASC(mep_name, asc_name, irrep);
}

double pcmsolver_compute_polarization_energy(pcmsolver_context_t * context,
    const char * mep_name,
    const char * asc_name)
{
  return (AS_TYPE(pcm::Meddle, context)->computePolarizationEnergy(mep_name, asc_name));
}

void pcmsolver_get_surface_function(pcmsolver_context_t * context,
    size_t size, double values[], const char * name)
{
  AS_TYPE(pcm::Meddle, context)->getSurfaceFunction(size, values, name);
}

void pcmsolver_set_surface_function(pcmsolver_context_t * context,
    size_t size, double values[], const char * name)
{
  AS_TYPE(pcm::Meddle, context)->setSurfaceFunction(size, values, name);
}

void pcmsolver_save_surface_functions(pcmsolver_context_t * context)
{
  AS_TYPE(pcm::Meddle, context)->saveSurfaceFunctions();
}

void pcmsolver_save_surface_function(pcmsolver_context_t * context, const char * name)
{
  AS_TYPE(pcm::Meddle, context)->saveSurfaceFunction(name);
}

void pcmsolver_load_surface_function(pcmsolver_context_t * context, const char * name)
{
  AS_TYPE(pcm::Meddle, context)->loadSurfaceFunction(name);
}

void pcmsolver_write_timings(pcmsolver_context_t * context)
{
  AS_TYPE(pcm::Meddle, context)->writeTimings();
}

namespace pcm {
  Meddle::Meddle(pcmsolver_reader_t input_reading, int nr_nuclei, double
      charges[], double coordinates[], int symmetry_info[], const PCMInput & host_input)
    : hasDynamic_(false)
  {
    initInput(input_reading, nr_nuclei, charges, coordinates, symmetry_info, host_input);
    initCavity();
    initStaticSolver();
    if (input_.isDynamic()) initDynamicSolver();
    // Reserve space for Tot-MEP/ASC, Nuc-MEP/ASC and Ele-MEP/ASC
    functions_.reserve(12);
  }

  Meddle::~Meddle()
  {
    delete cavity_;
    delete K_0_;
    if (hasDynamic_) delete K_d_;
  }

  size_t Meddle::getCavitySize() const
  {
    return cavity_->size();
  }

  size_t Meddle::getIrreducibleCavitySize() const
  {
    return cavity_->irreducible_size();
  }

  void Meddle::getCenters(double centers[]) const
  {
    Eigen::Map<Eigen::Matrix3Xd>(centers, 3, cavity_->size()) = cavity_->elementCenter();
  }

  void Meddle::getCenter(int its, double center[]) const
  {
    Eigen::Map<Eigen::Vector3d>(center, 3, 1) = cavity_->elementCenter(its-1);
  }

  double Meddle::computePolarizationEnergy(const char * mep_name, const char * asc_name) const
  {
    // Dot product of MEP and ASC surface function
    double energy = functions_[std::string(mep_name)].dot(functions_[std::string(asc_name)]);
    return (energy / 2.0);
  }

  void Meddle::computeASC(const char * mep_name, const char * asc_name, int irrep) const
  {
    std::string MEP(mep_name);
    std::string ASC(asc_name);

    // Get the proper iterators
    SurfaceFunctionMapConstIter iter_pot = functions_.find(MEP);
    Eigen::VectorXd asc = K_0_->computeCharge(iter_pot->second, irrep);
    // Renormalize
    asc /= double(cavity_->pointGroup().nrIrrep());
    // Insert it into the map
    if (functions_.count(ASC) == 1) { // Key in map already
      functions_[ASC] = asc;
    } else { // Create key-value pair
      functions_.insert(std::make_pair(ASC, asc));
    }
  }

  void Meddle::computeResponseASC(const char * mep_name, const char * asc_name, int irrep) const
  {
    std::string MEP(mep_name);
    std::string ASC(asc_name);

    // Get the proper iterators
    SurfaceFunctionMapConstIter iter_pot = functions_.find(MEP);
    Eigen::VectorXd asc(cavity_->size());
    if (hasDynamic_) {
      asc = K_d_->computeCharge(iter_pot->second, irrep);
    } else {
      asc = K_0_->computeCharge(iter_pot->second, irrep);
    }
    // Renormalize
    asc /= double(cavity_->pointGroup().nrIrrep());
    if (functions_.count(ASC) == 1) { // Key in map already
      functions_[ASC] = asc;
    } else { // Create key-value pair
      functions_.insert(std::make_pair(ASC, asc));
    }
  }

  void Meddle::getSurfaceFunction(size_t size, double values[], const char * name) const
  {
    if (cavity_->size() != size)
      PCMSOLVER_ERROR("You are trying to access a SurfaceFunction bigger than the cavity!");

    std::string functionName(name);

    SurfaceFunctionMapConstIter iter = functions_.find(functionName);
    if (iter == functions_.end())
      PCMSOLVER_ERROR("You are trying to access a non-existing SurfaceFunction.");

    Eigen::Map<Eigen::VectorXd>(values, size, 1) = iter->second;
  }

  void Meddle::setSurfaceFunction(size_t size, double values[], const char * name) const
  {
    if (cavity_->size() != size)
      PCMSOLVER_ERROR("You are trying to allocate a SurfaceFunction bigger than the cavity!");

    std::string functionName(name);
    Eigen::VectorXd func = Eigen::Map<Eigen::VectorXd>(values, size, 1);
    if (functions_.count(functionName) == 1) { // Key in map already
      functions_[functionName] = func;
    } else {
      functions_.insert(std::make_pair(functionName, func));
    }
  }

  void Meddle::saveSurfaceFunctions() const
  {
    printer("\nDumping surface functions to .npy files");
    BOOST_FOREACH(SurfaceFunctionPair pair, functions_) {
      unsigned int dim = static_cast<unsigned int>(pair.second.size());
      const unsigned int shape[] = {dim};
      std::string fname = pair.first + ".npy";
      cnpy::npy_save(fname, pair.second.data(), shape, 1, "w", true);
    }
  }

  void Meddle::saveSurfaceFunction(const char * name) const
  {
    std::string functionName(name);
    std::string fname = functionName + ".npy";

    SurfaceFunctionMapConstIter it = functions_.find(functionName);
    unsigned int dim = static_cast<unsigned int>(it->second.size());
    const unsigned int shape[] = {dim};
    cnpy::npy_save(fname, it->second.data(), shape, 1, "w", true);
  }

  void Meddle::loadSurfaceFunction(const char * name) const
  {
    std::string functionName(name);
    printer("\nLoading surface function " + functionName + " from .npy file");
    std::string fname = functionName + ".npy";
    cnpy::NpyArray raw_surfFunc = cnpy::npy_load(fname);
    unsigned int dim = raw_surfFunc.shape[0];
    if (dim != cavity_->size()) {
      PCMSOLVER_ERROR("Inconsistent dimension of loaded surface function!");
    } else {
      Eigen::VectorXd values = getFromRawBuffer<double>(dim, 1, raw_surfFunc.data);
      // Append to global map
      if (functions_.count(functionName) == 1) { // Key in map already
        functions_[functionName] = values;
      } else {
        functions_.insert(std::make_pair(functionName, values));
      }
    }
  }

  void Meddle::writeTimings() const
  {
    TIMER_DONE("pcmsolver.timer.dat");
  }

  void Meddle::initInput(pcmsolver_reader_t input_reading, int nr_nuclei, double charges[], double coordinates[], int symmetry_info[], const PCMInput & host_input)
  {
    if (input_reading) {
      input_ = Input(host_input);
    } else {
      input_ = Input("@pcmsolver.inp");
    }

    // 2. position and charges of atomic centers
    Eigen::VectorXd chg  = Eigen::Map<Eigen::VectorXd>(charges, nr_nuclei, 1);
    Eigen::Matrix3Xd centers = Eigen::Map<Eigen::Matrix3Xd>(coordinates, 3, nr_nuclei);

    if (input_.mode() != "EXPLICIT") {
      Molecule molec;
      Symmetry pg = buildGroup(symmetry_info[0], symmetry_info[1], symmetry_info[2], symmetry_info[3]);
      initMolecule(input_, pg, nr_nuclei, chg, centers, molec);
      input_.molecule(molec);
    }

    infoStream_ << std::endl;
    infoStream_ << "~~~~~~~~~~ PCMSolver ~~~~~~~~~~" << std::endl;
    infoStream_ << "Using CODATA " << input_.CODATAyear() << " set of constants." << std::endl;
    infoStream_ << "Input parsing done " << input_.providedBy() << std::endl;
  }

  void Meddle::initCavity()
  {
    cavity_ = Factory<Cavity, cavityData>::TheFactory().create(input_.cavityType(), input_.cavityParams());
    cavity_->saveCavity();

    infoStream_ << "========== Cavity " << std::endl;
    infoStream_ << *cavity_ << std::endl;
  }

  void Meddle::initStaticSolver()
  {
    IGreensFunction * gf_i = Factory<IGreensFunction, greenData>::TheFactory().create(input_.greenInsideType(),
        input_.insideGreenParams());
    IGreensFunction * gf_o = Factory<IGreensFunction, greenData>::TheFactory().create(input_.greenOutsideType(),
        input_.outsideStaticGreenParams());
    std::string modelType = input_.solverType();
    K_0_ = Factory<PCMSolver, solverData>::TheFactory().create(modelType, input_.solverParams());
    K_0_->buildSystemMatrix(*cavity_, *gf_i, *gf_o);

    infoStream_ << "========== Static solver " << std::endl;
    infoStream_ << *K_0_ << std::endl;
    mediumInfo(gf_i, gf_o);
    delete gf_o;
    delete gf_i;
  }

  void Meddle::initDynamicSolver()
  {
    IGreensFunction * gf_i = Factory<IGreensFunction, greenData>::TheFactory().create(input_.greenInsideType(),
        input_.insideGreenParams());
    IGreensFunction * gf_o = Factory<IGreensFunction, greenData>::TheFactory().create(input_.greenOutsideType(),
        input_.outsideDynamicGreenParams());
    std::string modelType = input_.solverType();
    K_d_ = Factory<PCMSolver, solverData>::TheFactory().create(modelType, input_.solverParams());
    K_d_->buildSystemMatrix(*cavity_, *gf_i, *gf_o);
    hasDynamic_ = true;

    infoStream_ << "========== Dynamic solver " << std::endl;
    infoStream_ << *K_d_ << std::endl;
    mediumInfo(gf_i, gf_o);
    delete gf_o;
    delete gf_i;
  }

  void Meddle::mediumInfo(IGreensFunction * gf_i, IGreensFunction * gf_o) const
  {
    infoStream_ << "============ Medium " << std::endl;
    if (input_.fromSolvent()) {
      infoStream_ << "Medium initialized from solvent built-in data." << std::endl;
      Solvent solvent = input_.solvent();
      infoStream_ << solvent << std::endl;
    }
    std::stringstream tmp;
    tmp << ".... Inside " << std::endl;
    tmp << *gf_i << std::endl;
    tmp << ".... Outside " << std::endl;
    tmp << *gf_o;
    infoStream_ << tmp.str() << std::endl;
  }

  void Meddle::printInfo() const
  {
    printer(citation_message());
    printer(infoStream_);
  }

  void printer(const std::string & message)
  {
    const char * message_C = message.c_str();
    host_writer(message_C, std::strlen(message_C));
  }

  void printer(const std::ostringstream & stream)
  {
    std::string message = stream.str();
    const char * message_C = message.c_str();
    host_writer(message_C, std::strlen(message_C));
  }

  void initMolecule(const Input & inp, const Symmetry & pg,
      int nuclei, const Eigen::VectorXd & charges, const Eigen::Matrix3Xd & centers,
      Molecule & molecule)
  {
    bool scaling = inp.scaling();
    std::string set = inp.radiiSet();
    double factor = angstromToBohr(inp.CODATAyear());
    std::vector<Atom> radiiSet, atoms;
    if ( set == "UFF" ) {
      radiiSet = initUFF();
    } else {
      radiiSet = initBondi();
    }
    std::vector<Sphere> spheres;
    for (int i = 0; i < charges.size(); ++i) {
      int index = int(charges(i)) - 1;
      atoms.push_back(radiiSet[index]);
      double radius = radiiSet[index].radius * factor;
      if (scaling) {
        radius *= radiiSet[index].radiusScaling;
      }
      spheres.push_back(Sphere(centers.col(i), radius));
    }
    Eigen::VectorXd masses = Eigen::VectorXd::Zero(nuclei);
    for (int i = 0; i < masses.size(); ++i) {
      masses(i) = atoms[i].mass;
    }
    // Based on the creation mode (Implicit or Atoms)
    // the spheres list might need postprocessing
    if (inp.mode() == "ATOMS") {
      initSpheresAtoms(inp, centers, spheres);
    }

    // OK, now get molecule
    molecule = Molecule(nuclei, charges, masses, centers, atoms, spheres, pg);
    // Check that all atoms have a radius attached
    std::vector<Atom>::const_iterator res =
      std::find_if(atoms.begin(), atoms.end(), invalid);
    if (res != atoms.end()) {
      std::ostringstream print_mol;
      print_mol << molecule << std::endl;
      printer(print_mol);
      PCMSOLVER_ERROR("Some atoms do not have a radius attached. Please specify a radius for all atoms!");
    }
  }

  void initSpheresAtoms(const Input & inp, const Eigen::Matrix3Xd & sphereCenter_,
      std::vector<Sphere> & spheres_)
  {
    // Loop over the atomsInput array to get which atoms will have a user-given radius
    for (size_t i = 0; i < inp.atoms().size(); ++i) {
      size_t index = inp.atoms(i) - 1; // -1 to go from human readable to machine readable
      // Put the new Sphere in place of the implicit-generated one
      spheres_[index] = Sphere(sphereCenter_.col(index), inp.radii(i));
    }
  }

  unsigned int pcmsolver_get_version(void)
  {
    return PCMSOLVER_VERSION;
  }

  void print(const PCMInput & inp)
  {
    std::cout << "cavity type " << std::string(inp.cavity_type) << std::endl;
    std::cout << "patch level " << inp.patch_level << std::endl;
    std::cout << "coarsity " << inp.coarsity << std::endl;
    std::cout << "area " << inp.area << std::endl;
    std::cout << "min distance " << inp.min_distance << std::endl;
    std::cout << "der order " << inp.der_order << std::endl;
    std::cout << "scaling " << inp.scaling << std::endl;
    std::cout << "radii set " << std::string(inp.radii_set) << std::endl;
    std::cout << "restart name " << std::string(inp.restart_name) << std::endl;
    std::cout << "min radius " << inp.min_radius << std::endl;
    std::cout << "solver type " << std::string(inp.solver_type) << std::endl;
    std::cout << "solvent " << std::string(inp.solvent) << std::endl;
    std::cout << "equation type " << std::string(inp.equation_type) << std::endl;
    std::cout << "correction " << inp.correction << std::endl;
    std::cout << "probe_radius " << inp.probe_radius << std::endl;
    std::cout << "inside type " << std::string(inp.inside_type) << std::endl;
    std::cout << "outside type " << std::string(inp.outside_type) << std::endl;
    std::cout << "epsilon outside " << inp.outside_epsilon << std::endl;
  }
} /* end namespace pcm */
