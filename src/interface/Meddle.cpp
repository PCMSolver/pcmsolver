/*
 * PCMSolver, an API for the Polarizable Continuum Model
 * Copyright (C) 2020 Roberto Di Remigio, Luca Frediani and contributors.
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

#include "Meddle.hpp"
#include "PCMInput.h"
#include "pcmsolver.h"

#include <string>
#include <vector>

#include "Config.hpp"

#include <Eigen/Core>

#include "bi_operators/BIOperatorData.hpp"
#include "bi_operators/BoundaryIntegralOperator.hpp"
#include "cavity/Cavity.hpp"
#include "cavity/CavityData.hpp"
#include "green/Green.hpp"
#include "green/GreenData.hpp"
#include "mmfq/FQOhno.hpp"
#include "solver/Solver.hpp"
#include "solver/SolverData.hpp"

#include "Citation.hpp"
#include "VersionInfo.hpp"
#include "utils/Atom.hpp"
#include "utils/Factory.hpp"
#include "utils/Solvent.hpp"
#include "utils/Sphere.hpp"
#include "utils/cnpy.hpp"

#ifndef AS_TYPE
#define AS_TYPE(Type, Obj) reinterpret_cast<Type *>(Obj)
#endif

#ifndef AS_CTYPE
#define AS_CTYPE(Type, Obj) reinterpret_cast<const Type *>(Obj)
#endif

pcmsolver_context_t * pcmsolver_new(pcmsolver_reader_t input_reading,
                                    int nr_nuclei,
                                    double charges[],
                                    double coordinates[],
                                    int symmetry_info[],
                                    PCMInput * host_input,
                                    HostWriter writer) {
  return pcmsolver_new_v1112(input_reading,
                             nr_nuclei,
                             charges,
                             coordinates,
                             symmetry_info,
                             "@pcmsolver.inp",
                             host_input,
                             writer);
}

pcmsolver_context_t * pcmsolver_new_v1112(pcmsolver_reader_t input_reading,
                                          int nr_nuclei,
                                          double charges[],
                                          double coordinates[],
                                          int symmetry_info[],
                                          const char * parsed_fname,
                                          PCMInput * host_input,
                                          HostWriter writer) {
  if (input_reading) {
    // Use host_input data structured as passed from host
    return AS_TYPE(
        pcmsolver_context_t,
        new pcm::Meddle(
            nr_nuclei, charges, coordinates, symmetry_info, *host_input, writer));
  } else {
    // Use the parsed PCMSolver input file parsed_fname, as found on disk
    return AS_TYPE(
        pcmsolver_context_t,
        new pcm::Meddle(
            nr_nuclei, charges, coordinates, symmetry_info, writer, parsed_fname));
  }
}

pcmsolver_context_t * pcmsolver_new_read_host(int nr_nuclei,
                                              double charges[],
                                              double coordinates[],
                                              int symmetry_info[],
                                              HostWriter writer) {
  return AS_TYPE(
      pcmsolver_context_t,
      new pcm::Meddle(nr_nuclei, charges, coordinates, symmetry_info, writer));
}

PCMInput pcmsolver_default_input() {
  PCMInput host_input;

  // These parameters would be set by the host input reading
  // Length and area parameters are all assumed to be in Angstrom,
  // the module will convert to Bohr internally
  strcpy(host_input.cavity_type, "gepol");
  host_input.patch_level = 2;
  host_input.coarsity = 0.5;
  host_input.area = 0.2;
  host_input.min_distance = 0.1;
  host_input.der_order = 4;
  host_input.scaling = true;
  strcpy(host_input.radii_set, "bondi");
  strcpy(host_input.restart_name, "cavity.npz");
  host_input.min_radius = 100.0;

  strcpy(host_input.solver_type, "iefpcm");
  strcpy(host_input.solvent, "water");
  strcpy(host_input.equation_type, "secondkind");
  host_input.correction = 0.0;
  host_input.probe_radius = 1.0;

  strcpy(host_input.inside_type, "vacuum");
  host_input.outside_epsilon = 1.0;
  strcpy(host_input.outside_type, "uniformdielectric");

  return host_input;
}

namespace pcm {
void Meddle::CTORBody() {
  // Write PCMSolver output header
  infoStream_ << "~~~~~~~~~~ PCMSolver ~~~~~~~~~~\n";
  infoStream_ << "Using CODATA " << input_.CODATAyear() << " set of constants."
              << std::endl;
  infoStream_ << "Input parsing done " << input_.providedBy() << std::endl;

  if (input_.isFQ()) { /* MMFQ calculation */
    hasFQ_ = true;
    TIMER_ON("Meddle::initMMFQ");
    initMMFQ();
    TIMER_OFF("Meddle::initMMFQ");
  } else { /* Pure PCM calculation */
    TIMER_ON("Meddle::initCavity");
    initCavity();
    TIMER_OFF("Meddle::initCavity");

    TIMER_ON("Meddle::initStaticSolver");
    initStaticSolver();
    TIMER_OFF("Meddle::initStaticSolver");

    if (input_.isDynamic()) {
      TIMER_ON("Meddle::initDynamicSolver");
      initDynamicSolver();
      TIMER_OFF("Meddle::initDynamicSolver");
    }
  }
}

Meddle::Meddle(const Input & input, const HostWriter & write)
    : hostWriter_(write),
      input_(input),
      host_input_(pcmsolver_default_input()),
      cavity_(nullptr),
      K_0_(nullptr),
      K_d_(nullptr),
      FQ_(nullptr),
      hasDynamic_(false),
      hasFQ_(false) {
  input_.initMolecule();
  CTORBody();
}

Meddle::Meddle(const std::string & inputFileName, const HostWriter & write)
    : hostWriter_(write),
      input_(Input(inputFileName)),
      host_input_(pcmsolver_default_input()),
      cavity_(nullptr),
      K_0_(nullptr),
      K_d_(nullptr),
      FQ_(nullptr),
      hasDynamic_(false),
      hasFQ_(false) {
  input_.initMolecule();
  CTORBody();
}

Meddle::Meddle(int nr_nuclei,
               double charges[],
               double coordinates[],
               int symmetry_info[],
               const HostWriter & write,
               const std::string & inputFileName)
    : hostWriter_(write),
      input_(Input(inputFileName.empty() ? "@pcmsolver.inp" : inputFileName)),
      host_input_(pcmsolver_default_input()),
      cavity_(nullptr),
      K_0_(nullptr),
      K_d_(nullptr),
      FQ_(nullptr),
      hasDynamic_(false),
      hasFQ_(false) {
  TIMER_ON("Meddle::initInput");
  initInput(nr_nuclei, charges, coordinates, symmetry_info);
  TIMER_OFF("Meddle::initInput");

  CTORBody();
}

Meddle::Meddle(int nr_nuclei,
               double charges[],
               double coordinates[],
               int symmetry_info[],
               const PCMInput & host_input,
               const HostWriter & write)
    : hostWriter_(write),
      input_(Input(host_input)),
      host_input_(host_input),
      cavity_(nullptr),
      K_0_(nullptr),
      K_d_(nullptr),
      FQ_(nullptr),
      hasDynamic_(false),
      hasFQ_(false) {
  TIMER_ON("Meddle::initInput");
  initInput(nr_nuclei, charges, coordinates, symmetry_info);
  TIMER_OFF("Meddle::initInput");

  CTORBody();
}

Meddle::Meddle(int nr_nuclei,
               double charges[],
               double coordinates[],
               int symmetry_info[],
               const HostWriter & write)
    : hostWriter_(write),
      input_(Input(pcmsolver_default_input())),
      host_input_(pcmsolver_default_input()),
      cavity_(nullptr),
      K_0_(nullptr),
      K_d_(nullptr),
      FQ_(nullptr),
      hasDynamic_(false),
      hasFQ_(false) {
  // This one does a deferred initialization:
  // The CTORBody function is not called at this point, but during refresh
  TIMER_ON("Meddle::initInput");
  initInput(nr_nuclei, charges, coordinates, symmetry_info, true);
  TIMER_OFF("Meddle::initInput");
}
} // namespace pcm

void pcmsolver_delete(pcmsolver_context_t * context) {
  if (!context)
    return;
  delete AS_TYPE(pcm::Meddle, context);
}
pcm::Meddle::~Meddle() {
  if (hasFQ_) { /* MMFQ calculation */
    delete FQ_;
  } else { /* Pure PCM calculation */
    delete cavity_;
    delete K_0_;
    if (hasDynamic_)
      delete K_d_;
  }
}

PCMSolverIndex pcmsolver_get_cavity_size(pcmsolver_context_t * context) {
  return (AS_CTYPE(pcm::Meddle, context)->getCavitySize());
}
PCMSolverIndex pcm::Meddle::getCavitySize() const { return std::get<0>(size_); }

PCMSolverIndex pcmsolver_get_irreducible_cavity_size(pcmsolver_context_t * context) {
  return (AS_CTYPE(pcm::Meddle, context)->getIrreducibleCavitySize());
}
PCMSolverIndex pcm::Meddle::getIrreducibleCavitySize() const {
  return std::get<1>(size_);
}

void pcmsolver_get_centers(pcmsolver_context_t * context, double centers[]) {
  TIMER_ON("pcmsolver_get_centers");
  AS_CTYPE(pcm::Meddle, context)->getCenters(centers);
  TIMER_OFF("pcmsolver_get_centers");
}
void pcm::Meddle::getCenters(double centers[]) const {
  TIMER_ON("Meddle::getCenters");
  if (hasFQ_) { /* MMFQ calculation */
    PCMSolverIndex size = input_.fragments().sites.cols();
    Eigen::Map<Eigen::Matrix3Xd>(centers, 3, size) = input_.fragments().sites;
  } else { /* Pure PCM calculation */
    Eigen::Map<Eigen::Matrix3Xd>(centers, 3, cavity_->size()) =
        cavity_->elementCenter();
  }
  TIMER_OFF("Meddle::getCenters");
}

void pcmsolver_get_center(pcmsolver_context_t * context, int its, double center[]) {
  AS_CTYPE(pcm::Meddle, context)->getCenter(its, center);
}
void pcm::Meddle::getCenter(int its, double center[]) const {
  if (hasFQ_) { /* MMFQ calculation */
    Eigen::Map<Eigen::Matrix3Xd>(center, 3, 1) =
        input_.fragments().sites.col(its - 1);
  } else { /* Pure PCM calculation */
    Eigen::Map<Eigen::Vector3d>(center, 3, 1) = cavity_->elementCenter(its - 1);
  }
}

void pcmsolver_get_areas(pcmsolver_context_t * context, double areas[]) {
  AS_CTYPE(pcm::Meddle, context)->getAreas(areas);
}
void pcm::Meddle::getAreas(double areas[]) const {
  Eigen::Map<Eigen::VectorXd>(areas, cavity_->size(), 1) = cavity_->elementArea();
}

double pcmsolver_compute_polarization_energy(pcmsolver_context_t * context,
                                             const char * mep_name,
                                             const char * asc_name) {
  return (
      AS_CTYPE(pcm::Meddle, context)
          ->computePolarizationEnergy(std::string(mep_name), std::string(asc_name)));
}
double pcm::Meddle::computePolarizationEnergy(const std::string & mep_name,
                                              const std::string & asc_name) const {
  double energy = 0.0;
  if (hasFQ_) { /* MMFQ calculation */
    // Dot product of MEP + Electronegativities times Fluctuating charges
    energy = (functions_.at(mep_name) + input_.fragments().chi)
                 .dot(functions_.at(asc_name));
    if (input_.isNonPolarizable()) { /* HACK Nonpolarizable doesn't need 1/2 */
      energy *= 2.0;
    }
  } else { /* Pure PCM calculation */
    // Dot product of MEP and ASC surface function
    energy = functions_.at(mep_name).dot(functions_.at(asc_name));
  }
  return (energy / 2.0);
}

double pcmsolver_get_asc_dipole(pcmsolver_context_t * context,
                                const char * asc_name,
                                double dipole[]) {
  return (
      AS_CTYPE(pcm::Meddle, context)->getASCDipole(std::string(asc_name), dipole));
}
double pcm::Meddle::getASCDipole(const std::string & asc_name,
                                 double dipole[]) const {
  Eigen::Vector3d asc_dipole = cavity_->elementCenter() * functions_.at(asc_name);
  // Bind to host-allocated array
  Eigen::Map<Eigen::Vector3d>(dipole, 3, 1) = asc_dipole;
  return asc_dipole.norm();
}

void pcmsolver_compute_asc(pcmsolver_context_t * context,
                           const char * mep_name,
                           const char * asc_name,
                           int irrep) {
  TIMER_ON("pcmsolver_compute_asc");
  AS_TYPE(pcm::Meddle, context)
      ->computeASC(std::string(mep_name), std::string(asc_name), irrep);
  TIMER_OFF("pcmsolver_compute_asc");
}
void pcm::Meddle::computeASC(const std::string & mep_name,
                             const std::string & asc_name,
                             int irrep) {
  // Get the proper iterators
  SurfaceFunctionMapConstIter iter_pot = functions_.find(mep_name);
  Eigen::VectorXd asc = Eigen::VectorXd::Zero(iter_pot->second.size());
  if (hasFQ_) {                      /* MMFQ calculation */
    if (input_.isNonPolarizable()) { /* HACK We store point charges in the eta vector
                                      */
      asc = input_.fragments().eta;
    } else {
      asc = FQ_->computeCharge(iter_pot->second);
    }
  } else { /* Pure PCM calculation */
    asc = K_0_->computeCharge(iter_pot->second, irrep);
    // Renormalize
    asc /= double(cavity_->pointGroup().nrIrrep());
  }
  // Insert it into the map
  if (functions_.count(asc_name) == 1) { // Key in map already
    functions_[asc_name] = asc;
  } else { // Create key-value pair
    functions_.insert(std::make_pair(asc_name, asc));
  }
}

void pcmsolver_compute_response_asc(pcmsolver_context_t * context,
                                    const char * mep_name,
                                    const char * asc_name,
                                    int irrep) {
  TIMER_ON("pcmsolver_compute_response_asc");
  AS_TYPE(pcm::Meddle, context)
      ->computeResponseASC(std::string(mep_name), std::string(asc_name), irrep);
  TIMER_OFF("pcmsolver_compute_response_asc");
}
void pcm::Meddle::computeResponseASC(const std::string & mep_name,
                                     const std::string & asc_name,
                                     int irrep) {
  // Get the proper iterators
  SurfaceFunctionMapConstIter iter_pot = functions_.find(mep_name);
  Eigen::VectorXd asc = Eigen::VectorXd::Zero(iter_pot->second.size());
  if (hasFQ_) {                       /* MMFQ calculation */
    if (!input_.isNonPolarizable()) { /* HACK Can we do it more cleanly/clearly? */
      // Do NOT add classical (electronegativities) contributions to RHS
      asc = FQ_->computeCharge(iter_pot->second, false);
    }
  } else { /* Pure PCM calculation */
    if (hasDynamic_) {
      asc = K_d_->computeCharge(iter_pot->second, irrep);
    } else {
      asc = K_0_->computeCharge(iter_pot->second, irrep);
    }
    // Renormalize
    asc /= double(cavity_->pointGroup().nrIrrep());
  }
  if (functions_.count(asc_name) == 1) { // Key in map already
    functions_[asc_name] = asc;
  } else { // Create key-value pair
    functions_.insert(std::make_pair(asc_name, asc));
  }
}

void pcmsolver_get_surface_function(pcmsolver_context_t * context,
                                    PCMSolverIndex size,
                                    double values[],
                                    const char * name) {
  TIMER_ON("pcmsolver_get_surface_function");
  AS_CTYPE(pcm::Meddle, context)
      ->getSurfaceFunction(size, values, std::string(name));
  TIMER_OFF("pcmsolver_get_surface_function");
}
void pcm::Meddle::getSurfaceFunction(PCMSolverIndex size,
                                     double values[],
                                     const std::string & name) const {
  if (std::get<0>(size_) != size)
    PCMSOLVER_ERROR("The " + name + " SurfaceFunction is bigger than the cavity!");

  SurfaceFunctionMapConstIter iter = functions_.find(name);
  if (iter == functions_.end())
    PCMSOLVER_ERROR("The " + name + " SurfaceFunction does not exist.");

  Eigen::Map<Eigen::VectorXd>(values, size, 1) = iter->second;
}

void pcmsolver_set_surface_function(pcmsolver_context_t * context,
                                    PCMSolverIndex size,
                                    double values[],
                                    const char * name) {
  TIMER_ON("pcmsolver_set_surface_function");
  AS_TYPE(pcm::Meddle, context)->setSurfaceFunction(size, values, std::string(name));
  TIMER_OFF("pcmsolver_set_surface_function");
}
void pcm::Meddle::setSurfaceFunction(PCMSolverIndex size,
                                     double values[],
                                     const std::string & name) {
  if (std::get<0>(size_) != size)
    PCMSOLVER_ERROR("The " + name + " SurfaceFunction is bigger than the cavity!");

  Eigen::VectorXd func = Eigen::Map<Eigen::VectorXd>(values, size, 1);
  if (functions_.count(name) == 1) { // Key in map already
    functions_[name] = func;
  } else {
    functions_.insert(std::make_pair(name, func));
  }
}

void pcmsolver_print_surface_function(pcmsolver_context_t * context,
                                      const char * name) {
  AS_CTYPE(pcm::Meddle, context)->printSurfaceFunction(std::string(name));
}
void pcm::Meddle::printSurfaceFunction(const std::string & name) const {
  if (functions_.count(name) == 1) { // Key in map already
    std::ostringstream print_sf;
    Eigen::IOFormat fmt(Eigen::FullPrecision);
    print_sf << functions_.at(name).format(fmt) << std::endl;
    hostWriter_(print_sf);
  } else {
    PCMSOLVER_ERROR("You are trying to print a nonexistent SurfaceFunction!");
  }
}

void pcmsolver_save_surface_functions(pcmsolver_context_t * context) {
  AS_CTYPE(pcm::Meddle, context)->saveSurfaceFunctions();
}
void pcm::Meddle::saveSurfaceFunctions() const {
  hostWriter_("\nDumping surface functions to .npy files");
  for (auto sf_pair : functions_) {
    cnpy::custom::npy_save(sf_pair.first + ".npy", sf_pair.second);
  }
}

void pcmsolver_save_surface_function(pcmsolver_context_t * context,
                                     const char * name) {
  AS_CTYPE(pcm::Meddle, context)->saveSurfaceFunction(std::string(name));
}
void pcm::Meddle::saveSurfaceFunction(const std::string & name) const {
  SurfaceFunctionMapConstIter it = functions_.find(name);
  cnpy::custom::npy_save(name + ".npy", it->second);
}

void pcmsolver_load_surface_function(pcmsolver_context_t * context,
                                     const char * name) {
  AS_TYPE(pcm::Meddle, context)->loadSurfaceFunction(std::string(name));
}
void pcm::Meddle::loadSurfaceFunction(const std::string & name) {
  hostWriter_("\nLoading surface function " + name + " from .npy file");
  Eigen::VectorXd values = cnpy::custom::npy_load<double>(name + ".npy");
  if (values.size() != std::get<0>(size_))
    PCMSOLVER_ERROR("The loaded " + name +
                    " surface function is bigger than the cavity!");
  // Append to global map
  if (functions_.count(name) == 1) { // Key in map already
    functions_[name] = values;
  } else {
    functions_.insert(std::make_pair(name, values));
  }
}

void pcmsolver_write_timings(pcmsolver_context_t * context) {
  AS_CTYPE(pcm::Meddle, context)->writeTimings();
}
void pcm::Meddle::writeTimings() const { TIMER_DONE("pcmsolver.timer.dat"); }

void pcmsolver_refresh(pcmsolver_context_t * context) {
  AS_TYPE(pcm::Meddle, context)->refresh();
}
void pcm::Meddle::refresh() {
  // Gather info to refresh Molecule
  // FIXME I need to refresh it because scaling of radii is done there -.-
  Symmetry pg = molecule().pointGroup();
  size_t nuclei = molecule().nAtoms();
  Eigen::VectorXd chgs = molecule().charges();
  Eigen::MatrixXd cnts = molecule().geometry();
  // Refresh input_
  input_ = Input(host_input_);
  // Refresh Molecule
  input_.molecule(detail::initMolecule(input_, pg, nuclei, chgs, cnts, false));
  CTORBody();
}

void pcmsolver_set_bool_option(pcmsolver_context_t * context,
                               const char * parameter,
                               bool value) {
  AS_TYPE(pcm::Meddle, context)->setInputOption(std::string(parameter), value);
}
void pcm::Meddle::setInputOption(std::string parameter, bool value) {
  detail::PCMInputFields p = detail::string_to_enum(parameter);
  switch (p) {
    case detail::PCMInputFields::scaling:
      host_input_.scaling = value;
      break;
    default:
      std::ostringstream errmsg;
      errmsg << "Unknown parameter " << parameter << std::endl;
      errmsg << " See http://pcmsolver.readthedocs.org/en/latest/users/input.html"
             << std::endl;
      PCMSOLVER_ERROR(errmsg.str());
  }
}

void pcmsolver_set_int_option(pcmsolver_context_t * context,
                              const char * parameter,
                              int value) {
  AS_TYPE(pcm::Meddle, context)->setInputOption(std::string(parameter), value);
}
void pcm::Meddle::setInputOption(std::string parameter, int value) {
  detail::PCMInputFields p = detail::string_to_enum(parameter);
  switch (p) {
    case detail::PCMInputFields::patch_level:
      host_input_.patch_level = value;
      break;
    case detail::PCMInputFields::der_order:
      host_input_.der_order = value;
      break;
    default:
      std::ostringstream errmsg;
      errmsg << "Unknown parameter " << parameter << std::endl;
      errmsg << " See http://pcmsolver.readthedocs.org/en/latest/users/input.html"
             << std::endl;
      PCMSOLVER_ERROR(errmsg.str());
  }
}

void pcmsolver_set_double_option(pcmsolver_context_t * context,
                                 const char * parameter,
                                 double value) {
  AS_TYPE(pcm::Meddle, context)->setInputOption(std::string(parameter), value);
}
void pcm::Meddle::setInputOption(std::string parameter, double value) {
  detail::PCMInputFields p = detail::string_to_enum(parameter);
  switch (p) {
    case detail::PCMInputFields::coarsity:
      host_input_.coarsity = value;
      break;
    case detail::PCMInputFields::area:
      host_input_.area = value;
      break;
    case detail::PCMInputFields::min_distance:
      host_input_.min_distance = value;
      break;
    case detail::PCMInputFields::min_radius:
      host_input_.min_radius = value;
      break;
    case detail::PCMInputFields::correction:
      host_input_.correction = value;
      break;
    case detail::PCMInputFields::probe_radius:
      host_input_.probe_radius = value;
      break;
    case detail::PCMInputFields::outside_epsilon:
      host_input_.outside_epsilon = value;
      break;
    default:
      std::ostringstream errmsg;
      errmsg << "Unknown parameter " << parameter << std::endl;
      errmsg << " See http://pcmsolver.readthedocs.org/en/latest/users/input.html"
             << std::endl;
      PCMSOLVER_ERROR(errmsg.str());
  }
}

void pcmsolver_set_string_option(pcmsolver_context_t * context,
                                 const char * parameter,
                                 const char * value) {
  AS_TYPE(pcm::Meddle, context)
      ->setInputOption(std::string(parameter), std::string(value));
}
void pcm::Meddle::setInputOption(std::string parameter, std::string value) {
  detail::PCMInputFields p = detail::string_to_enum(parameter);
  switch (p) {
    case detail::PCMInputFields::cavity_type:
      strcpy(host_input_.cavity_type, value.c_str());
      break;
    case detail::PCMInputFields::radii_set:
      strcpy(host_input_.radii_set, value.c_str());
      break;
    case detail::PCMInputFields::restart_name:
      strcpy(host_input_.restart_name, value.c_str());
      break;
    case detail::PCMInputFields::solver_type:
      strcpy(host_input_.solver_type, value.c_str());
      break;
    case detail::PCMInputFields::solvent:
      strcpy(host_input_.solvent, value.c_str());
      break;
    case detail::PCMInputFields::equation_type:
      strcpy(host_input_.equation_type, value.c_str());
      break;
    case detail::PCMInputFields::inside_type:
      strcpy(host_input_.inside_type, value.c_str());
      break;
    case detail::PCMInputFields::outside_type:
      strcpy(host_input_.outside_type, value.c_str());
      break;
    default:
      std::ostringstream errmsg;
      errmsg << "Unknown parameter " << parameter << std::endl;
      errmsg << " See http://pcmsolver.readthedocs.org/en/latest/users/input.html"
             << std::endl;
      PCMSOLVER_ERROR(errmsg.str());
  }
}

void pcmsolver_print(pcmsolver_context_t * context) {
  AS_CTYPE(pcm::Meddle, context)->printInfo();
}
void pcm::Meddle::printInfo() const { hostWriter_(infoStream_); }

void pcmsolver_citation(HostWriter writer) { writer(citation_message().c_str()); }

std::string pcm::Meddle::printCitation() const { return citation_message(); }

bool pcmsolver_is_compatible_library(void) {
  unsigned int major = (pcm::pcmsolver_get_version() >> 16);
  return (major == PROJECT_VERSION_MAJOR);
}
unsigned int pcm::pcmsolver_get_version(void) { return PCMSOLVER_VERSION; }

namespace pcm {
Molecule Meddle::molecule() const { return input_.molecule(); }

Eigen::Matrix3Xd Meddle::getCenters() const { return cavity_->elementCenter(); }

void Meddle::initInput(int nr_nuclei,
                       double charges[],
                       double coordinates[],
                       int symmetry_info[],
                       bool deferred_init) {
  // Position and charges of atomic centers
  Eigen::VectorXd chg = Eigen::Map<Eigen::VectorXd>(charges, nr_nuclei, 1);
  Eigen::Matrix3Xd centers = Eigen::Map<Eigen::Matrix3Xd>(coordinates, 3, nr_nuclei);

  if (input_.mode() != "EXPLICIT") {
    auto pg = Symmetry(
        symmetry_info[0], symmetry_info[1], symmetry_info[2], symmetry_info[3]);
    input_.molecule(
        detail::initMolecule(input_, pg, nr_nuclei, chg, centers, deferred_init));
  }
}

void Meddle::initCavity() {
  cavity_ = cavity::bootstrapFactory().create(input_.cavityParams().cavityType,
                                              input_.cavityParams());
  size_ = std::make_tuple(cavity_->size(), cavity_->irreducible_size());
  cavity_->saveCavity();

  infoStream_ << "========== Cavity " << std::endl;
  infoStream_ << "Atomic radii set: " << input_.radiiSetName() << std::endl;
  infoStream_ << *cavity_;
}

void Meddle::initStaticSolver() {
  IGreensFunction * gf_i = green::bootstrapFactory().create(
      input_.insideGreenParams().greensFunctionType, input_.insideGreenParams());
  IGreensFunction * gf_o = green::bootstrapFactory().create(
      input_.outsideStaticGreenParams().greensFunctionType,
      input_.outsideStaticGreenParams());

  K_0_ = solver::bootstrapFactory().create(input_.solverParams().solverType,
                                           input_.solverParams());

  IBoundaryIntegralOperator * biop = bi_operators::bootstrapFactory().create(
      input_.integratorParams().integratorType, input_.integratorParams());
  K_0_->buildSystemMatrix(*cavity_, *gf_i, *gf_o, *biop);
  delete biop;

  // Perform Gauss' theorem check for nuclear charges
  if (!hasFQ_)
    GaussCheck();

  infoStream_ << "========== Static solver " << std::endl;
  infoStream_ << *K_0_ << std::endl;
  mediumInfo(gf_i, gf_o);
  delete gf_o;
  delete gf_i;
}

void Meddle::initDynamicSolver() {
  IGreensFunction * gf_i = green::bootstrapFactory().create(
      input_.insideGreenParams().greensFunctionType, input_.insideGreenParams());
  IGreensFunction * gf_o = green::bootstrapFactory().create(
      input_.outsideDynamicGreenParams().greensFunctionType,
      input_.outsideDynamicGreenParams());

  K_d_ = solver::bootstrapFactory().create(input_.solverParams().solverType,
                                           input_.solverParams());

  IBoundaryIntegralOperator * biop = bi_operators::bootstrapFactory().create(
      input_.integratorParams().integratorType, input_.integratorParams());
  K_d_->buildSystemMatrix(*cavity_, *gf_i, *gf_o, *biop);
  hasDynamic_ = true;
  delete biop;

  infoStream_ << "========== Dynamic solver " << std::endl;
  infoStream_ << *K_d_ << std::endl;
  mediumInfo(gf_i, gf_o);
  delete gf_o;
  delete gf_i;
}

void Meddle::initMMFQ() {
  FQ_ = new mmfq::FQOhno(input_.fragments(), input_.isNonPolarizable());
  size_ = std::make_tuple(input_.fragments().sites.cols(),
                          input_.fragments().sites.cols());
  hasFQ_ = true;

  infoStream_ << "========== MMFQ solver " << std::endl;
  infoStream_ << *FQ_ << std::endl;
}

void Meddle::mediumInfo(IGreensFunction * gf_i, IGreensFunction * gf_o) {
  using utils::Solvent;
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

void Meddle::GaussCheck() const {
  Eigen::VectorXd nuclear_mep = computeMEP(input_.molecule(), cavity_->elements());
  Eigen::VectorXd nuclear_asc = K_0_->computeCharge(nuclear_mep);
  double total_nuclear_asc = nuclear_asc.sum() * cavity_->pointGroup().nrIrrep();
  double gauss_nuclear_asc = GaussEstimate(input_.molecule().charges(),
                                           input_.outsideStaticGreenParams().epsilon,
                                           input_.correction());
  double abs_rel_diff =
      std::abs((total_nuclear_asc - gauss_nuclear_asc) / gauss_nuclear_asc);
  std::stringstream tmp;
  if (!utils::isZero(abs_rel_diff, 1.0e-2)) {
    std::ostringstream errmsg;
    errmsg
        << "Absolute value of the relative difference between the Gauss' theorem ("
        << gauss_nuclear_asc << ") ";
    errmsg << "and computed (" << total_nuclear_asc << ") values ";
    errmsg << "of the total nuclear ASC higher than threshold (" << abs_rel_diff
           << ")." << std::endl;
    errmsg << "Consider changing the average area of the cavity finite elements."
           << std::endl;
    errmsg
        << "Please report this issue: https://github.com/PCMSolver/pcmsolver/issues"
        << std::endl;
    PCMSOLVER_ERROR(errmsg.str());
  }
}

namespace detail {
Molecule initMolecule(const Input & inp,
                      const Symmetry & pg,
                      int nuclei,
                      const Eigen::VectorXd & charges,
                      const Eigen::Matrix3Xd & centers,
                      bool deferred_init) {
  bool scaling = inp.scaling();
  std::string set = inp.radiiSet();
  std::vector<Atom> radiiSet;
  std::vector<Atom> atoms;
  // FIXME Code duplication in function initMolecule in interface/Input.cpp
  tie(ignore, radiiSet) = utils::bootstrapRadiiSet().create(set);
  std::vector<Sphere> spheres;
  atoms.reserve(nuclei);
  spheres.reserve(nuclei);
  for (int i = 0; i < charges.size(); ++i) {
    int index = int(charges(i)) - 1;
    atoms.push_back(radiiSet[index]);
    if (scaling)
      atoms[i].radiusScaling = 1.2;
    // Convert to Bohr and multiply by scaling factor (alpha)
    double radius = atoms[i].radius * angstromToBohr() * atoms[i].radiusScaling;
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
  Molecule molecule(nuclei, charges, masses, centers, atoms, spheres, pg);
  // Check that all atoms have a radius attached
  std::vector<Atom>::const_iterator res =
      std::find_if(atoms.begin(), atoms.end(), invalid);
  if (res != atoms.end() && !deferred_init) {
    std::ostringstream errmsg;
    errmsg << "In the molecule:\n" << molecule << std::endl;
    errmsg << "Some atoms do not have a radius attached." << std::endl;
    errmsg << "Please specify a radius for all atoms!" << std::endl;
    errmsg << " See http://pcmsolver.readthedocs.org/en/latest/users/input.html"
           << std::endl;
    PCMSOLVER_ERROR(errmsg.str());
  }
  return molecule;
}

void initSpheresAtoms(const Input & inp,
                      const Eigen::Matrix3Xd & sphereCenter_,
                      std::vector<Sphere> & spheres_) {
  // Loop over the atomsInput array to get which atoms will have a user-given radius
  for (size_t i = 0; i < inp.atoms().size(); ++i) {
    size_t index =
        inp.atoms(i) - 1; // -1 to go from human readable to machine readable
    // Put the new Sphere in place of the implicit-generated one
    spheres_[index] = Sphere(sphereCenter_.col(index), inp.radii(i));
  }
}

void print(const PCMInput & inp) {
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
} // namespace detail
} // namespace pcm
