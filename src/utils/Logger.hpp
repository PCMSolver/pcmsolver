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

#include <ctime>
#include <mutex>
#include <sstream>
#include <string>

#include "LoggerImpl.hpp"

namespace logging {
std::string getTime();

template <typename logPolicy> class logger {
private:
  printLevel globalPrintLevel_;
  std::stringstream logStream_;
  logPolicy * policy_;
  std::mutex writeMutex_;

  /*! @name Core printing functionality
   *
   *  A variadic template is used, we specify the version
   *  with a variable number of parameters/types and the version
   *  with no parameters/types.
   *  The variadic template is called recursively.
   *  The type for the first parameter is resolved and streamed
   *  to logStream_. When all the parameters have been streamed
   *  the version with no arguments is called.
   */
  /// @{
  /*! */
  void printImpl() {
    policy_->write(logStream_.str());
    logStream_.str("");
  }
  template <typename First, typename... Rest>
  void printImpl(First parm1, Rest... parm) {
    logStream_.precision(std::numeric_limits<double>::digits10);
    logStream_ << parm1 << std::endl;
    printImpl(parm...);
  }
  /// @}
public:
  /*! Constructor
   *  \param[in] name name for the log file
   *  The build parameters are logged first
   */
  logger(const std::string & name, printLevel print = coarse)
      : globalPrintLevel_(print), policy_(new logPolicy) {
    if (!policy_) {
      PCMSOLVER_ERROR("LOGGER: Unable to create the logger instance");
    }
    policy_->open_ostream(name);
    // Write the logfile header
    logStream_ << "\t\tPCMSolver execution log\n"
               << buildInfo() << "\n\t\tLog started : " << getTime() << std::endl;
  }
  /// Destructor
  ~logger() {
    if (policy_) {
      policy_->close_ostream();
      delete policy_;
    }
  }

  void globalPrintLevel(int printLvl) { globalPrintLevel_ = printLvl; }

  /// User interface for the logger class
  template <printLevel printLvl, typename... Args> void print(Args... args) {
    if (globalPrintLevel_ >= printLvl) {
      writeMutex_.lock();
      printImpl(args...);
      writeMutex_.unlock();
    }
  }
};

/*! \brief Returns date and time */
inline std::string getTime() {
  std::string time_str;
  time_t raw_time;

  std::time(&raw_time);
  time_str = std::ctime(&raw_time);

  // Without the newline character
  return time_str;
}

} // close namespace logging
