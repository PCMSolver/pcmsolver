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

#include <fstream>
#include <memory>
#include <string>

namespace logging {
enum printLevel { timings, coarse, fine, everything };

/*! \class ILogPolicy
 *  \brief ABC for logging policy classes
 */
class ILogPolicy {
public:
  /*! \brief Opens an output stream with the given name
   *  \param[in] name name of the stream
   */
  virtual void open_ostream(const std::string & name) = 0;
  /*! \brief Closes an output stream with the given name
   */
  virtual void close_ostream() = 0;
  /*! \brief Writes to stream
   *  \param[in] msg message to be written to stream
   */
  virtual void write(const std::string & msg) = 0;
};

/*! \class FileLogPolicy
 *  \brief Implementation which allows to write into a file
 */
class FileLogPolicy : public ILogPolicy {
private:
  std::unique_ptr<std::ofstream> outStream_;

public:
  /// Constructor
  FileLogPolicy() : outStream_(new std::ofstream) {}
  /// Destructor
  virtual ~FileLogPolicy() {
    if (outStream_) {
      close_ostream();
    }
  }
  /*! \brief Opens an output stream with the given name
   *  \param[in] name name of the stream
   */
  virtual void open_ostream(const std::string & name) {
    outStream_->open(name.c_str(), std::ios_base::binary | std::ios_base::out);
    if (!outStream_->is_open()) {
      PCMSOLVER_ERROR("LOGGER: Unable to open an output stream");
    }
  }
  /*! \brief Closes an output stream with the given name
   */
  virtual void close_ostream() {
    if (outStream_) {
      outStream_->close();
    }
  }
  /*! \brief Writes to stream
   *  \param[in] msg message to be written to stream
   */
  virtual void write(const std::string & msg) { (*outStream_) << msg << std::endl; }
};

} // close namespace logging
