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
#include <string>
#include <utility>

#include "Cxx11Workarounds.hpp"

#include <boost/container/flat_map.hpp>
#include <boost/foreach.hpp>

namespace timer {
typedef pcm::tuple<double, double> timing;
typedef boost::container::flat_map<std::string, timing> TimingsMap;
typedef std::pair<std::string, timing> TimingsPair;

timing get_timing();
double get_wall_time();
double get_cpu_time();

/*! \file Timer.hpp
 *  \class Timer
 *  \brief A wrapper around the basic C standard library timing facilities
 *  \author Roberto Di Remigio
 *  \date 2015
 */

class Timer {
private:
  /// Checkpoint-timing map
  TimingsMap timings_;

  std::ostream & printObject(std::ostream & os) const {
    os << "            PCMSolver API timing results            " << std::endl;
    os << "----------------------------------------------------" << std::endl;
    BOOST_FOREACH (TimingsPair t_pair, timings_) {
      os << "Checkpoint:  " << t_pair.first << std::endl;
      os << "   Wall time: " << pcm::get<0>(t_pair.second) << " ms" << std::endl;
      os << "   CPU time:  " << pcm::get<1>(t_pair.second) << " ms" << std::endl;
    }
    os << "----------------------------------------------------" << std::endl;
    return os;
  }

public:
  static Timer & TheTimer() {
    static Timer obj;
    return obj;
  }
  friend std::ostream & operator<<(std::ostream & os, const Timer & timer) {
    return timer.printObject(os);
  }

  /*! \brief Inserts a checkpoint-timing pair in the timings_ map
   *  \param[in] p checkpoint-timing to be stored
   */
  void insertTiming(const TimingsPair & p) { timings_.insert(p); }

  /*! \brief Calculates and registers elapsed time
   *  \param[in] chkpt_name timing checkpoint
   *  \param[in] t_stop stopping time
   */
  void registerElapsed(const std::string & chkpt_name, timing t_stop) {
    double wall_elapsed = pcm::get<0>(t_stop) - pcm::get<0>(timings_[chkpt_name]);
    double cpu_elapsed = pcm::get<1>(t_stop) - pcm::get<1>(timings_[chkpt_name]);
    timings_[chkpt_name] = pcm::make_tuple(wall_elapsed, cpu_elapsed);
  }
};

/*! \fn inline void timerON(const std::string & chkpt_name)
 *  \param[in] chkpt_name name of the checkpoint to be timed
 *  \brief Starts a timer with the given name
 *
 *  A timer is added to the timers_ map of the Timer object.
 */
inline void timerON(const std::string & chkpt_name) {
  Timer::TheTimer().insertTiming(std::make_pair(chkpt_name, get_timing()));
}

/*! \fn inline void timerOFF(const std::string & chkpt_name)
 *  \param[in] chkpt_name name of the checkpoint to be timed
 *  \brief Stops a timer with the given name
 *
 *  The timing results associated with the given timer
 *  are added to the timings_ map of the Timer object,
 *  the timer is then removed from the timers_ map.
 */
inline void timerOFF(const std::string & chkpt_name) {
  Timer::TheTimer().registerElapsed(chkpt_name, get_timing());
}

/*! \fn inline void timerDONE(const std::string & fname)
 *  \param[in] fname name for the timing report file
 */
inline void timerDONE(const std::string & fname) {
  std::ofstream timing_report;
  timing_report.open(fname.c_str(), std::ios::out);
  timing_report << Timer::TheTimer() << std::endl;
  timing_report.close();
}

/*! Returns wall and CPU times in milliseconds */
inline timing get_timing() {
  return pcm::make_tuple(get_wall_time() / 1000.0, get_cpu_time() / 1000.0);
}

// This code was taken from:
// http://stackoverflow.com/a/17440673/2528668
// It is most likely not going to work properly with OpenMP
#ifdef _WIN32
#include <Windows.h>
inline double get_wall_time() {
  LARGE_INTEGER time, freq;
  if (!QueryPerformanceFrequency(&freq)) {
    //  Handle error
    return 0;
  }
  if (!QueryPerformanceCounter(&time)) {
    //  Handle error
    return 0;
  }
  return (double)time.QuadPart / freq.QuadPart;
}
inline double get_cpu_time() {
  FILETIME a, b, c, d;
  if (GetProcessTimes(GetCurrentProcess(), &a, &b, &c, &d) != 0) {
    //  Returns total user time.
    //  Can be tweaked to include kernel times as well.
    return (double)(d.dwLowDateTime | ((unsigned long long)d.dwHighDateTime << 32)) *
           0.0000001;
  } else {
    //  Handle error
    return 0;
  }
}
#else /* _WIN32 */
#include <sys/time.h>
#include <time.h>

inline double get_wall_time() {
  struct timeval time;
  if (gettimeofday(&time, NULL)) {
    //  Handle error
    return 0;
  }
  return (double)time.tv_sec + (double)time.tv_usec * .000001;
}

inline double get_cpu_time() { return (double)clock() / CLOCKS_PER_SEC; }
#endif /* _WIN32 */
} // namespace timer
