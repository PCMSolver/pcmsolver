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
 *     PCMSolver API, see: <https://repo.ctcc.no/projects/pcmsolver>
 */
/* pcmsolver_copyright_end */

#ifndef TIMER_HPP
#define TIMER_HPP

#include <fstream>
#include <map>
#include <string>
#include <utility>

#include "Config.hpp"

#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/timer/timer.hpp>

typedef std::map<std::string, boost::timer::cpu_timer> timersMap;
typedef std::pair<std::string, boost::timer::cpu_timer> timersPair;

typedef std::map<std::string, boost::timer::cpu_times> timingsMap;
typedef std::pair<std::string, boost::timer::cpu_times> timingsPair;

/*! \file Timer.hpp
 *  \class Timer
 *  \brief A wrapper around Boost.Timer
 *  \author Roberto Di Remigio
 *  \date 2014
 *
 *  It is implemented as a Singleton. To time a code snippet:
 *  \code{.cpp}
 *  timerON("code-snippet");
 *  // code-snippet
 *  timerOFF("code-snippet");
 *  \endcode
 *  The timings are printed out when all the checkpoints created
 *  by a call to timerON are closed by a call to timerOFF.
 */

class Timer
{
private:
    /// Checkpoint-timer map
    timersMap timers_;
    /// Checkpoint-timing map
    timingsMap timings_;

    Timer() : timers_(), timings_() {}
    Timer(const Timer & other);
    Timer& operator=(const Timer & other);
    ~Timer() {}
    std::ostream & printObject(std::ostream & os) const;
public:
    static Timer& TheTimer() {
        static Timer obj;
        return obj;
    }
    friend std::ostream & operator<<(std::ostream & os, const Timer & timer) {
        return timer.printObject(os);
    }

    /*! \brief Inserts a checkpoint-timer pair in the timers_ map
     */
    void insertTimer(const timersPair & checkpoint);
    /*! \brief Erases a checkpoint-timer pair in the timers_ map
     */
    void eraseTimer(const std::string & checkpoint_name);

    /*! \brief Inserts a checkpoint-timing pair in the timings_ map
     */
    void insertTiming(const std::string & checkpoint_name);

    /*! \brief Returns number of active timers
     */
    int activeTimers() {
        return timers_.size();
    }
};

/*! \fn void timerON(const std::string & checkpoint_name)
 *  \param[in] checkpoint_name name of the checkpoint to be timed
 *  \brief Starts a timer with the given name
 *
 *  A timer is added to the timers_ map of the Timer object.
 */
void timerON(const std::string & checkpoint_name);

/*! \fn void timerOFF(const std::string & checkpoint_name)
 *  \param[in] checkpoint_name name of the checkpoint to be timed
 *  \brief Stops a timer with the given name
 *
 *  The timing results associated with the given timer
 *  are added to the timings_ map of the Timer object,
 *  the timer is then removed from the timers_ map.
 *  If no timers are left in the timers_ map, the final results
 *  are written to pcmsolver.timer.dat
 */
void timerOFF(const std::string & checkpoint_name);

/*! \fn printTimings(const std::string & fname)
 *  \param[in] fname timers report filename
 *  \brief Writes timing results to given filename
 *  \warning fname is removed if already existent
 *
 *  This function is invoked by timerOFF when there
 *  are no more active timers in the timers_ map.
 */
void printTimings(const std::string & fname);

#endif // TIMER_HPP
