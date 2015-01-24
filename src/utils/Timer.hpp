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

#include <map>
#include <sstream>
#include <string>
#include <utility>

#include "Config.hpp"

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
    /// Number of active timers
    size_t active_;
    /// Final timings report
    std::ostringstream report_;

    Timer();
    Timer(const Timer & other);
    Timer& operator=(const Timer & other);
    ~Timer();
    std::ostream & printObject(std::ostream & os) const;
    /*! Returns names of active timers */
    std::string activeTimers(); 
public:
    static Timer& TheTimer() {
        static Timer obj;
        return obj;
    }
    friend std::ostream & operator<<(std::ostream & os, const Timer & timer) {
        return timer.printObject(os);
    }

    /*! \brief Inserts a checkpoint-timer pair in the timers_ map
     *  \param[in] checkpoint timer to be stored 
     */
    void insertTimer(const timersPair & checkpoint);

    /*! \brief Erases a checkpoint-timer pair in the timers_ map
     *  \param[in] checkpoint_name timer to be erased 
     */
    void eraseTimer(const std::string & checkpoint_name);

    /*! \brief Inserts a checkpoint-timing pair in the timings_ map
     *  \param[in] checkpoint_name timing to be stored 
     */
    void insertTiming(const std::string & checkpoint_name);

    /*! \brief Appends timing to report_ stream
     *  \param[in] checkpoint_name timings to be reported
     */
    void reportTiming(const std::string & checkpoint_name);

    /*! Returns number of active timers */
    int active() { return active_; }
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
 */
void timerOFF(const std::string & checkpoint_name);

#endif // TIMER_HPP
