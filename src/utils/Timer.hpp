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

#include <boost/foreach.hpp>
#include <boost/timer/timer.hpp>

namespace timer {
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

    std::ostream & printObject(std::ostream & os) const {
        os << "            PCMSolver API timing results            " << std::endl;
        os << "----------------------------------------------------" << std::endl;
        timingsPair t_pair;
        BOOST_FOREACH(t_pair, timings_) {
            os << t_pair.first << " : " << boost::timer::format(t_pair.second);
        }
        if (active_ != 0) {
            os << " These timers were not shut down:" << std::endl;
            timersPair t_pair;
            BOOST_FOREACH(t_pair, timers_) {
                os << " - " << t_pair.first << std::endl;
            }
            os << " Reported timings might be unreliable!" << std::endl;
        }
        os << "----------------------------------------------------" << std::endl;
        return os;
    }
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
    void insertTimer(const timersPair & checkpoint) {
        timers_.insert(checkpoint);
        ++active_;
    }

    /*! \brief Erases a checkpoint-timer pair in the timers_ map
     *  \param[in] checkpoint_name timer to be erased
     */
    void eraseTimer(const std::string & checkpoint_name) {
        timers_.erase(checkpoint_name);
        --active_;
    }

    /*! \brief Inserts a checkpoint-timing pair in the timings_ map
     *  \param[in] checkpoint_name timing to be stored
     */
    void insertTiming(const std::string & checkpoint_name) {
        // Find timer associated with given checkpoint_name
        timersMap::iterator checkpoint_timer = timers_.find(checkpoint_name);
        // Get elapsed time, create a timingsPair and insert it into timings_ map
        boost::timer::cpu_times const checkpoint_times((checkpoint_timer->second).elapsed());
        timingsPair checkpoint = timingsMap::value_type(checkpoint_name, checkpoint_times);
        timings_.insert(checkpoint);
    }

    /*! Returns number of active timers */
    int active() const { return active_; }
};

/*! \fn inline void timerON(const std::string & checkpoint_name)
 *  \param[in] checkpoint_name name of the checkpoint to be timed
 *  \brief Starts a timer with the given name
 *
 *  A timer is added to the timers_ map of the Timer object.
 */
inline void timerON(const std::string & checkpoint_name)
{
    boost::timer::cpu_timer checkpoint_timer;
    timersPair checkpoint = timersMap::value_type(checkpoint_name, checkpoint_timer);
    Timer::TheTimer().insertTimer(checkpoint);
}

/*! \fn inline void timerOFF(const std::string & checkpoint_name)
 *  \param[in] checkpoint_name name of the checkpoint to be timed
 *  \brief Stops a timer with the given name
 *
 *  The timing results associated with the given timer
 *  are added to the timings_ map of the Timer object,
 *  the timer is then removed from the timers_ map.
 */
inline void timerOFF(const std::string & checkpoint_name)
{
    Timer::TheTimer().insertTiming(checkpoint_name);
    Timer::TheTimer().eraseTimer(checkpoint_name);
}

/*! \fn inline void timerDONE(const std::string & fname)
 *  \param[in] fname name for the timing report file
 */
inline void timerDONE(const std::string & fname)
{
    std::ofstream timing_report;
    timing_report.open(fname.c_str(), std::ios::out);
    timing_report << Timer::TheTimer() << std::endl;
    timing_report.close();
}
} // namespace timer

#endif // TIMER_HPP
