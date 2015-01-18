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

/*! \file Timer.hpp
 *  \class Timer
 *  \brief A wrapper around Boost.Timer
 *  \author Roberto Di Remigio
 *  \date 2014
 */

class Timer
{
private:
    typedef std::map<std::string, boost::timer::cpu_times> timingsMap;
    typedef std::pair<std::string, boost::timer::cpu_times> timingsPair;
private:
    timersMap timers_;
    timingsMap timings_;

    Timer() : timers_(), timings_() {}
    Timer(const Timer & other);
    Timer& operator=(const Timer & other);
    ~Timer() {}
    std::ostream & printObject(std::ostream & os) const {
        typedef timingsMap::const_iterator it_type;
        for (it_type iter = timings_.begin(); iter != timings_.end(); ++iter) {
            os << iter->first << " : " << boost::timer::format(iter->second);
        }
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
    timersMap & timers() { return timers_; }
    void appendTimings(const std::string & checkpoint_name,
                       const boost::timer::cpu_times & checkpoint_times) {
        timingsPair checkpoint = timingsMap::value_type(checkpoint_name, checkpoint_times);
        timings_.insert(checkpoint);
    }
};

inline void printTimings() {
	namespace fs = boost::filesystem;
	fs::path file("pcmsolver.timer.dat");
        std::ofstream timing_report;
	if (fs::exists(file)) {
        	timing_report.open("pcmsolver.timer.dat", std::ios::out | std::ios::app);
	} else {
        	timing_report.open("pcmsolver.timer.dat", std::ios::out);
        	timing_report << "            PCMSolver API timing results            " << std::endl;
	        timing_report << "----------------------------------------------------" << std::endl;
	}
        timing_report << Timer::TheTimer() << std::endl;
        timing_report.close();
}
    
inline void timerON(const std::string & checkpoint_name) {
        boost::timer::cpu_timer checkpoint_timer;
        timersPair checkpoint = timersMap::value_type(checkpoint_name, checkpoint_timer);
	Timer::TheTimer().timers().insert(checkpoint);
}
    
inline void timerOFF(const std::string & checkpoint_name) {
        timersMap::iterator checkpoint_timer = Timer::TheTimer().timers().find(checkpoint_name);
        boost::timer::cpu_times const checkpoint_times((checkpoint_timer->second).elapsed());
	Timer::TheTimer().appendTimings(checkpoint_name, checkpoint_times);
	// Remove checkpoint_timer from the timers_ map
	Timer::TheTimer().timers().erase(checkpoint_name);
	// If all timers are turned OFF write results to file
	if (Timer::TheTimer().timers().size() == 0) {
		printTimings();
	}
}

#endif // TIMER_HPP
