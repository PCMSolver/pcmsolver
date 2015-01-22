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

#include "Timer.hpp"

#include "Config.hpp"

#include <boost/foreach.hpp>

std::ostream & Timer::printObject(std::ostream & os) const
{
    timingsPair t_pair;
    BOOST_FOREACH(t_pair, timings_) {
        os << t_pair.first << " : " << boost::timer::format(t_pair.second);
    }
    return os;
}

void Timer::insertTimer(const timersPair & checkpoint)
{
    timers_.insert(checkpoint);
}

void Timer::eraseTimer(const std::string & checkpoint_name)
{
    timers_.erase(checkpoint_name);
}

void Timer::insertTiming(const std::string & checkpoint_name)
{
    // Find timer associated with given checkpoint_name
    timersMap::iterator checkpoint_timer = timers_.find(checkpoint_name);
    // Get elapsed time, create a timingsPair and insert it into timings_ map
    boost::timer::cpu_times const checkpoint_times((checkpoint_timer->second).elapsed());
    timingsPair checkpoint = timingsMap::value_type(checkpoint_name, checkpoint_times);
    timings_.insert(checkpoint);
}

void printTimings(const std::string & fname)
{
    namespace fs = boost::filesystem;
    fs::path file(fname);
    std::ofstream timing_report;
    if (fs::exists(file)) {
        fs::remove(fname);
        timing_report.open(fname, std::ios::out);
        timing_report << "            PCMSolver API timing results            " << std::endl;
        timing_report << "----------------------------------------------------" << std::endl;
    }
    timing_report << Timer::TheTimer() << std::endl;
    timing_report.close();
}

void timerON(const std::string & checkpoint_name)
{
    boost::timer::cpu_timer checkpoint_timer;
    timersPair checkpoint = timersMap::value_type(checkpoint_name, checkpoint_timer);
    Timer::TheTimer().insertTimer(checkpoint);
}

void timerOFF(const std::string & checkpoint_name)
{
    Timer::TheTimer().insertTiming(checkpoint_name);
    Timer::TheTimer().eraseTimer(checkpoint_name);
    // If all timers are turned OFF write results to file
    if (Timer::TheTimer().activeTimers() == 0) {
        printTimings("pcmsolver.timer.dat");
    }
}
