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

#include <sstream>
#include <string>

#include "Config.hpp"

#include <boost/foreach.hpp>
    
Timer::Timer() : timers_(), timings_(), active_(0) 
{
    report_ << "            PCMSolver API timing results            " << std::endl; 
    report_ << "----------------------------------------------------" << std::endl;
}

Timer::~Timer() 
{
    if (active_ != 0) {
       report_ << activeTimers() << std::endl;
    }
}

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
    ++active_;
}

void Timer::eraseTimer(const std::string & checkpoint_name)
{
    timers_.erase(checkpoint_name);
    --active_;
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

void Timer::reportTiming(const std::string & checkpoint_name)
{
    timingsMap::iterator checkpoint = timings_.find(checkpoint_name);
    report_ << checkpoint->first << " : " << boost::timer::format(checkpoint->second) << std::endl; 
}
    
std::string Timer::activeTimers() 
{
    std::ostringstream tmp;
    timersPair t_pair;                                   
    tmp << " These timers were not shut down:" << std::endl;
    BOOST_FOREACH(t_pair, timers_) {                     
        tmp << " - " << t_pair.first << std::endl;       
    }                                                    
    tmp << " Reported timings might be unreliable!\n" << std::endl;
    // Remove last newline                               
    return tmp.str().substr(0, tmp.str().size() - 1);    
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
    // Write to report_ **after** timer has been erased
    Timer::TheTimer().reportTiming(checkpoint_name);
}
