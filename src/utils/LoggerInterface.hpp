#ifndef LOGGERINTERFACE_HPP
#define LOGGERINTERFACE_HPP

#ifdef HAS_CXX11

#include "Logger.hpp"
#include "Timer.hpp"

static logging::logger<logging::FileLogPolicy> loggerInstance("pcmsolver.execution.log");

#define LOG loggerInstance.print<logging::printLevel::coarse>
#define LOG_FINE loggerInstance.print<logging::printLevel::fine>
#define LOG_ALL loggerInstance.print<logging::printLevel::everything>
#define LOG_TIME loggerInstance.print<logging::printLevel::timings>(Timer::TheTimer());

#else // HAS_CXX11

#define LOG(...)
#define LOG_FINE(...)
#define LOG_ALL(...)
#define LOG_TIME(...)

#endif // HAS_CXX11
 
#endif // LOGGERINTERFACE_HPP
