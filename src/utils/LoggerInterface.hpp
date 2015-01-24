#ifndef LOGGERINTERFACE_HPP
#define LOGGERINTERFACE_HPP

#ifdef HAS_CXX11

#include "Logger.hpp"
#include "Timer.hpp"

static logging::logger<logging::FileLogPolicy> loggerInstance("pcmsolver.execution.log");

#define LOG loggerInstance.print<logging::printLevel::coarse>
#define LOG_FINE loggerInstance.print<logging::printLevel::fine>
#define LOG_ALL loggerInstance.print<logging::printLevel::everything>

#else // HAS_CXX11

#define LOG(...)
#define LOG_FINE(...)
#define LOG_ALL(...)

#endif // HAS_CXX11
 
#endif // LOGGERINTERFACE_HPP
