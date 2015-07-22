#ifndef LOGGERINTERFACE_HPP
#define LOGGERINTERFACE_HPP

#ifdef ENABLE_LOGGER

#include "Logger.hpp"
#include "Timer.hpp"

static logging::logger<logging::FileLogPolicy> loggerInstance("pcmsolver.execution.log");

#define LOG loggerInstance.print<logging::printLevel::coarse>
#define LOG_FINE loggerInstance.print<logging::printLevel::fine>
#define LOG_ALL loggerInstance.print<logging::printLevel::everything>
#define LOG_TIME loggerInstance.print<logging::printLevel::timings>(timer::Timer::TheTimer());

#else // ENABLE_LOGGER

#define LOG(...)
#define LOG_FINE(...)
#define LOG_ALL(...)
#define LOG_TIME

#endif // ENABLE_LOGGER

#endif // LOGGERINTERFACE_HPP
