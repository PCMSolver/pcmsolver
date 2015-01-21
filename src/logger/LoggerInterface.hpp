#ifndef LOGGERINTERFACE_HPP
#define LOGGERINTERFACE_HPP

#ifdef HAS_CXX11_SUPPORT

#include "Logger.hpp"

static logging::logger<logging::FileLogPolicy> loggerInstance("execution.log");

#ifdef LOGGING_LEVEL_1
#define LOG loggerInstance.print<logging::severityType::debug>
#define LOG_ERR loggerInstance.print<logging::severityType::error>
#define LOG_WARN loggerInstance.print<logging::severityType::warning>
#else
#define LOG(...)
#define LOG_ERR(...)
#define LOG_WARN(...)
#endif

#ifdef LOGGING_LEVEL_2
#define ELOG loggerInstance.print<logging::severityType::debug>
#define ELOG_ERR loggerInstance.print<logging::severityType::error>
#define ELOG_WARN loggerInstance.print<logging::severityType::warning>
#else
#define ELOG(...)
#define ELOG_ERR(...)
#define ELOG_WARN(...)
#endif

#else // HAS_CXX11_SUPPORT

#ifdef LOGGING_LEVEL_1
#define LOG(...)
#define LOG_ERR(...)
#define LOG_WARN(...)
#endif

#ifdef LOGGING_LEVEL_2
#define ELOG(...)
#define ELOG_ERR(...)
#define ELOG_WARN(...)
#endif

#endif // HAS_CXX11_SUPPORT
 
#endif // LOGGERINTERFACE_HPP
