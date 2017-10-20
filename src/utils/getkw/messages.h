/** @file messages.h
 *
 * @brief Collection of assertions and a standard error/warn/info/debug
 * message interface.
 *
 * Written by Jonas Juselius <jonas.juselius@chem.uit.no>
 * CTCC, University of Tromsø, July 2009
 *
 */
#ifndef MESSAGES_H
#define MESSAGES_H

#include <cassert>
#include <iostream>
#include <sstream>
#include <string>

#include "GetkwError.h"

namespace GetkwMessageStream {
struct msg {
  static int DebugLevel;
  static std::ostream * out;
  static void setOutputStream(std::ostream & o) { out = &o; }
  static void setDebugLevel(int i) { DebugLevel = i; }
};
}

#define STR_DEBUG(S, X)                                                             \
  {                                                                                 \
    std::ostringstream _str;                                                        \
    _str << "Debug: " << __func__ << ",  line " << __LINE__ << " in  " << __FILE__  \
         << ": " << X << endl;                                                      \
    S = _str.str();                                                                 \
  }
#define STR_INFO(S, X)                                                              \
  {                                                                                 \
    std::ostringstream _str;                                                        \
    _str << "Info: " << __func__ << ",  line " << __LINE__ << " in  " << __FILE__   \
         << ": " << X << endl;                                                      \
    S = _str.str();                                                                 \
  }
#define STR_WARN(S, X)                                                              \
  {                                                                                 \
    std::ostringstream _str;                                                        \
    _str << "Warning: " << __func__ << ",  line " << __LINE__ << " in  "            \
         << __FILE__ << ": " << X << endl;                                          \
    S = _str.str();                                                                 \
  }
#define STR_ERROR(S, X)                                                             \
  {                                                                                 \
    std::ostringstream _str;                                                        \
    _str << "Error: " << __func__ << ",  line " << __LINE__ << " in  " << __FILE__  \
         << ": " << X << endl;                                                      \
    S = _str.str();                                                                 \
  }

#ifdef NDEBUG
#define PRINT(STR)
#define SPRINT(STR)
#else
#define PRINT(STR) *GetkwMessageStream::msg::out << STR;
#define SPRINT(_out, STR) _out << STR;
#endif

#define debug(level, STR)                                                           \
  if (level < GetkwMessageStream::msg::DebugLevel)                                  \
    *GetkwMessageStream::msg::out << STR;
#define sdebug(level, _out, STR)                                                    \
  if (level < GetkwMessageStream::msg::DebugLevel)                                  \
    _out << STR;

#define SET_DEBUG_LEVEL(a) GetkwMessageStream::msg::setDebugLevel(a);
#define SET_MESSAGE_STREAM(s) GetkwMessageStream::msg::setOutputStrem(s);
#define DEBUG_LEVEL GetkwMessageStream::msg::DebugLevel

#define MSG_DEBUG(X)                                                                \
  {                                                                                 \
    *GetkwMessageStream::msg::out << "Debug: " << __func__ << "(), line "           \
                                  << __LINE__ << "in " << __FILE__ << ": " << X     \
                                  << endl;                                          \
  }
#define MSG_INFO(X)                                                                 \
  {                                                                                 \
    *GetkwMessageStream::msg::out << "Info: " << __FILE__ << ": " << __func__       \
                                  << "(), line " << __LINE__ << ": " << X << endl;  \
  }
#define MSG_NOTE(X)                                                                 \
  {                                                                                 \
    *GetkwMessageStream::msg::out << "Note: " << __FILE__ << ": " << __func__       \
                                  << "(), line " << __LINE__ << ": " << X << endl;  \
  }
#define MSG_WARN(X)                                                                 \
  {                                                                                 \
    *GetkwMessageStream::msg::out << "Warning: " << __func__ << "(), line "         \
                                  << __LINE__ << ": " << X << endl;                 \
  }
#define MSG_ERROR(X)                                                                \
  {                                                                                 \
    *GetkwMessageStream::msg::out << "Error: " << __func__ << "(), line "           \
                                  << __LINE__ << ": " << X << endl;                 \
  }
#define MSG_FATAL(X)                                                                \
  {                                                                                 \
    *GetkwMessageStream::msg::out << "Error: " << __FILE__ << ": " << __func__      \
                                  << "(), line " << __LINE__ << ": " << X << endl;  \
    abort();                                                                        \
  }

#define MSG_INVALID_ARG(X)                                                          \
  {                                                                                 \
    *GetkwMessageStream::msg::out << "Error, invalid argument passed: " << __func__ \
                                  << "(), line " << __LINE__ << ": " << X << endl;  \
  }
#define INVALID_ARG_ABORT                                                           \
  {                                                                                 \
    *GetkwMessageStream::msg::out << "Error, invalid argument passed: " << __func__ \
                                  << "(), line " << __LINE__ << endl;               \
    abort();                                                                        \
  }
#define NOT_REACHED_ABORT                                                           \
  {                                                                                 \
    *GetkwMessageStream::msg::out << "Error, should not be reached: " << __func__   \
                                  << "(), line " << __LINE__ << endl;               \
    abort();                                                                        \
  }
#define INTERNAL_INCONSISTENCY                                                      \
  {                                                                                 \
    *GetkwMessageStream::msg::out                                                   \
        << "Internal inconsistency! You have found a bug: " << __func__             \
        << "(), line " << __LINE__ << endl;                                         \
    abort();                                                                        \
  }

#define NEEDS_TESTING                                                               \
  {                                                                                 \
    static bool __once = true;                                                      \
    if (__once) {                                                                   \
      __once = false;                                                               \
      *GetkwMessageStream::msg::out << "NEEDS TESTING: " << __FILE__ << ", "        \
                                    << __func__ << "(), line " << __LINE__ << endl; \
    }                                                                               \
  }

#define ASSERT_FILE(A, B)                                                           \
  {                                                                                 \
    if (A == NULL) {                                                                \
      *GetkwMessageStream::msg::out << "Error: " << __func__ << "(), line "         \
                                    << __LINE__ << ": No such file, " << B << endl; \
      abort();                                                                      \
    }                                                                               \
  }

#define NOT_IMPLEMENTED_ABORT                                                       \
  {                                                                                 \
    *GetkwMessageStream::msg::out << "Error: Not implemented, " << __FILE__ ", "    \
                                  << __func__ << "(), line " << __LINE__ << endl;   \
    abort();                                                                        \
  }

#define NOTE(X)                                                                     \
  {                                                                                 \
    static bool __once = true;                                                      \
    if (__once) {                                                                   \
      __once = false;                                                               \
      *GetkwMessageStream::msg::out << "NOTE: " << __FILE__ << ", " << __func__     \
                                    << "(), line " << __LINE__ << ": " << X         \
                                    << endl;                                        \
    }                                                                               \
  }

#define NEEDS_FIX(X)                                                                \
  {                                                                                 \
    static bool __once = true;                                                      \
    if (__once) {                                                                   \
      __once = false;                                                               \
      *GetkwMessageStream::msg::out << "NEEDS FIX: " << __FILE__ << ", "            \
                                    << __func__ << "(), line " << __LINE__ << ": "  \
                                    << X << endl;                                   \
    }                                                                               \
  }

#define WRONG(X)                                                                    \
  {                                                                                 \
    *GetkwMessageStream::msg::out << "WRONG: " << __FILE__ << ", " << __func__      \
                                  << "(), line " << __LINE__ << ": " << X << endl;  \
    abort();                                                                        \
  }

#define STR_DEBUG(S, X)                                                             \
  {                                                                                 \
    std::ostringstream _str;                                                        \
    _str << "Debug: " << __func__ << ",  line " << __LINE__ << " in  " << __FILE__  \
         << ": " << X << endl;                                                      \
    S = _str.str();                                                                 \
  }
#define STR_INFO(S, X)                                                              \
  {                                                                                 \
    std::ostringstream _str;                                                        \
    _str << "Info: " << __func__ << ",  line " << __LINE__ << " in  " << __FILE__   \
         << ": " << X << endl;                                                      \
    S = _str.str();                                                                 \
  }
#define STR_WARN(S, X)                                                              \
  {                                                                                 \
    std::ostringstream _str;                                                        \
    _str << "Warning: " << __func__ << ",  line " << __LINE__ << " in  "            \
         << __FILE__ << ": " << X << endl;                                          \
    S = _str.str();                                                                 \
  }
#define STR_ERROR(S, X)                                                             \
  {                                                                                 \
    std::ostringstream _str;                                                        \
    _str << "Error: " << __func__ << ",  line " << __LINE__ << " in  " << __FILE__  \
         << ": " << X << endl;                                                      \
    S = _str.str();                                                                 \
  }

/* The quiet versions...
 #define SET_DEBUG_LEVEL(a)
 #define SET_MESSAGE_STREAM(s)
 #define DEBUG_LEVEL

 #define MSG_DEBUG(X)
 #define MSG_INFO(X)
 #define MSG_WARN(X)
 #define MSG_ERROR(X)
 #define MSG_FATAL(X) abort();


 #define MSG_INVALID_ARG
 #define INVALID_ARG_ABORT abort();
 #define NOT_REACHED_ABORT abort();

 #define NEEDS_TESTING

 #define ASSERT_FILE(A,B)
 #define NOT_IMPLEMENTED_ABORT abort();
 */
#endif
