#ifndef TIMERINTERFACE_HPP
#define TIMERINTERFACE_HPP

/*  To time a code snippet:
 *  \code{.cpp}
 *  TIMER_ON("code-snippet");
 *  // code-snippet
 *  TIMER_OFF("code-snippet");
 *  \endcode
 *  The timings are printed out by a call to the TIMER_DONE function.
 */
#ifdef ENABLE_TIMER

#include "Timer.hpp"

#define TIMER_ON timer::timerON
#define TIMER_OFF timer::timerOFF
#define TIMER_DONE timer::timerDONE

#else // ENABLE_TIMER

#define TIMER_ON(...)
#define TIMER_OFF(...)
#define TIMER_DONE(...)

#endif // ENABLE_TIMER

#endif // TIMERINTERFACE_HPP
