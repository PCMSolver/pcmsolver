/*
 * \file messages.cpp
 *
 *  \date Jul 24, 2009
 *  \author Jonas Juselius <jonas.juselius@uit.no> \n
 *          CTCC, University of Tromsø
 *
 * \breif
 */

#include "messages.h"

int GetkwMessageStream::msg::DebugLevel = 0;
std::ostream * GetkwMessageStream::msg::out = &std::cout;
