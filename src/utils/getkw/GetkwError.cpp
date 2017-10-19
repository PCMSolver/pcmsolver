/*
 * \file GetkwError.cpp
 *
 *  \date Jul 2, 2009
 *  \author Jonas Juselius <jonas.juselius@uit.no> \n
 *          CTCC, University of Tromsø
 *
 * \breif
 */

#include "GetkwError.h"

bool GetkwError::strict = false;
bool GetkwError::verbose = true;

GetkwError::GetkwError() {}

GetkwError::GetkwError(const string & err) : msg(err) {
  if (verbose or strict) {
    cout << "Error: " << msg << endl;
  }
  if (strict) {
    cout << "Exiting..." << endl;
    exit(1);
  }
}

GetkwError::GetkwError(ostringstream & err) {
  this->msg = err.str();
  if (verbose or strict) {
    cout << "Error: " << msg << endl;
  }
  if (strict) {
    cout << "Exiting..." << endl;
    exit(1);
  }
}
GetkwError::~GetkwError() throw() {
  // TODO Auto-generated destructor stub
}

void GetkwError::trigger(const string & err) {
  msg = err;
  if (verbose or strict) {
    cout << "Error: " << msg << endl;
  }
  if (strict) {
    cout << "Exiting..." << endl;
    exit(1);
  }
  throw * this;
}

void GetkwError::setVerbose(bool flag) { verbose = flag; }

void GetkwError::setStrict(bool flag) { strict = flag; }
