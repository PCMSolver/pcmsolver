/*
 * \file GetkwError.h
 *
 *  \date Jul 2, 2009
 *  \author Jonas Juselius <jonas.juselius@uit.no> \n
 *          CTCC, University of Tromsø
 *
 * \breif
 */

#ifndef GETKWERROR_H_
#define GETKWERROR_H_

#include <cstdlib>
#include <exception>
#include <iostream>
#include <sstream>
#include <string>

// TODO: define ABORT_EXCEPTION and redefine THROW not to throw, but abort.

#define THROW_GETKW(X)                                                              \
  {                                                                                 \
    std::ostringstream _err;                                                        \
    _err << "Error: " << __func__ << ",  line " << __LINE__ << " in  " << __FILE__  \
         << ": " << X << endl;                                                      \
    throw GetkwError(_err);                                                         \
  }

using namespace std;

class GetkwError : public exception {
public:
  GetkwError();
  GetkwError(const string & err);
  GetkwError(ostringstream & err);
  virtual ~GetkwError() throw();
  void trigger(const string & msg);
  static void setVerbose(bool flag);
  static void setStrict(bool flag);
  friend ostream & operator<<(ostream & o, const GetkwError & e) {
    return o << e.msg;
  }
  //	virtual const char *what() {
  //		return err.c_str();
  //	}
private:
  string msg;
  static bool verbose;
  static bool strict;
};

#endif /* GETKWERROR_H_ */
