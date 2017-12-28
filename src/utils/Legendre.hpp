//==================================================================
/*
 *  Legendre.hpp -- C++ functions to evaluate Legendre polynomials
 *
 *  Copyright (C) 2014 by James A. Chappell
 *
 *  Permission is hereby granted, free of charge, to any person
 *  obtaining a copy of this software and associated documentation
 *  files (the "Software"), to deal in the Software without
 *  restriction, including without limitation the rights to use,
 *  copy, modify, merge, publish, distribute, sublicense, and/or
 *  sell copies of the Software, and to permit persons to whom the
 *  Software is furnished to do so, subject to the following
 *  condition:
 *
 *  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 *  EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 *  OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 *  NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 *  HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 *  WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 *  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 *  OTHER DEALINGS IN THE SOFTWARE.
 */
//=================================================================
/*
 * legendre.h:  Version 0.02
 * Created by James A. Chappell <rlrrlrll@gmail.com>
 * http://www.storage-b.com/math-numerical-analysis/18
 * Created 29 September 2005
 *
 * History:
 * 29-sep-2005  created
 * 14-nov-2014  templates
 */
//==============

#pragma once

/*
 *  Function calculates Legendre Polynomials Pn(x)
 */

namespace Legendre {
// n = 0
template <typename T> inline T P0(const T & /* x */) { return static_cast<T>(1.0); }

// n = 1
template <typename T> inline T P1(const T & x) { return x; }

// n = 2
template <typename T> inline T P2(const T & x) {
  return ((static_cast<T>(3.0) * x * x) - static_cast<T>(1.0)) * static_cast<T>(0.5);
}

/*
 *  Pn(x)
 */
template <typename T> inline T Pn(unsigned int n, const T & x) {
  if (n == 0) {
    return P0<T>(x);
  } else if (n == 1) {
    return P1<T>(x);
  } else if (n == 2) {
    return P2<T>(x);
  }

  T pnm1(P2<T>(x));
  T pnm2(P1<T>(x));
  T pn(pnm1);

  for (unsigned int l = 3; l <= n; l++) {
    pn = (((static_cast<T>(2.0) * static_cast<T>(l)) - static_cast<T>(1.0)) * x *
              pnm1 -
          ((static_cast<T>(l) - static_cast<T>(1.0)) * pnm2)) /
         static_cast<T>(l);
    pnm2 = pnm1;
    pnm1 = pn;
  }

  return pn;
}
} // namespace Legendre
