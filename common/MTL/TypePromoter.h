/***************************************************************************
 *   Copyright (C) 2006 by Florian Goth   *
 *   CaptainSifff@gmx.de   *
 *                                                                         *
 *   Permission is hereby granted, free of charge, to any person obtaining *
 *   a copy of this software and associated documentation files (the       *
 *   "Software"), to deal in the Software without restriction, including   *
 *   without limitation the rights to use, copy, modify, merge, publish,   *
 *   distribute, sublicense, and/or sell copies of the Software, and to    *
 *   permit persons to whom the Software is furnished to do so, subject to *
 *   the following conditions:                                             *
 *                                                                         *
 *   The above copyright notice and this permission notice shall be        *
 *   included in all copies or substantial portions of the Software.       *
 *                                                                         *
 *   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,       *
 *   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF    *
 *   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. *
 *   IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR     *
 *   OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, *
 *   ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR *
 *   OTHER DEALINGS IN THE SOFTWARE.                                       *
 ***************************************************************************/
#include <limits>
#include "meta.h"
#ifndef TYPEPROMOTER_H
#define TYPEPROMOTER_H

template < class A , class B >
struct TypePromoter
{//Selects between the different Scalar Types 
   enum {
       max_exp_A = std::numeric_limits<A>::max_exponent10,
       max_exp_B = std::numeric_limits<B>::max_exponent10,
       digits_A  = std::numeric_limits<A>::digits,
       digits_B  = std::numeric_limits<B>::digits
   };
   typedef typename Select<
      // The type is considered to be smaller if the exponent is smaller
      // (integer types always have exponent == 0).
      max_exp_A < max_exp_B ||

      // If the exponents are equal the type with the smaller # of digits
      // is the smaller one.
      // This comparison will never happen between an integral and a non-integral
      // type since this case is already handled by the previous condition
      ( max_exp_A == max_exp_B && digits_A < digits_B ),
      B,
      A>::Result RET;
    //FIXME: Handle more Complicated Cases
};
#endif
