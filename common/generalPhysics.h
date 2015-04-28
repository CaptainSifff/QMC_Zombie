/***************************************************************************
 *   Copyright (C) 2008 by Florian Goth   *
 *   fgoth@physik.uni-wuerzburg.de   *
 *                                                                         *
 *   All rights reserved.                                                  *
 *                                                                         *
 *   Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met: *
 *     * Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer. *
 *     * Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution. *
 *     * Neither the name of the <ORGANIZATION> nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission. *
 *                                                                         *
 *   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS   *
 *   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT     *
 *   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR *
 *   A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR *
 *   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, *
 *   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,   *
 *   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR    *
 *   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF *
 *   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING  *
 *   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS    *
 *   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.          *
 ***************************************************************************/
#ifndef GENERALPHYSICS_H
#define GENERALPHYSICS_H
#include <cmath>
#include "MTL/MTL_Macros.h"
/**
@file generalPhysics.h file to store some common physics things
*/

/**
an enum for spins
*/
enum SPINS
{
    UP = 1,
    DOWN =  -1
};

SPINS operator!(const SPINS& rhs)
{
  if (rhs == UP) return DOWN;
    else
      return UP;
}

/**
The fermi function
@return An occupation propability
*/
template <typename T>
inline T fermi(const T& e) throw()
{
    T xhalf = -e/T(2.0);
    T retval;
    if(unlikely(xhalf > std::log(0.01*std::numeric_limits<T>::max())))
      retval = 1.0;
    else
      if(unlikely(xhalf < std::log(10*std::numeric_limits<T>::min())))
	retval = 0.0;//FIXME: not yet decided whether setting this to zero is better than setting it to epsilon
	else
	  retval = 0.5*std::exp(xhalf)/std::cosh(xhalf);
    return retval;
}
#endif
