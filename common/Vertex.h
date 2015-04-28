/***************************************************************************
 *   Copyright (C) 2009 - 2013 by Florian Goth   *
 *   fgoth@wthp095   *
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
#ifndef VERTEX_H
#define VERTEX_H
#include <iostream>
#include <cmath>
#include <limits>
#include "MTL/precision.h"
#include "generalPhysics.h"

template <class Vertex, class Move, typename FPType>
class VertexDomain
{
};

template <typename FPType>
class Basic_Vertex
{
public:
    inline Basic_Vertex() throw() : tau(0), spin(UP)  {}
    inline Basic_Vertex(FPType t1, SPINS s) throw();
    inline Basic_Vertex(uint, FPType t1, SPINS s) throw();//HACK!!
    FPType tau;///< the imaginary time of the Vertex
    SPINS spin;///< the Ising-spin of the Vertex
    inline bool operator<(const Basic_Vertex& rhs) const throw();
    inline bool operator==(const Basic_Vertex& rhs) const throw();
    inline bool operator!=(const Basic_Vertex& rhs) const throw();
private:
};

template <class Move, typename FPType>
class VertexDomain<Basic_Vertex<FPType>, Move, FPType>
{
public:
    template<class PRNG>
//FIXME!!!!!!!!!! this function should not need the number of sites!!!!!!!!!!!!!!!!!!
    static inline Basic_Vertex<FPType> generateRandomVertex(FPType contourlen, int sites, PRNG& prng) throw()
    {
        return Basic_Vertex<FPType>(static_cast<FPType>(prng.rndfloat(contourlen)), (prng.template rndInteger<64>() < 32 ? DOWN : UP));
    }
};

template <typename FPType>
Basic_Vertex<FPType>::Basic_Vertex(FPType t1, SPINS s) throw() : tau(t1), spin(s) {}

template <typename FPType>
Basic_Vertex<FPType>::Basic_Vertex(uint, FPType t1, SPINS s) throw() :tau(t1), spin(s) {}

template <typename FPType>
bool Basic_Vertex<FPType>::operator==(const Basic_Vertex& rhs) const throw()
{   
    return (spin == rhs.spin) && fpequal(tau, rhs.tau);
}

template <typename FPType>
bool Basic_Vertex<FPType>::operator!=(const Basic_Vertex& rhs) const throw()
{
    return !(*this == rhs);
}

template <typename FPType>
bool Basic_Vertex<FPType>::operator<(const Basic_Vertex& rhs) const throw()
{
  //the following two volatile statements are necessary if the code is compiled with i387 FP math.
  //the volatile statements hopefully force the compiler to write the data back to the memory and hence truncate from the 80bit
  // FPU representation to the 64bit double representation. This issue is not triggered if the XMM registers(2* 64bit) are
  // used for floating point math.
  if(spin == rhs.spin)//try to defer the fp comparison as far back as possible...
  {
    volatile FPType tau1 = tau;
    volatile FPType tau2 = rhs.tau;
    return tau1 < tau2;
  }
  return spin < rhs.spin;
}

/**
The Basic Vertex of the Hubbard-Model
*/
template <typename FPType>
class Hubbard_Vertex : public Basic_Vertex<FPType>
{
public:
    /**
    Constructor of the Vertex.
    @param s1 the first lattice site    
    @param t1 the first time.
    @param s the Ising spin
    */
    inline Hubbard_Vertex(int s1, FPType t1, SPINS s) throw();
    /**
    An operator that compares two vertices lexicographically.
    @param rhs the other vertex to compare against.
    @return true if this vertex is smaller than the other else false
    */
    inline bool operator<(const Hubbard_Vertex& rhs) const throw();
    /**
    An operator that probes the identity of two Vertices
    @param rhs the Vertex to compare against
    @return true if the Vertices are equal in all of its members, else false.
    */
    inline bool operator==(const Hubbard_Vertex& rhs) const throw();
    inline bool operator!=(const Hubbard_Vertex& rhs) const throw();
    int site;///< the site of the Vertex
private:
};

template <typename FPType>
bool Hubbard_Vertex<FPType>::operator!=(const Hubbard_Vertex<FPType>& rhs) const throw()
{
    return  !(*this == rhs);
}

template <typename FPType>
bool Hubbard_Vertex<FPType>::operator==(const Hubbard_Vertex<FPType>& rhs) const throw()
{
    return  (site == rhs.site) && Basic_Vertex<FPType>::operator==(rhs);
}

template <typename FPType>
bool Hubbard_Vertex<FPType>::operator<(const Hubbard_Vertex<FPType>& rhs) const throw()
{
    if (site == rhs.site)
        return static_cast<Basic_Vertex<FPType>* >(*this)->operator<(rhs);
    return site < rhs.site;
}

/**
An overloaded stream operator for printing the Hubbard vertex to stdout
*/
template <typename FPType>
inline std::ostream& operator<<(std::ostream& out, const Hubbard_Vertex<FPType>& rhs )
{
    out<<"("<<rhs.site<<", "<<rhs.tau<<", "<<rhs.spin<<")";
    return out;
}

/**
An overloaded stream operator for printing the Basic vertex to stdout
*/
template <typename FPType>
inline std::ostream& operator<<(std::ostream& out, const Basic_Vertex<FPType>& rhs )
{
    out<<"("<<rhs.tau<<", "<<rhs.spin<<")";
    return out;
}

template <typename FPType>
Hubbard_Vertex<FPType>::Hubbard_Vertex(int si, FPType t, SPINS s) throw() : Basic_Vertex<FPType>(t, s), site(si)
{
}
#endif
