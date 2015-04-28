/***************************************************************************
 *   Copyright (C) 2007-2013 by Florian Goth   *
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
#ifndef MTL_COMMON_H
#define MTL_COMMON_H
#include "MTL_Macros.h"
#if (GCC_VERSION > GCC_VER(4,5,0))
#include <x86intrin.h>
#elif defined(__SSE2__)
#include <emmintrin.h>
#endif
#include "meta.h"

namespace MTLICCL
{
template < class A , class B >
struct MatchingSizes
{
    typedef Int2Type<A::Config::columns> o;
    typedef Int2Type<B::Config::rows> p;
    int z[ ( o::value == p::value  ) ? 1 : -1 ];
    //HACK!! THIS GIVES A COMPILE TIME ERROR IF THE SIZES DON'T MATCH!
};

template< class M >
class SubMatrix
{
private:
    const M& A;//A Reference to the Matrix, that will be a SubMatrix
    unsigned int i_0, k_0;//The Origins(upper left) of the Matrix.
    unsigned int deltai, deltak;//how far should this Matrix extend;
public:
    SubMatrix(const M& a, unsigned int b, unsigned int c, unsigned int di, unsigned int dk) : A(a), i_0(b), k_0(c), deltai(di), deltak(dk)
    {}
    bool HasIndex(unsigned int i, unsigned int k) const
    {
        return ( i >= i_0 ) && ( k >= k_0 ) && ( i < (i_0 + deltai) ) && ( k < (k_0 + deltak) );
    }
    Scalar(M) operator() (unsigned int i, unsigned int k) const
    {
        return A(i-i_0,k-k_0);
    }
    inline unsigned int Length()
    {
        return A.Rows();
    }//Annahme: A ist quadratisch
};

template< class M >
class InitFromMatrix
{
private:
    typedef Scalar(M) ScalarType;
    const M& A;
    unsigned int Origini,Origink;
public:
    inline InitFromMatrix(const M& arg ,unsigned int a, unsigned int b) : A(arg), Origini(a), Origink(b)
    {}
    inline ScalarType operator() (unsigned int i, unsigned int k )
    {
        return A(Origini + i, Origink + k);
    }
};

typedef double v2sd __attribute__ ( ( vector_size ( 16 ) ) );
class Vec2
{// a class for a vector of two thingy that is vectorization enabled
public:
    inline const Vec2& operator= ( const Vec2& rhs )
    {
        v = rhs.v;
        return *this;
    }
    inline void loadpd ( const double *const p )
    {
        v = _mm_load_pd ( p );
    }
    inline void loadupd ( const double *const p )
    {
        v = _mm_loadu_pd ( p );
    }
    inline void storepd ( double *const p ) const
    {
        _mm_store_pd ( p, v );
    }
    inline Vec2 operator+= ( const Vec2& rhs )
    {
        v += rhs.v;
        return *this;
    }
    inline void zero()
    {
        v = _mm_setzero_pd();
    }
    union
    {
        double x[2];
        v2sd v;
    };
};

}
#endif
