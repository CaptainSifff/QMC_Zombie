/***************************************************************************
 *   Copyright (C) 2007 by Florian Goth   *
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
#ifndef GENERATORS_H
#define GENERATORS_H
namespace MTLICCL
{
template < typename T >
class None
{
private:
    T Null;
public:
    inline None(T N=0 ) : Null(N)
    {}
    inline T operator() (unsigned int i , unsigned int j )
    {
        return Null;
    }
};

template < typename T >
class Identity
{
private:
    T Eins;
    T Null;
public:
    inline Identity(T E = 1 , T N = 0) : Eins(E) , Null(N)
    {}
    inline T operator() (unsigned int i , unsigned int j )
    {
        return ( (i == j ) ?  Eins  : Null );
    }
};

/**
A Generator to initialize from a continous range of memory
*/
template< typename T >
class InitFromMem
{
private:
public:
    const T* ptr;
    unsigned int lines;
    inline InitFromMem(const T *const p , unsigned int l ) : ptr(p) , lines(l)
    {}
    inline T operator() (unsigned int i, unsigned int k )
    {
        return ptr[ k + i * lines ];
    }
};

template < typename T >
class Enumerator
{
public:
    inline Enumerator(unsigned int s): size(s)
    {}
    inline T operator() (unsigned int i, unsigned int k)
    {
        return i * size + k;
    }
private:
    unsigned int size;
};

template < typename T >
class TriGen
{
private:
    T L;
    T M;
    T R;
public:
    inline TriGen(T E , T N, T z ) : L(E) , M(N) , R(z)
    {}
    inline T operator() (unsigned int i , unsigned int j )
    {
        if ( i == j )
            return M;
        if ( j+1 == i )
            return R;
        if ( j-1 == i )
            return L;
        return 0;
    }
};

template <typename T>
class ScalingMatrixGenerator
{
public:
    inline ScalingMatrixGenerator(T* vals, unsigned int c) : diagonalelements(vals), cnt(c)
    {}
    inline T operator() (unsigned int i, unsigned int j)
    {
        if (i == j) return diagonalelements[i];
        return 0;
    }
private:
    T* diagonalelements;
    unsigned int cnt;
};
}
#endif
