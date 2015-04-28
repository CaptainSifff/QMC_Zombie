/***************************************************************************
 *   Copyright (C) 2007,2009 by Florian Goth   *
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
#ifndef DETERMINANT_H
#define DETERMINANT_H
#include "Matrix_Header.h"
#include "Matrix.h"
#include <fstream>
#include <iostream>

namespace MTLICCL
{
template <class T>
inline T det( const Matrix< Static<Config<T, 2, 2> > > &MTL_RESTRICT m) throw()
{
    return m(0,0) * m(1,1) - m(0,1) * m(1,0);
}

template <class T>
inline T det( const Matrix< Static<Config<T, 3, 3> > > &MTL_RESTRICT m) throw()
{
    return m(0,0) * m(1,1) * m(2,2) + m(0,1) * m(1,2) * m(2,0) + m(0,2) * m(1,0) * m(2,1) - m(0,0) * m(1,2) * m(2,1) - m(0,1) * m(1,0) * m(2,2) - m(0,2) * m(1,1) * m(2,0);
}

/**
This function determines the determinant of a Matrix, if you add one column: u2 and one row : v1 and the new Element in the lower right corner mnn
@param m an LR decomposed Matrix
@param idx the permutation array
@param detm the previous determinant of m
@param v1 the new row
@param u2 the new column
@param mnn the new Element in the lower right corner
*/
template <class A>
inline typename Matrix<A>::Elementtype det( const Matrix<A> &MTL_RESTRICT m, const unsigned int *const MTL_RESTRICT idx, const typename Matrix<A>::Elementtype &MTL_RESTRICT detm, const typename Matrix<A>::Elementtype *const v1, const typename Matrix<A>::Elementtype *const u2, const typename Matrix<A>::Elementtype mnn)
{
  typename Matrix<A>::Elementtype RetVal(mnn);
  if(likely(m.Rows() != 0))
  {
  typename Matrix<A>::Elementtype utemp[m.Rows()];
  for(unsigned int k = 0; k < m.Rows(); ++k)//copy Elements because lubacksubstitute overwrites u2 with the solution vector
    utemp[k] = u2[k];
  lubacksubstitute(m, idx, utemp);
  for(unsigned int k = 0; k < m.Rows(); ++k) RetVal -= utemp[k]*v1[k];
    RetVal *= detm;
  }
  else
  {//Nothing, it's only a 1x1 Matrix in this case
  }
  return RetVal;
}

template<class A>
inline typename Matrix<A>::Elementtype det( Matrix<A> m)//whole copy is necessary, because ludecompose destroys it
{
    typedef typename Matrix<A>::Elementtype ScalarType;
    ScalarType RetVal;
    if ( likely(m.Rows() > 1))
    {
        unsigned int *MTL_RESTRICT idx = new unsigned int[m.Rows()];
        RetVal = ludecompose(m, idx);
        delete [] idx;
        for (unsigned int j = 0; j < m.Rows(); ++j)
            RetVal *= m(j,j);
    }
    else
    {
        if (m.Rows() == 1)
            RetVal = m(0,0);
        else
            RetVal = 0;//the empty Matrix Case... has an empty Matrix a determinant of 0 ???
    }
    return RetVal;
}
/*
template <class T>
inline T det( const Matrix< Static< Config<T, 4, 4> > >& m);//Not yet Implemented
*/
}
#endif
