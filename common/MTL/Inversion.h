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
#include <cmath>
#include "Matrix_Header.h"
#include "Determinant.h"

namespace MTLICCL
{

template <typename T>
inline Matrix< Static<Config<T, 2, 2> > > inverse( const Matrix< Static<Config<T, 2, 2> > >& m) throw()
{
    typedef Matrix< Static< Config<T, 2, 2> > > ReturnType;
    ReturnType retval;
    T factor = 1.0/det(m);
    retval(1, 1, factor * m(0, 0));
    retval(0, 1, -factor * m(0, 1));
    retval(1, 0, -factor * m(1, 0));
    retval(0, 0, factor * m(1, 1));
    return retval;
}

template <template <class> class T, typename M>
inline Matrix< T<Dynamic<M> > > inverse(const Matrix<T<Dynamic<M> > >& m)
{
    typedef Matrix<T<Dynamic<M> > > RetType;
    const unsigned int size = m.Rows();
    RetType retVal(m);
    RetType lu(m);//get's overwritten anyway
    unsigned int idx[size];

    typename Matrix<T<Dynamic<M> > >::value_type col[size];
    ludecompose(lu, idx);
    for (unsigned int j = 0; j < size; j++)
    {
        for (unsigned int i = 0; i < size; i++) col[i] = 0.0;
        col[j] = 1.0;
        lubacksubstitute( lu, idx, col );
        for (unsigned int i = 0; i < size; i++)
            retVal(i,j, col[i]);
    }
    return retVal;
}

template <typename M>
inline Matrix< /*Dynamic<M>*/M > inversefromLU(const Matrix</*Dynamic<M>*/M >& lu, const unsigned int *const idx)
{
    typedef Matrix</*Dynamic<M>*/M > RetType;
    const unsigned int size = lu.Rows();
    RetType retVal(lu);
    typename Matrix</*Dynamic<M>*/ M>::value_type col[size];

    for (unsigned int j = 0; j < size; j++)
    {
        for (unsigned int i = 0; i < size; i++) col[i] = 0.0;
        col[j] = 1.0;

        lubacksubstitute(lu, idx, col );

        for (unsigned int i = 0; i < size; i++)
            retVal(i,j) = col[i];
    }

    return retVal;
}

template <typename M>
inline Matrix<M> inverse(const Matrix<M>& m)
{
    typedef Matrix<M> RetType;
    const unsigned int size = m.Rows();
    RetType retVal(m);
    RetType lu(m);
    unsigned int idx[size];
    typename Matrix<M>::value_type col[size];
    ludecompose(lu, idx);

    for (unsigned int j = 0; j < size; j++)
    {
        for (unsigned int i = 0; i < size; i++) col[i] = 0.0;
        col[j] = 1.0;

        lubacksubstitute(lu, idx, col );

        for (unsigned int i = 0; i < size; i++)
            retVal(i,j) = col[i];
    }

    return retVal;
}

inline Matrix< Static<Config<float, 4, 4> > >inverse( const Matrix< Static< Config<float, 4, 4> > >& m)
{
    //from GGems II
    typedef Matrix< Static< Config<float, 4, 4> > > ReturnType;
    MatrixCreator<ReturnType, InitFromMatrix<  ReturnType > > M;
    ReturnType RetVal = M.Create(4, 4, InitFromMatrix<ReturnType>(m, 0 ,0) );
    int pvt_i[4], pvt_j[4];            // Locations of pivot elements
    float pvt_val;               // Value of current pivot element
    float hold;                  // Temporary storage
    float determinat = 1.0f;
    for (unsigned int k = 0; k < 4; k++)
    {
        // Locate k'th pivot element
        pvt_val = RetVal(k,k);            // Initialize for search
        pvt_i[k]=k;
        pvt_j[k]=k;
        for (unsigned int i = k; i < 4; i++)
        {
            for ( unsigned int j = k; j < 4; j++)
            {
                if (fabs(RetVal(i,j)) > fabs(pvt_val))
                {
                    pvt_i[k]=i;
                    pvt_j[k]=j;
                    pvt_val = RetVal(i,j);
                }
            }
        }

        // Product of pivots, gives determinant when finished
        determinat *= pvt_val;
        if (fabs(determinat)< (1e-9) )
        {
            throw("Matrix not inverted!! It is singular!");  // Matrix is singular (zero determinant)
            //return;
        }

        // "Interchange" rows (with sign change stuff)
        unsigned int i = pvt_i[k];
        if (i != k)
        {
            // If rows are different
            for ( unsigned int j = 0; j < 4 ; j++)
            {
                hold = (-1)*RetVal(k,j);
                RetVal(k, j, RetVal(i,j) );
                RetVal(i, j, hold   );
            }
        }

        // "Interchange" columns
        unsigned int j = pvt_j[k];
        if (j != k)
        {
            // If columns are different
            for (i = 0; i < 4; i++)
            {
                hold =- RetVal(i,k);
                RetVal(i, k, RetVal(i,j));
                RetVal(i, j, hold);
            }
        }

        // Divide column by minus pivot value
        for (i = 0; i < 4; i++)
        {
            if (i != k)
                RetVal(i,k, RetVal(i,k) / ( -pvt_val) );
        }

        // Reduce the matrix
        for (i = 0; i < 4; i++)
        {
            hold = RetVal(i,k);
            for (j = 0; j < 4; j++)
            {
                if (i != k && j != k)
                    RetVal(i,j, RetVal(i,j) + hold * RetVal(k,j) );
            }
        }

        // Divide row by pivot
        for (j = 0; j < 4; j++)
        {
            if (j != k)
                RetVal(k, j, RetVal(k,j) / pvt_val );
        }

        // Replace pivot by reciprocal (at last we can touch it).
        RetVal(k,k, 1.0f / pvt_val );
    }

    // That was most of the work, one final pass of row/column interchange
    // to finish
    for ( int k = 4-2; k >= 0; k--)
    {
        // Don't need to work with 1 by 1 corner
        unsigned int i=pvt_j[k];            // Rows to swap correspond to pivot COLUMN
        if (static_cast<int>(i) != k)
        {
            // If rows are different
            for (unsigned int j = 0; j < 4; j++)
            {
                hold = RetVal(k,j);
                RetVal(k, j, -RetVal(i,j) );
                RetVal(i, j, hold);
            }
        }

        unsigned int j = pvt_i[k];           // Columns to swap correspond to pivot ROW
        if (static_cast<int>(j) != k)             // If columns are different
            for (i = 0; i < 4; i++)
            {
                hold=RetVal(i,k);
                RetVal(i,k, -RetVal(i,j));
                RetVal(i, j, hold);
            }
    }
    return RetVal;
}
}
