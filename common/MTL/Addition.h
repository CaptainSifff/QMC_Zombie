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
#ifndef ADDITION_H
#define ADDITION_H
#include "Matrix_Header.h"

namespace MTLICCL
{

template < class A, class B >
struct AdditionSizer
{
    enum {
        Rows = B::Config::rows,
        Columns = A::Config::columns
    };
    /*char z[(Columns == B::Config::columns)? 1 : 0];
    char u[(A::Config::rows == Rows)? 1 : 0];*/
    //HACK!!
};

template < class A, class B >
inline Matrix< typename ConfigParser<Matrix<A>, Matrix<B>, AdditionSizer >::RET > operator+ ( const Matrix<A>& lhs , const Matrix<B>& rhs)
{//Matrix Addition
    typedef typename TypePromoter<Scalar(Matrix<A>), Scalar(Matrix<B>)>::RET ScalarType;
    //typedef typename TypeTraits<A>::NonConstType RetValType;
    typedef Matrix< typename ConfigParser<Matrix<A>, Matrix<B>, AdditionSizer >::RET > ReturnType;
    const unsigned int NrOfRows = lhs.Rows();
    const unsigned int NrOfColumns = rhs.Columns();
    MatrixCreator<ReturnType, None<ScalarType> > M;
    ReturnType RetVal = M.Create( NrOfRows , NrOfColumns , None<ScalarType>(0) );
    if (unlikely( (lhs.Rows() != rhs.Rows()) || (lhs.Columns() != rhs.Columns() ) ))
        throw(" Matrices incompatibel!!");
    for ( unsigned int r = 0 ; r < NrOfRows ; ++r )
        for ( unsigned int c = 0 ; c < NrOfColumns ; ++c)
            RetVal( r , c , lhs.GetElement( r , c ) + rhs.GetElement( r , c ) );
    return RetVal;
}

template < class A , class B >
inline Matrix< typename ConfigParser<Matrix<A>, Matrix<B>, AdditionSizer >::RET > operator- ( const Matrix<A>& lhs , const Matrix<B>& rhs)
{//Matrix-Subtraction
    typedef typename TypePromoter<Scalar(Matrix<A>), Scalar(Matrix<B>)>::RET ScalarType;
    typedef Matrix< typename ConfigParser<Matrix<A>, Matrix<B>, AdditionSizer >::RET > ReturnType;
    const unsigned int NrOfRows = rhs.Rows();
    const unsigned int NrOfColumns = rhs.Columns();
    MatrixCreator<ReturnType, None<ScalarType> > M;
    ReturnType RetVal = M.Create( NrOfRows , NrOfColumns , None<ScalarType>(0) );
    if (unlikely( (lhs.Rows() != rhs.Rows()) || (lhs.Columns() != rhs.Columns() ) ))
        throw(" Matrices incompatibel!!");
    for ( unsigned int r = 0 ; r < lhs.Rows() ; ++r )
        for ( unsigned int c = 0 ; c < lhs.Columns() ; ++c)
            RetVal( r , c , lhs.GetElement( r , c ) - rhs.GetElement( r , c ) );
    return RetVal;
}
}
#endif
