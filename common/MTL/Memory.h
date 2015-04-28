/***************************************************************************
 *   Copyright (C) 2006,2009 by Florian Goth   *
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
#ifndef MEMORY_H
#define MEMORY_H
#include "meta.h"
#include "TypePromoter.h"
#include "MTL_Common.h"
#include "MemArray1D.h"

namespace MTLICCL
{

template <
typename T_Elementtype , int _rows , int _columns >
struct Config
{
    typedef T_Elementtype ElementType;
    enum {
        Is_Static = true,
        Is_Dynamic = false,
        Is_OneD = false,
        rows = _rows,
        columns = _columns
    };
};

template< class A, class B, template<typename, typename> class Sizer >//A and B are Matrices
struct OnStatic
{
    typedef Sizer<A,B> S;
    typedef Config< typename TypePromoter< Scalar(A) , Scalar(B)>::RET ,S::Rows ,S::Columns > RET;
};

template < class Generator >
class Dynamic
{
public:
    typedef Generator Config;
    typedef typename Config::ElementType Elementtype;

private:
    unsigned int _rows;
    unsigned int _columns;
    Elementtype** arr;

public:
    inline unsigned int  Rows() MTL_RESTRICT const throw()
    {
        return _rows;
    }
    inline unsigned int Columns() const MTL_RESTRICT throw()
    {
        return _columns;
    }
    inline Dynamic(unsigned int r , unsigned int c) : _rows(r), _columns(c)
    {
        arr = new Elementtype* [Rows()];
        for ( unsigned int k = 0 ; k < Rows() ; ++ k)
            arr[k] = new Elementtype[Columns()];
    }
    inline Dynamic( const Dynamic&MTL_RESTRICT arg ) : _rows(arg.Rows()), _columns(arg.Columns())
    {
        arr = new Elementtype* [Rows()];
        for ( unsigned int k = 0 ; k < Rows() ; ++ k)
            arr[k] = new Elementtype[Columns()];
        for ( unsigned int k = 0 ; k < Rows() ; ++k )
            for ( unsigned int i = 0 ; i < Columns() ; ++i )
                SetElement( k , i , arg(k,i) );
        return;
    }
    inline Dynamic& operator=(const Dynamic& rhs)
    {
        if (this != &rhs)
        {
            if ((this->Rows() != rhs.Rows()) || (this->Columns() != rhs.Columns()))
            {//If we copy around matrices of the same size we don't need that step
                for ( unsigned int k = 0 ; k < Rows() ; ++ k)
                    delete [] arr[k];
                delete [] arr;
                _rows = rhs.Rows();
                _columns = rhs.Columns();
                arr = new Elementtype* [Rows()];
                for ( unsigned int k = 0 ; k < Rows() ; ++ k)
                    arr[k] = new Elementtype[Columns()];
            }
            for ( unsigned int k = 0 ; k < Rows() ; ++k )
                for ( unsigned int i = 0 ; i < Columns() ; ++i )
                    SetElement( k , i , rhs(k,i) );
        }
        return *this;
    }
#ifdef HAS_RVALUE_REFERENCES
    inline Dynamic( Dynamic&& arg ) : _rows(arg.Rows()), _columns(arg.Columns())
    {
        arr = new Elementtype* [Rows()];
        for ( unsigned int k = 0 ; k < Rows() ; ++ k)
        {
            arr[k] = arg.arr[k];
            arg.arr[k] = NULL;
        }
        return;
    }
    inline Dynamic& operator=(Dynamic&& rhs)
    {
        if (this != &rhs)
        {
            if ((this->Rows() != rhs.Rows()) || (this->Columns() != rhs.Columns()))
            {//If we copy around matrices of the same size we don't need that step
                for ( unsigned int k = 0; k < Rows(); ++ k)
                    delete [] arr[k];
                delete [] arr;
                _rows = rhs.Rows();
                _columns = rhs.Columns();
                arr = new Elementtype* [Rows()];
            }
            for ( unsigned int k = 0 ; k < Rows() ; ++k )
            {
                arr[k] = rhs.arr[k];
                rhs.arr[k] = NULL;
            }
        }
        return *this;
    }
#endif
    inline ~Dynamic()
    {
        for ( unsigned int k = 0 ; k < Rows() ; ++ k)
            delete [] arr[k];
        delete [] arr;
    }
    inline const Elementtype& operator() (unsigned int i , unsigned int k) const throw()
    {//const Version
        return arr[i][k];
    }
    inline Elementtype& operator() (unsigned int i , unsigned int k) throw()
    {//Non-const Version
        return arr[i][k];
    }
    inline const Elementtype& GetElement (unsigned int i , unsigned int k) const throw()
    {
        return arr[i][k];
    }
    inline void SetElement (unsigned int i , unsigned int k, Elementtype val) throw()
    {
        arr[i][k] = val;
        return;
    }
};

template < class Generator >
class Static
{
public:
    typedef Generator Config;
    typedef typename Config::ElementType Elementtype;

    inline unsigned int Rows() const throw()
    {
        return _rows;
    }
    inline unsigned int Columns() const throw()
    {
        return _columns;
    }
    inline Elementtype& operator() (unsigned int i , unsigned int k) throw()
    {//Non-const Version
        return arr[i*_columns + k];
    }
    inline const Elementtype& operator() (unsigned int i , unsigned int k) const throw()
    {//Const Version
        return arr[i*_columns + k];
    }
    inline const Elementtype& GetElement (unsigned int i , unsigned int k) const throw()
    {
        return arr[i*_columns + k];
    }
    inline void SetElement (unsigned int i , unsigned int k, Elementtype val ) throw()
    {
        arr[i*_columns + k] = val;
        return;
    }
    inline Static() throw()
    {}
    inline Static( const Static& arg ) throw()
    {
        for ( unsigned int k = 0 ; k < Rows() ; ++k )
            for ( unsigned int i = 0 ; i < Columns() ; ++i )
                SetElement( k , i , arg(k,i) );
        return;
    }
    inline Static& operator=(const Static& arg) throw()
    {
    if(this != &arg)
    {
        for ( unsigned int k = 0 ; k < Rows() ; ++k )
            for ( unsigned int i = 0 ; i < Columns() ; ++i )
                arr[k*_columns + i] = arg(k,i);
    }
    return *this;
    }
private:
    enum { _rows = Config::rows,
           _columns = Config::columns
         };
    Elementtype arr[_rows*_columns];
};

template <
typename T_Elementtype >
struct Config_Dynamic
{
    typedef T_Elementtype ElementType;
    enum {
        Is_Static = false,
        Is_Dynamic = true,
        Is_OneD = false,
        rows = 0,
        columns = 0,
    };
};

template< class A, class B>
struct ChooseDynamic
{
    typedef typename Select< (A::Config::Is_OneD == false) && (B::Config::Is_OneD == false),// if A and b are 2D Dynamic
    Dynamic< Config_Dynamic<typename TypePromoter< Scalar(A) , Scalar(B) >::RET> >,//return a 2D dynamic matrix
    MemArray1D< Config1D<typename TypePromoter< Scalar(A) , Scalar(B) >::RET> >
    >::Result RET;
};

template< class A , class B, template<typename, typename> class Sizer >
struct ConfigParser
{//This decides which Configurator to use and builts MatrixTypes According to this.
    typedef typename Select<( A::Config::Is_Static == true ) && ( B::Config::Is_Static == true ),/*Are both Matrices static ones?*/
    Static< typename OnStatic<A , B, Sizer >::RET >,/*If so then return a new static Matrix Type*/
    typename ChooseDynamic<A, B>::RET
    //Dynamic< Config_Dynamic<typename TypePromoter< Scalar(A) , Scalar(B) >::RET> >/*Else return a Dynamic One*/
    >::Result RET;//This contains now a fitting Config struct
    //FIXME: BoundsChecker??
};

template<class A, class Initializer>
struct MatrixCreator
{//Select the appropriate Constructor
private:
    inline A p_Create(unsigned int r, unsigned int c, Initializer G, Int2Type<true>)
    {//static
//        std::cout<<"returning static Matrix"<<std::endl;
        return A(G);
    }
    inline A p_Create(unsigned int r, unsigned int c, Initializer G, Int2Type<false>)
    {//dynamic
//        std::cout<<"returning dynamic Matrix"<<std::endl;
        return A(r,c,G);
    }
public:
    inline A Create(unsigned int r, unsigned int c, Initializer G)
    {
        return p_Create( r, c, G, Int2Type<A::Config::Is_Static>() );
    }
};

template < class A , class B >
struct IdentitySizer
{
    enum {
        Rows = A::Config::rows,
        Columns = A::Config::columns
    };
};

};
#endif
