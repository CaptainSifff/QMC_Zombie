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
#ifndef MULTIPLY_H
#define MULTIPLY_H
#include "Matrix_Header.h"
#include "MTL_Common.h"
#include "Matrix.h"
namespace MTLICCL
{

template < class A, class B >
struct MultiplicationSizer
{
    typedef MatchingSizes< A , B > SizeDummy;
    enum {
        Rows = A::Config::rows,
        Columns = B::Config::columns,
        o = A::Config::columns,
        p = B::Config::rows,
        SizesMatch = (o==p)
    };
    //int z[ o == p  ? 1 : -1 ];
};

template<class M>
class Default_Patcher
{
public:
    inline Scalar(M) Get (unsigned int i, unsigned int k, SubMatrix<M>* S, unsigned int cnt) const
    {
        for (unsigned int c = 0; c < cnt; ++c)
            if (S[c].HasIndex(i,k) == true)
                return S[c](i,k);
    }
};

template<class M>
class Strassen_Patcher
{
public:
    inline Scalar(M) Get (unsigned int i, unsigned int k, SubMatrix<M>* S, unsigned int cnt) const
    {
        return S[(i / S[0].Length()) + 2 * (k / S[0].Length())](i,k);
    }
};

template< class M, template<class> class Patching_Trait = Default_Patcher >
class PatchWorkfromMatrices : public Patching_Trait<M>
{
private:
    typedef Scalar(M) ScalarType;
    SubMatrix<M>* Sub;
    unsigned int cnt;
public:
    PatchWorkfromMatrices(SubMatrix<M> S[],unsigned int c): Sub(S), cnt(c)
    {}
    ScalarType operator() (unsigned int i, unsigned int k) const
    {
        return Get(i,k,Sub,cnt);
    }
};

template < class A, class B >
struct StrassenSizer
{
    typedef MatchingSizes< A , B > SizeDummy;
    enum {
        Rows = A::Config::rows/2,
        Columns = B::Config::columns/2,
        o = A::Config::columns,
        p = B::Config::rows,
        SizesMatch = (o==p)
    };
    //int z[ o == p  ? 1 : -1 ];
};

template < class A, class B >
inline Matrix< typename ConfigParser<A,B, MultiplicationSizer >::RET > StrassensMultiply( const A& lhs, const B& rhs)
{/*Matrix Multplication according to Strassen
                                                Complexity: O( n^2.8 )
                                                See: H�mmerlin/Hoffmann: Numerische Mathematik
                                                */
    typedef typename TypePromoter<Scalar(A), Scalar(B)>::RET ScalarType; //The used type of the scalars
    typedef Matrix< typename ConfigParser<A,B, MultiplicationSizer >::RET > ReturnType;
    typedef Matrix< typename ConfigParser<A,B, StrassenSizer >::RET > HilfsType;
    MatrixCreator<HilfsType, InitFromMatrix<A> > HA;
    MatrixCreator<HilfsType, InitFromMatrix<B> > HB;
    //Die Matrizen A und B zerlegen
    HilfsType A11 = HA.Create(lhs.Rows() / 2 , lhs.Columns() / 2 , InitFromMatrix<A>(lhs, 0 ,0 ) );
    HilfsType A12 = HA.Create(lhs.Rows() / 2 , lhs.Columns() / 2 , InitFromMatrix<A>(lhs, 0 ,lhs.Columns()/2 ) );
    HilfsType A21 = HA.Create(lhs.Rows() / 2 , lhs.Columns() / 2 , InitFromMatrix<A>(lhs, lhs.Rows()/2 ,0 ) );
    HilfsType A22 = HA.Create(lhs.Rows() / 2 , lhs.Columns() / 2 , InitFromMatrix<A>(lhs, lhs.Rows()/2 ,lhs.Columns() / 2 ) );

    HilfsType B11 = HB.Create(rhs.Rows() / 2 , rhs.Columns() / 2 , InitFromMatrix<B>(rhs, 0 ,0 ) );
    HilfsType B12 = HB.Create(rhs.Rows() / 2 , rhs.Columns() / 2 , InitFromMatrix<B>(rhs, 0 ,rhs.Columns()/2 ) );
    HilfsType B21 = HB.Create(rhs.Rows() / 2 , rhs.Columns() / 2 , InitFromMatrix<B>(rhs, rhs.Rows() / 2, 0 ) );
    HilfsType B22 = HB.Create(rhs.Rows() / 2 , rhs.Columns() / 2 , InitFromMatrix<B>(rhs, rhs.Rows() / 2 ,rhs.Columns()/2 ) );
    //Die HilfsMatrizen berechnen
    HilfsType M1 = (A12 - A22) * (B21 + B22);
    HilfsType M2 = (A11 + A22) * (B11 + B22);
    HilfsType M3 = (A11 - A21) * (B11 + B12);
    HilfsType M4 = (A11 + A12) * B22;
    HilfsType M5 = A11 * (B12 - B22);
    HilfsType M6 = A22 * (B21 - B11);
    HilfsType M7 = (A21 + A22) * B11;
    //Zusammensetzen
    //Erzeugen der ErgebnisTeil-Matrizen
    HilfsType D1 = M1+M2-M4+M6;
    HilfsType D2 = M4+M5;
    HilfsType D3 = M6+M7;
    HilfsType D4 = M2-M3+M5-M7;

    SubMatrix<HilfsType> C[] = {
                                   SubMatrix<HilfsType>(D1,         0,      0     ,M1.Rows(),M1.Columns()),//C11 //M1 should have the correct Size;
                                   SubMatrix<HilfsType>(D2,         0,M1.Columns(),M4.Rows(),M4.Columns()),//C12 //M4 should have the correct Size
                                   SubMatrix<HilfsType>(D3, M1.Rows(),   0        ,M6.Rows(),M6.Columns()),//C21 //M6 should have the correct Size
                                   SubMatrix<HilfsType>(D4, M1.Rows(),M1.Columns(),M2.Rows(),M2.Columns()),//C22 //M2 should have the correct Size
                               };
    MatrixCreator<ReturnType, PatchWorkfromMatrices<HilfsType> > MC;
    return MC.Create(lhs.Rows(), rhs.Columns(), PatchWorkfromMatrices<HilfsType/*, Strassen_Patcher*/ >(C, 4));
}

template < class A >
inline A Multiply ( typename A::Elementtype scalar, const A& rhs, Int2Type<true> )
{//Scalar Multiply
    typedef typename A::Elementtype ScalarType;
    typedef Matrix< typename ConfigParser<A, A, IdentitySizer >::RET > ReturnType;
    unsigned int NrOfRows = rhs.Rows();
    unsigned int NrOfColumns = rhs.Columns();
    MatrixCreator<ReturnType, None<ScalarType> > M;
    ReturnType RetVal = M.Create( NrOfRows , NrOfColumns , None<ScalarType>(0) );
    for ( unsigned int r = 0 ; r < rhs.Rows() ; ++r )
        for ( unsigned int c = 0 ; c < rhs.Columns() ; ++c)
            RetVal( r , c , scalar * rhs.GetElement(r,c) );
    return RetVal;
}

template< class A, class B >
struct MultiplyReturnValueWrapper
{//if(A ist ScalarTyp von B) then Typ(skalare Multiplikation) = B else Matrixmultiplikation
    enum {Sametype = Conversion<A, typename B::Elementtype >::sametype,};
    typedef typename Select<Sametype,B ,A >::Result lht;//A further Level of indirection was needed to avoid a compiler-error
    typedef B ScalarReturnType;
    typedef Matrix< typename ConfigParser< lht, B, MultiplicationSizer >::RET > MatrixReturntype;
    /*Explanation for the above Line:
    In a scalar Multiplication A would be e.g. a double, thus ConfigParser would not compile. We now do the following:
    A is the Scalartype of B => we call Configparser<B,B> (not very useful, but we don't use it's Returnvalue) and  return B
    A is NOT the Scalartype of B => We suspect A is a Matrix and use normal Matrix-Matrix-Multiplication 
    */
    typedef typename Select< Sametype, ScalarReturnType, MatrixReturntype>::Result RET;
};

template <class RetType, class srcA, class srcB>
class Multiplication_specialization
{
public:
  static inline void mult(RetType& retval, const srcA& lhs, const srcB& rhs)
  {//generic multiplication, a bit optimized for superscalar CPUs
    typedef typename TypePromoter<Scalar(srcA), Scalar(srcB)>::RET ScalarType;
    const uint maxk2 = rhs.Columns()/4;
    const uint maxk = maxk2*4;
    const uint maxj2 = rhs.Rows()/4;
    const uint maxj = maxj2*4;
    for ( unsigned int i = 0 ; i < lhs.Rows(); ++i)
    {
        for ( unsigned int k2 = 0 ; k2 < maxk2; /*k+=4*/ ++k2 )//let's assume we can schedule four additions of anything
        {
	    uint k = 4*k2;//removes gcc warning about not being able to optimize the loop
            ScalarType temp0[4] = {0};
	    ScalarType temp1[4] = {0};
	    ScalarType temp2[4] = {0};
	    ScalarType temp3[4] = {0};
            for (unsigned int j2 = 0; j2 < maxj2; ++j2/*j+=4*/)
	    {
	      uint j = 4*j2;//removes gcc warning about not being able to optimize the loop
	      ScalarType lhs_ij0 = lhs.GetElement(i,j + 0);
	      ScalarType lhs_ij1 = lhs.GetElement(i,j + 1);
	      ScalarType lhs_ij2 = lhs.GetElement(i,j + 2);
	      ScalarType lhs_ij3 = lhs.GetElement(i,j + 3);
	      
              temp0[0] += lhs_ij0 * rhs.GetElement(j + 0,k + 0);
	      temp0[1] += lhs_ij0 * rhs.GetElement(j + 0,k + 1);
	      temp0[2] += lhs_ij0 * rhs.GetElement(j + 0,k + 2);
	      temp0[3] += lhs_ij0 * rhs.GetElement(j + 0,k + 3);
	      
	      temp1[0] += lhs_ij1 * rhs.GetElement(j + 1,k + 0);
	      temp1[1] += lhs_ij1 * rhs.GetElement(j + 1,k + 1);
	      temp1[2] += lhs_ij1 * rhs.GetElement(j + 1,k + 2);
	      temp1[3] += lhs_ij1 * rhs.GetElement(j + 1,k + 3);
	      
	      temp2[0] += lhs_ij2 * rhs.GetElement(j + 2,k + 0);
	      temp2[1] += lhs_ij2 * rhs.GetElement(j + 2,k + 1);
	      temp2[2] += lhs_ij2 * rhs.GetElement(j + 2,k + 2);
	      temp2[3] += lhs_ij2 * rhs.GetElement(j + 2,k + 3);
	      
	      temp3[0] += lhs_ij3 * rhs.GetElement(j + 3,k + 0);
	      temp3[1] += lhs_ij3 * rhs.GetElement(j + 3,k + 1);
	      temp3[2] += lhs_ij3 * rhs.GetElement(j + 3,k + 2);
	      temp3[3] += lhs_ij3 * rhs.GetElement(j + 3,k + 3);
	    }
	    for(uint t = 0; t < 4; ++t)
	    {
	      temp0[t] += temp1[t] + temp2[t] + temp3[t];
	    }
	    for (unsigned int j = maxj; j < rhs.Rows() ; ++j)
	    {
	      ScalarType lhs_ij = lhs.GetElement(i,j);
              temp0[0] += lhs_ij * rhs.GetElement(j,k + 0);
	      temp0[1] += lhs_ij * rhs.GetElement(j,k + 1);
	      temp0[2] += lhs_ij * rhs.GetElement(j,k + 2);
	      temp0[3] += lhs_ij * rhs.GetElement(j,k + 3);
	    }
            retval(i,k + 0,temp0[0]);
	    retval(i,k + 1,temp0[1]);
	    retval(i,k + 2,temp0[2]);
	    retval(i,k + 3,temp0[3]);
        }
        for ( unsigned int k = maxk; k < rhs.Columns(); ++k )
        {
            ScalarType temp[4] = {0};
            for (unsigned int j = 0; j < maxj; j+=4)
	    {
                temp[0] += lhs.GetElement(i,j + 0) * rhs.GetElement(j + 0,k);
		temp[1] += lhs.GetElement(i,j + 1) * rhs.GetElement(j + 1,k);
		temp[2] += lhs.GetElement(i,j + 2) * rhs.GetElement(j + 2,k);
		temp[3] += lhs.GetElement(i,j + 3) * rhs.GetElement(j + 3,k);
	    }
	    ScalarType temp2 = (temp[0] + temp[1]) + (temp[2] + temp[3]);
	    for (unsigned int j = maxj ; j < rhs.Rows() ; ++j)
                temp2 += lhs.GetElement(i,j) * rhs.GetElement(j,k);
            retval(i,k, temp2);
        }
    }
  }
};

#if 0
//__SSE2__
template <>
class Multiplication_specialization<MTLICCL::Matrix<MTLICCL::MemArray1D<MTLICCL::Config1D<double> > >, MTLICCL::Matrix<MTLICCL::MemArray1D<MTLICCL::Config1D<double> > >, MTLICCL::Matrix<MTLICCL::MemArray1D<MTLICCL::Config1D<double> > > >
{
public:
  static inline void mult(MTLICCL::Matrix<MTLICCL::MemArray1D<MTLICCL::Config1D<double> > >& retval, const MTLICCL::Matrix<MTLICCL::MemArray1D<MTLICCL::Config1D<double> > >& lhs, const MTLICCL::Matrix<MTLICCL::MemArray1D<MTLICCL::Config1D<double> > >& rhs)
  {
    constexpr size_t CLS = 64;//Cachelinesize is assumed to be 64bytes. Valid for almost all CPUs(known exceptions: AMD K7, Intel Atom, AMD E350)
    constexpr size_t DPC = CLS/(sizeof(double));//DPC = Double per cacheline
    const uint maxk = (rhs.Columns()/(DPC/2))*(DPC/2);
    const uint maxj = (rhs.Rows()/DPC)*DPC;
    for ( unsigned int i = 0 ; i < lhs.Rows(); ++i)
    {
        for ( unsigned int k = 0; k < maxk; k+=(DPC/2))
        {
            Vec2 temp0[DPC/2] = {0};
	    Vec2 temp1[DPC/2] = {0};
	    Vec2 temp2[DPC/2] = {0};
	    Vec2 temp3[DPC/2] = {0};
	    const double* lhs_ptr = lhs.arr + i * lhs.Columns();
            for (unsigned int j = 0; j < maxj; j+=DPC, lhs_ptr += DPC)
	    {
	      Vec2 lhs_ij[4];
	      for(uint t = 0; t < DPC/2; ++t)
	        lhs_ij[t].loadpd(/*lhs.arr + i * lhs.Columns() +*/ lhs_ptr + 2*t);
	      Vec2 b[4];
	      double* rhs_ptr_a = rhs.arr + (j + 2*0)*rhs.Columns() + k;
	      double* rhs_ptr_b = rhs.arr + (j + 2*0 + 1)*rhs.Columns() + k;
	      b[0].v = _mm_setr_pd(*(rhs_ptr_a) , *(rhs_ptr_b));
	      b[1].v = _mm_setr_pd(*(rhs_ptr_a + 1), *(rhs_ptr_b + 1));
	      b[2].v = _mm_setr_pd(*(rhs_ptr_a + 2), *(rhs_ptr_b + 2));
	      b[3].v = _mm_setr_pd(*(rhs_ptr_a + 3), *(rhs_ptr_b + 3));
	      b[0].v *= lhs_ij[0].v;
	      b[1].v *= lhs_ij[0].v;
	      b[2].v *= lhs_ij[0].v;
	      b[3].v *= lhs_ij[0].v;
	      temp0[0] += b[0];
	      temp1[0] += b[1];
	      temp2[0] += b[2];
	      temp3[0] += b[3];

rhs_ptr_a += 2*rhs.Columns();
rhs_ptr_b += 2*rhs.Columns();
	      b[0].v = _mm_setr_pd(*(rhs_ptr_a) , *(rhs_ptr_b));
	      b[1].v = _mm_setr_pd(*(rhs_ptr_a + 1), *(rhs_ptr_b + 1));
	      b[2].v = _mm_setr_pd(*(rhs_ptr_a + 2), *(rhs_ptr_b + 2));
	      b[3].v = _mm_setr_pd(*(rhs_ptr_a + 3), *(rhs_ptr_b + 3));
		
	      	      b[0].v *= lhs_ij[1].v;
	      b[1].v *= lhs_ij[1].v;
	      b[2].v *= lhs_ij[1].v;
	      b[3].v *= lhs_ij[1].v;
	      temp0[1] += b[0];
	      temp1[1] += b[1];
	      temp2[1] += b[2];
	      temp3[1] += b[3];

rhs_ptr_a += 2*rhs.Columns();
rhs_ptr_b += 2*rhs.Columns();
	      b[0].v = _mm_setr_pd(*(rhs_ptr_a) , *(rhs_ptr_b));
	      b[1].v = _mm_setr_pd(*(rhs_ptr_a + 1), *(rhs_ptr_b + 1));
	      b[2].v = _mm_setr_pd(*(rhs_ptr_a + 2), *(rhs_ptr_b + 2));
	      b[3].v = _mm_setr_pd(*(rhs_ptr_a + 3), *(rhs_ptr_b + 3));
		
	      b[0].v *= lhs_ij[2].v;
	      b[1].v *= lhs_ij[2].v;
	      b[2].v *= lhs_ij[2].v;
	      b[3].v *= lhs_ij[2].v;
	      temp0[2] += b[0];
	      temp1[2] += b[1];
	      temp2[2] += b[2];
	      temp3[2] += b[3];

rhs_ptr_a += 2*rhs.Columns();
rhs_ptr_b += 2*rhs.Columns();
	      b[0].v = _mm_setr_pd(*(rhs_ptr_a) , *(rhs_ptr_b));
	      b[1].v = _mm_setr_pd(*(rhs_ptr_a + 1), *(rhs_ptr_b + 1));
	      b[2].v = _mm_setr_pd(*(rhs_ptr_a + 2), *(rhs_ptr_b + 2));
	      b[3].v = _mm_setr_pd(*(rhs_ptr_a + 3), *(rhs_ptr_b + 3));
	      
	      
	      b[0].v *= lhs_ij[3].v;
	      b[1].v *= lhs_ij[3].v;
	      b[2].v *= lhs_ij[3].v;
	      b[3].v *= lhs_ij[3].v;
	      temp0[3] += b[0];
	      temp1[3] += b[1];
	      temp2[3] += b[2];
	      temp3[3] += b[3];
	    }
	    temp0[0] += temp0[2];
	    temp0[1] += temp0[3];
	    
	    temp1[0] += temp1[2];
	    temp1[1] += temp1[3];
	    
	    temp2[0] += temp2[2];
	    temp2[1] += temp2[3];
	    
	    temp3[0] += temp3[2];
	    temp3[1] += temp3[3];
	    
	    double temp_2[4];
	    temp_2[0] = (temp0[0].x[0] + temp0[0].x[1]) + (temp0[1].x[0] + temp0[1].x[1]);
	    temp_2[1] = (temp1[0].x[0] + temp1[0].x[1]) + (temp1[1].x[0] + temp1[1].x[1]);
	    temp_2[2] = (temp2[0].x[0] + temp2[0].x[1]) + (temp2[1].x[0] + temp2[1].x[1]);
	    temp_2[3] = (temp3[0].x[0] + temp3[0].x[1]) + (temp3[1].x[0] + temp3[1].x[1]);
//            double temp2 = temp[0] + temp[1] + temp[2] + temp[3];
	    for (unsigned int j = maxj ; j < rhs.Rows() ; ++j)
	    {
	      double lhs_ij = lhs.GetElement(i,j);
                temp_2[0] += lhs_ij * rhs.GetElement(j, k);
		temp_2[1] += lhs_ij * rhs.GetElement(j, k+1);
		temp_2[2] += lhs_ij * rhs.GetElement(j, k+2);
		temp_2[3] += lhs_ij * rhs.GetElement(j, k+3);
	    }
            retval(i, k, temp_2[0]);
	    retval(i, k+1, temp_2[1]);
	    retval(i, k+2, temp_2[2]);
	    retval(i, k+3, temp_2[3]);
        }
        
        for ( unsigned int k = maxk; k < rhs.Columns(); ++k)
        {
            Vec2 temp[DPC/2] = {0};
	    const double* lhs_ptr = lhs.arr + i * lhs.Columns();
            for (unsigned int j = 0; j < maxj; j+=DPC, lhs_ptr += DPC)
	    {
	      Vec2 lhs_ij[4];
	      for(uint t = 0; t < DPC/2; ++t)
	        lhs_ij[t].loadpd(/*lhs.arr + i * lhs.Columns() +*/ lhs_ptr + 2*t);
	      Vec2 b;
	      b.v = _mm_set_pd(rhs.GetElement(j + 2*0 + 1,k), rhs.GetElement(j + 2*0,k));
	        lhs_ij[0].v *= b.v;
		temp[0] += lhs_ij[0];
		
		b.v = _mm_set_pd(rhs.GetElement(j + 2*1 + 1,k), rhs.GetElement(j + 2*1,k));
	        lhs_ij[1].v *= b.v;
		temp[1] += lhs_ij[1];
		
		b.v = _mm_set_pd(rhs.GetElement(j + 2*2 + 1,k), rhs.GetElement(j + 2*2,k));
	        lhs_ij[2].v *= b.v;
		temp[2] += lhs_ij[2];
		
		b.v = _mm_set_pd(rhs.GetElement(j + 2*3 + 1,k), rhs.GetElement(j + 2*3,k));
	        lhs_ij[3].v *= b.v;
		temp[3] += lhs_ij[3];
	    }
	    temp[0] += temp[2];
	    temp[1] += temp[3];
	    double temp2 = (temp[0].x[0] + temp[0].x[1]) + (temp[1].x[0] + temp[1].x[1]);
//            double temp2 = temp[0] + temp[1] + temp[2] + temp[3];
	    for (unsigned int j = maxj ; j < rhs.Rows() ; ++j)
                temp2 += lhs.GetElement(i,j) * rhs.GetElement(j,k);
            retval(i, k, temp2);
        }
        
/*        for ( unsigned int k = maxk; k < rhs.Columns(); ++k )
        {
            double temp[4] = {0};
            for (unsigned int j = 0; j < maxj; j+=4)
	    {
                temp[0] += lhs.GetElement(i,j + 0) * rhs.GetElement(j + 0,k);
		temp[1] += lhs.GetElement(i,j + 1) * rhs.GetElement(j + 1,k);
		temp[2] += lhs.GetElement(i,j + 2) * rhs.GetElement(j + 2,k);
		temp[3] += lhs.GetElement(i,j + 3) * rhs.GetElement(j + 3,k);
	    }
	    double temp2 = (temp[0] + temp[1]) + (temp[2] + temp[3]);
	    for (unsigned int j = maxj ; j < rhs.Rows() ; ++j)
                temp2 += lhs.GetElement(i,j) * rhs.GetElement(j,k);
            retval(i,k, temp2);
        }*/
    }
  }
};
#endif

template < class A , class B >
inline Matrix< typename ConfigParser<A,B, MultiplicationSizer >::RET > Multiply( const A& lhs , const B& rhs, Int2Type<false> )
{//Matrix Multiplication
 //Complexity: 2n^3 - n^2
    typedef typename TypePromoter<Scalar(A), Scalar(B)>::RET ScalarType;
    typedef Matrix< typename ConfigParser<A,B, MultiplicationSizer >::RET > ReturnType;
    MatrixCreator<ReturnType, None<ScalarType> > M;
    ReturnType RetVal = M.Create( lhs.Rows(), rhs.Columns(), None<ScalarType>(0));
    if(unlikely(lhs.Columns() != rhs.Rows()) )
        throw("Sizes don't match for multiplication!");
    Multiplication_specialization<ReturnType, A, B>::mult(RetVal, lhs, rhs);
    return RetVal;
}

//Matrix-Matrix - Multiply
template<class A, class B >
inline typename MultiplyReturnValueWrapper<MTLICCL::Matrix<A>, MTLICCL::Matrix<B> >::RET operator* (const MTLICCL::Matrix<A>& lhs, const MTLICCL::Matrix<B>& rhs) MTL_NEVER_NEEDED_INLINE;

template<class A, class B >
inline typename MultiplyReturnValueWrapper<MTLICCL::Matrix<A>, MTLICCL::Matrix<B> >::RET operator* (const MTLICCL::Matrix<A>& lhs, const MTLICCL::Matrix<B>& rhs)
{//A wrapper around Multiplication to support Scalar and Matrix arguments as lhs value
    return Multiply(lhs, rhs , Int2Type<Conversion<MTLICCL::Matrix<A>, typename MTLICCL::Matrix<B>::Elementtype >::sametype >() );
}

template<class A/*scalar*/, class B/*matrix*/ >
inline typename MultiplyReturnValueWrapper<A, MTLICCL::Matrix<B> >::RET operator* (const A& lhs, const MTLICCL::Matrix<B>& rhs) MTL_NEVER_NEEDED_INLINE;

//Scalar - Matrix
template<class A/*scalar*/, class B/*matrix*/ >
inline typename MultiplyReturnValueWrapper<A, MTLICCL::Matrix<B> >::RET operator* (const A& lhs, const MTLICCL::Matrix<B>& rhs)
{//A wrapper around Multiplication to support Scalar and Matrix arguments as lhs value
    return Multiply(lhs, rhs , Int2Type<Conversion<typename MTLICCL::Matrix<B>::Elementtype, A>::sametype >() );
}

}
#endif
