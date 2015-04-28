/***************************************************************************
 *   Copyright (C) 2006-2009 by Florian Goth   *
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


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <cstdlib>
#include <complex>

#include "Matrix.h"
#include "MatrixGenerator.h"
#include "CommonTypes.h"
#include "MemArray1D.h"

using namespace std;

typedef MTLICCL::Matrix< MTLICCL::BoundsChecker< MTLICCL::Static<MTLICCL::Config<float,3,3> > > > Mat_Static;

typedef MTLICCL::Matrix< MTLICCL::MemArray1D<MTLICCL::Config1D<double> > > Mat_1D;


template < typename Matrix_T>
float ZeilenSummenNorm(Matrix_T a)
{
    float maximum = 0;
    for (unsigned int i = 0 ; i< a.Rows() ; ++i )
    {
        float t = a(i,0);
        for (unsigned int k = 1 ; k< a.Columns() ; ++k )
            t += a(i,k);
        maximum = max(maximum,t);
    }
    return maximum;
}

template < typename Matrix_T>
float SpaltenSummenNorm(Matrix_T a)
{
    float maximum = 0;
    for (unsigned int i = 0 ; i < a.Columns() ; ++i )
    {
        float t = a(i,0);
        for (unsigned int k = 1 ; k< a.Rows() ; ++k )
            t += a(i,k);
        maximum = max(maximum,t);
    }
    return maximum;
}

template < class A >
typename A::Elementtype EuklidNorm(const A& a)
{//only usable for Vector-Types!!!!!!!
typedef typename A::Elementtype ScalarType;
    ScalarType n = a(0,0)*a(0,0);

    for (unsigned int k = 1 ; k< a.Rows() ; ++k )
        n += a(k,0) * a(k,0);

    return sqrt(n);
}


template <
class VektorTyp , //Der Typ der rechten Seite...
class Abstiegsmethode // THe Method to use for Iteration
>
pair< float, unsigned int > Abstiegsverfahren (VektorTyp X0 , unsigned int Imax , float epsilon , Abstiegsmethode Method )
{
    VektorTyp xvorher(X0.Rows() , X0.Columns() ) ;
    unsigned int t =0;
    VektorTyp xt = X0;
    VektorTyp rt = Method.Initial_rt(X0);
    float alphat = Method.Initial_alphat();
    //     for ( t = 1 ; t < Imax /*|| (MaxNorm(xt + (-1)*xvorher) < epsilon)*/; ++t)
    //     {
    //         alphat = Method.alphat(t);
    //         rt = Method.rt(t);
    //         xt = xt + alphat * rt ;
    //     }
    return make_pair(1,1);
}

template
<
class A // Der MatrixTyp
//FIXME: The Vectortype should be derived from the MatrixType
>
class GradientenVerfahren
{
private:
    typedef typename A::Elementtype elementtype;
    A Mat;
    A rechteSeite;
public:
    GradientenVerfahren(A l , A r) : Mat(l) , rechteSeite(r)
    {}
    elementtype Initial_alphat()
    {}
    A Initial_rt(A x0 )
    {
        return (Mat * x0)  - rechteSeite ;
    }
    elementtype alphat()
    {}
    A rt(A x0 )
    {}
}
;


class Jacobi
{
private:
public:
    template< class Type_A , class Type_B>
    Jacobi(Type_A A,Type_B B )
    {}
}
;

class PotenzMethode
{
private:
    Mat_Dynamic* It;
    Mat_Dynamic b;
public:
    template< class Type_A>
    PotenzMethode(Type_A& Alpha ) : b(Alpha.Columns(), 1,MTLICCL::None<float>(0.0) )
    {
        It = &Alpha;
    }
    //MTLICCL::Matrix< MTLICCL::BoundsChecker<MTLICCL::Dynamic< MTLICCL::Config_Dynamic<float> > > > B(unsigned int) { return *It;}
    Mat_Dynamic c(unsigned int)
    {
        return b; /*return an empty vector*/
    }

};

template< class Iterationsverfahren, class VektorTyp>
pair<float , unsigned int > Iterate(Iterationsverfahren Iterationmethod , VektorTyp z_0,  unsigned int t_Max , float epsilon)
{
    unsigned int t = 0;
    VektorTyp x_alt = z_0;
    while ( (t < t_Max) /* && (epsilon)*/)
    {
        Mat_Dynamic IterationsMatrix = Iterationmethod.B(t);
        Mat_Dynamic IterationsVektor = Iterationmethod.c(t);
        VektorTyp x_neu =  IterationsMatrix  * x_alt;
        x_neu = x_neu + IterationsVektor;
        t++;
        x_alt = x_neu;
    }
    return make_pair( 1 ,1 );
}

template < typename VektorTyp>
struct PotenzMethodenRetVal
{
    PotenzMethodenRetVal( typename VektorTyp::Elementtype A , VektorTyp B , unsigned int C) : lambda(A) , z(B) , I(C)
    {}
    typename VektorTyp::Elementtype lambda;
    VektorTyp z;
    unsigned int I;
};

template< typename Matrix_type >
typename Matrix_type::Elementtype RayleighQuotient( const Matrix_type& A , const Matrix_type& z )
{//z muss in der Euklid-Norm normiert sein!!
    typedef typename Matrix_type::Elementtype ScalarType;
    Matrix_type t(A*z);
    ScalarType ret = 0;
    for ( unsigned int k = 0 ; k < z.Rows() ; ++k )
        ret += t(k,0) * z(k,0);
    return ret;
}

template < class Matrix_type, class Vectortype >
PotenzMethodenRetVal<Vectortype> PotenzMethode( const Matrix_type& U , const Vectortype& z0 , unsigned int Imax , float epsilon)
{//z0 muss auf 1 normiert sein!!
    unsigned int t = 0;
    typedef typename Matrix_type::Elementtype ScalarType;
    ScalarType lambdat = 0;

    Matrix_type zalt ( z0 ) ;
    Matrix_type zneu ( z0 );

    while ( (t< Imax) &&( EuklidNorm( U * zneu - lambdat *  zneu ) > epsilon ) )
    {
        Matrix_type zneuschlange ( U*zalt );//Iterationsverfahren -> Potenzbildung
        ScalarType alpha = 1.0 / (EuklidNorm(zneuschlange) );//Normierungsfaktor
        zneu = alpha * zneuschlange; //Normierung
        lambdat = RayleighQuotient( U , zneu);// Berechnung der t-ten Eigenwert-Iterierten
        t+=1;
        //Zuweisung der neuen auf die Alten Vektoren
        zalt = zneu;
    }
    return PotenzMethodenRetVal<Vectortype> ( lambdat , zneu , t );
}

template < typename M >
pair<M , typename M::Elementtype * > QR_Zwischenschritt ( const M& A )
{
    typedef typename M::Elementtype ScalarType;
    unsigned int n = A.Rows();
    ScalarType* d = new ScalarType[2*n];//The Diagonalelement that is used later to reconstruct the Matrices
    //In the first n Elements are the diagonalelements stored, then come the betas to restore Q
    M R( A );
    for ( unsigned int j = 0 ; j < n ; ++j )
    {
        ScalarType sigma = 0;
        for ( unsigned int i = j ; i < n ; ++i)
            sigma += R(i,j)*R(i,j) ;
        if( sigma == 0.0 )
            throw ("Singulär!!");

        ScalarType s =   (R(j,j) < 0) ? sqrt(sigma) : -sqrt(sigma);
        d[j] = s;
        ScalarType beta = 1.0 / (s*R(j,j) - sigma);
        d[n+j] = beta;//Save it to ease Reconstruction
        R( j , j , R( j , j ) - s  );//Zuweisung
        for ( unsigned int k = j+1 ; k < n ; ++k)
        {
            ScalarType sum = 0;
            for ( unsigned int i = j ; i < n ; ++i )
                sum += R(i,j) * R(i,k);
            sum *= beta;
            for ( unsigned int i = j ; i < n ; ++i )
                R( i , k , R(i,k) + R(i,j) * sum ); //Zuweisung
        }
    }
    return make_pair(R , d );
}

template < typename M >
pair< M , M > QR_Zerlegung(M A )
{
    typedef typename M::Elementtype ScalarType;

    pair< M , typename M::Elementtype * > Ret( QR_Zwischenschritt(A)) ;
    //    cout<<Ret.first<<endl;
    const unsigned int n = Ret.first.Rows();
    M R(n, n, MTLICCL::None<ScalarType>(0.0f));
    for ( unsigned int i = 0 ; i < n ; ++i )
        for ( unsigned int k = i ; k < n; ++k )
            R( i , k , Ret.first(i,k) );
    for ( unsigned int i = 0 ;i < n ; ++i)
        R(i,i, Ret.second[i] );
    //R is now ready
    //cout<<R<<endl;
    //Get Q
    M Q(n, n, MTLICCL::Identity<ScalarType>(1 , 0 ));
    ScalarType* beta = &(Ret.second[n]);
    for ( int j = n-1  ; j >=0 ; --j ) // Diese Schleife multipliziert alle Matrizen P auf um Q zu erhalten
    {
        M P( n, n, MTLICCL::Identity<ScalarType>(1 , 0 ));//Die einzelne Householder - Transformationsmatrix
        const unsigned int Udim = n - j ;
        M u (Udim , 1 , MTLICCL::None<ScalarType> ( 0 ));
        for ( unsigned int cnt = 0 ; cnt < Udim ; ++cnt )
            u( cnt , 0 , Ret.first(j + cnt ,j) );
//u = 1.0 / EuklidNorm(u) * u;
        M v( ~u);

        M PTilde ( Udim , Udim ,MTLICCL::Identity<ScalarType>(1.0 ,0.0f) );
        PTilde = PTilde + (-beta[j]) * ( u * v /* == ~u*/ );
        for ( int o = 0 ; o < Udim ; ++o)
            for ( int p = 0 ; p < Udim ; ++p)
                //Now Copy the small Matrix PTilde into the larger one P
                P( j + o , j + p , PTilde(o , p ));
        Q = Q*P;
    }//Got Q

    delete [] Ret.second;
    return make_pair( Q , R );
}
/*
template <typename M >
M QR_Verfahren(const M& A, unsigned int t )
{
    M It ( A );
    for ( unsigned int j = 0 ; j < t ; ++j )
    {
        pair < M , M > QR = QR_Zerlegung(A);
        It = QR.second * QR.first;
    }
    return It;
}*/

int main(int argc, char *argv[])
{
    //MATRIX_GENERATOR< matrix<float> >::RET A = 3.9;//Should Resolve to a static 3x3 matrix of floats
    Mat_Static Static;
    Mat_Dynamic Dynamic(3,3, MTLICCL::None<double>(0.0f));
    Mat_Dynamic Dyn4x4(4,4,MTLICCL::None<double>(0.0f));
    Mat_Dynamic Vector3(3,1, MTLICCL::None<double>(0.0));
    Mat_Dynamic Vector4(4,1, MTLICCL::None<double>(0.0));
    Mat_Dynamic Dyn5x5(5,5, MTLICCL::TriGen<double>(1,2,1) );
//    MTLICCL::Matrix< MTLICCL::Dynamic< MTLICCL::Config_Dynamic<double> > > Dyn6x6(6,6, MTLICCL::None<float>(0.0f));
    Mat_1D Dyn6x6(6, 6, MTLICCL::None<float>(0.0f));
MTLICCL::Matrix< MTLICCL::Dynamic< MTLICCL::Config_Dynamic<float> > > DynLU3x3(3,3, MTLICCL::None<float>(0.0f));
    MTLICCL::Matrix< MTLICCL::BoundsChecker< MTLICCL::Static<MTLICCL::Config<float,4,4> > > > StaticBounds(MTLICCL::None<double>(0.0));
MTLICCL::Matrix< MTLICCL::BoundsChecker< MTLICCL::Static<MTLICCL::Config<float,4,4> > > > StaticIdentity(MTLICCL::Identity<double>(1.0, 0.0));
MTLICCL::Matrix< MTLICCL::BoundsChecker< MTLICCL::Dynamic<MTLICCL::Config<float,4,4> > > > DynamicIdentity(3,3, MTLICCL::Identity<double>(1.0, 0.0));
MTLICCL::Matrix<  MTLICCL::Static<MTLICCL::Config<float,6,1> >  > StaticVecLU(MTLICCL::None<double>(0.0));

Matrix4x4 Inverstest(MTLICCL::TriGen<float>(1.0,-4.0,1.0));
MTLICCL::Matrix< MTLICCL::Static<MTLICCL::Config<float,1,3> > > StaticVec(MTLICCL::None<double>(0.0));
StaticVec(0,0,4);
Vec4 vec;
vec(0,0, 1);
vec(1,0, 2);
vec(2,0, 3);
vec(3,0, 4);
//Initialize Vectors and Matrices
StaticBounds(0,1,3);
    Vector3(0,0,1);
    Vector4(0,0,1);

    Dynamic(0,0,2);
    Dynamic(0,1,1);
    Dynamic(0,2,0);
    Dynamic(1,0,1);
    Dynamic(1,1,6);
    Dynamic(1,2,1);
    Dynamic(2,0,0);
    Dynamic(2,1,1);
    Dynamic(2,2,10);

    Static(0,0,2);
    Static(0,1,1);
    Static(0,2,0);
    Static(1,0,1);
    Static(1,1,6);
    Static(1,2,1);
    Static(2,0,0);
    Static(2,1,1);
    Static(2,2,10);

    Dyn4x4(0,0,4);
    Dyn4x4(0,1,4);
    Dyn4x4(0,2,4);
    Dyn4x4(0,3,0);
    Dyn4x4(1,0,16);
    Dyn4x4(1,1,5);
    Dyn4x4(1,2,1);
    Dyn4x4(1,3,2);
    Dyn4x4(2,0,4);
    Dyn4x4(2,1,1);
    Dyn4x4(2,2,2);
    Dyn4x4(2,3,-1);
    Dyn4x4(3,0,0);
    Dyn4x4(3,1,2);
    Dyn4x4(3,2,1);
    Dyn4x4(3,3,6);
    Dyn6x6(0,0,1);Dyn6x6(0,1,3.3);Dyn6x6(0,2,2);
    Dyn6x6(1,0,2);Dyn6x6(1,1,4.5);Dyn6x6(1,2,1);
    Dyn6x6(2,0,-0.27);Dyn6x6(2,1,0);Dyn6x6(2,2,3);
    Dyn6x6(3,0,2);Dyn6x6(3,1,-0.1);Dyn6x6(3,2,2);Dyn6x6(3,3,4);
    Dyn6x6(4,0,2);Dyn6x6(4,1,2);Dyn6x6(4,2,2);Dyn6x6(4,3,2);Dyn6x6(4,4,5);
    Dyn6x6(5,0,-0.5);Dyn6x6(5,1,2);Dyn6x6(5,2,2);Dyn6x6(5,3,2);Dyn6x6(5,4,2);Dyn6x6(5,5,6);
    //Abstiegsverfahren(Vec0, 10u , 0.1f , GradientenVerfahren<Mat_Dynamic>(Dynamic , Vec0) );
    //Iterate( PotenzMethode(Dynamic) ,Vec0 , 10 , 0.01);
/*Mat_Dynamic AQ((QR_Zerlegung(Dyn5x5 )).first);
Mat_Dynamic AR ((QR_Zerlegung(Dyn5x5 )).second);
    cout<<"Q: "<<endl <<AQ<<endl;
    cout<<"R: "<<  endl<<AR<<endl;
 cout<<(AQ * AR)<<endl;*/

//    Mat_Dynamic Ret15 ( QR_Verfahren( Dyn5x5 , 15));
//    cout<<Ret15<<endl;
//for ( unsigned int k = 0 ; k < Ret15.Rows() ; ++ k )
//cout<<"Eigenwert : "<<k<<" "<<Ret15(k,k)<<endl;

//Tests:
cout<<"Testing Copy-Ctor's"<<endl;
cout<<Dyn4x4<<endl;
Mat_Dynamic Dyn4x4B(Dyn4x4);
cout<<Dyn4x4B<<endl;
Dyn4x4B(2,2,1000);
cout<<"Testing Assignment Operators"<<endl;
cout<<Dyn4x4<<endl;
Dyn4x4B = Dyn4x4;
cout<<Dyn4x4B<<endl;
cout<<"Testing dynamic Matrix Vector Multiplication"<<endl;
cout<<Dynamic<<endl;
cout<<Vector3<<endl;
cout<<Dynamic * Vector3<<endl;//Dynamic - Dynamic - Multiplication
cout<<"Testing static Matrix Vector Multiplication"<<endl;
cout<<Inverstest<<endl;
cout<<vec<<endl;
cout<<Inverstest * vec<<endl;//static - static - Multiplication
/*cout<<"Testing static Matrix-Matrix Multiplication"<<endl;
cout<<StaticBounds<<endl;
cout<<StaticIdentity<<endl;
cout<<StaticBounds * StaticIdentity<<endl;//Static - Static - Mutiplication
cout<<"Testing Multiplication between a Static Matrix and a dynamic Vector"<<endl;
cout<<Static<<endl;
cout<<Vector3<<endl;
cout<<Static * Vector3<<endl;//Static-Dynamic-Multiplication(Matrix - Vektor)
cout<<"Testing Multiplication between a dynamic Matrix and a static Matrix"<<endl;
cout<<Dynamic<<endl;
cout<<Static<<endl;
cout<<Dynamic * Static<<endl;//Dynamic - Static-Multiplication(Matrix-Matrix)
cout<<"Testing dynamic Matrix Transposition"<<endl;
cout<<Dyn4x4<<endl;
cout<<~Dyn4x4<<endl;
cout<<"Testing static Matrix Transposition"<<endl;
cout<<StaticBounds<<endl;
cout<<~StaticBounds<<endl;
cout<<"Testing dynamic Vector Transposition"<<endl;
cout<<Vector3<<endl;
cout<<~Vector3<<endl;
cout<<"Testing static Vector Transposition"<<endl;
cout<<StaticVec<<endl;
cout<<~StaticVec<<endl;
cout<<"Testing Euclidean Product"<<endl;
cout<<(~Vector3)*(~StaticVec)<<endl;
cout<<"Testing more complicated things: x(dyn) * A(stat) * x(stat)"<<endl;
cout<<(~Vector3)*(Static * (~StaticVec))<<endl;
cout<<"Testing static - Static Matrix Addition"<<endl;
cout<<StaticBounds<<endl;
cout<<StaticIdentity<<endl;
cout<<StaticIdentity+StaticBounds<<endl;
cout<<"Testing dynamic-dynamic Matrix Addition"<<endl;
cout<<Dyn4x4<<endl;
cout<<Dyn4x4+Dyn4x4<<endl;
cout<<"Testing static -static Matrix Subtraction"<<endl;
cout<<StaticBounds<<endl;
cout<<StaticIdentity<<endl;
cout<<StaticIdentity - StaticBounds<<endl;
cout<<"Testing dynamic-dynamic Matrix Subtraction"<<endl;
cout<<Dyn4x4<<endl;
cout<<Dyn4x4 - Dyn4x4<<endl;
cout<<"Testing Dynamic - Static Matrix Subtraction"<<endl;
cout<<Dynamic<<endl;
cout<<Static<<endl;
cout<<Dynamic - Static<<endl;

cout<<"Testing Scalar(5.0) * Matrix(dyn) Multiplication"<<endl;
cout<<Dyn4x4<<endl;
cout<<5.0 * Dyn4x4<<endl;*/
/*cout<<"Testing Matrix(dyn) * Scalar(5.0) Multiplication"<<endl;
cout<<Dyn4x4<<endl;
cout<<Dyn4x4 * 5.0f<<endl;*/
cout<<"Testing Scalar(5.0) * Matrix(stat) Multiplication"<<endl;
cout<<Static<<endl;
cout<<5.0f * Static<<endl;//FIXME: The f is necessary, because Static is a float Matrix and 5.0 is a double
//Tests that should fail:
//cout<<Static * StaticBounds<<endl;
try {
cout<<StrassensMultiply(Dyn4x4, Dyn4x4)<<endl;
}
catch(const char* e)
{
cout<<e<<endl;
}
cout<<"Something more fancy: Potenzmethode"<<endl;
    //PotenzMethodenRetVal< MTLICCL::Matrix< MTLICCL::BoundsChecker< MTLICCL::Static<Config<float,1,3> > > > > Erg3 = PotenzMethode(Static, StaticVec , 10 , 0.001);
    PotenzMethodenRetVal<Mat_Dynamic> Erg4 = PotenzMethode(Dyn4x4, Vector4 , 10 , 0.001);

    cout<<"A = " <<endl<<Dyn4x4<<endl;
    cout<<"lambda_Max: "<<Erg4.lambda<<endl;
    cout<<"Eigenvektor: "<<endl<<Erg4.z<<endl;
    cout<<"benötigte Iterationen : "<<Erg4.I<<endl;

    /*cout<<"B = "<<endl<<Dynamic<<endl;
    cout<<"lambda_Max: "<<Erg3.lambda<<endl;
    cout<<"Eigenvektor: "<<endl<<Erg3.z<<endl;
    cout<<"benoetigte Iterationen : "<<Erg3.I<<endl;*/

cout<<"Testing Inversion"<<endl;
cout<<Inverstest<<endl;
cout<<inverse(Inverstest) * Inverstest<<endl;
cout<<"Testing Determinant"<<endl;
cout<<det(Dyn6x6)<<endl;
cout<<"testing LU Decomposition"<<endl;
cout<<Dyn6x6<<endl;
Mat_1D invDyn6x6(inverse(Dyn6x6));
unsigned int* idx = new unsigned int[Dyn6x6.Rows()];
int sign = ludecompose(Dyn6x6, idx);
cout<<"Sign: "<<sign<<endl;
cout<<Dyn6x6<<endl;
cout<<"testing Inversion"<<endl;
StaticVecLU(0,0,0);StaticVecLU(1,0,0);StaticVecLU(2,0,0);StaticVecLU(3,0,0);StaticVecLU(4,0,0);StaticVecLU(5,0,1);
cout<<invDyn6x6<<endl;
cout<<( invDyn6x6 * StaticVecLU)<<endl;
cout<<"Testing solving linear systems"<<endl;
double tv[] = {1.0, 1.0, 2.0,1.0, 1.0, 1.0};
double* testvec = tv;
lubacksubstitute(Dyn6x6, idx, StaticVecLU);
for(unsigned int k = 0; k < 6; ++k)
cout<<StaticVecLU(k,0)<<endl;
cout<<"Testing Determinant Update formula:"<<endl;
double u2[6] = {1,0,0,0,0,0};
double v1[6] = {1,0,0,0,0,0};
cout<<det(Dyn6x6, idx, -571.32, v1, u2, 1)<<endl;
delete [] idx;
std::cout<<"Proceeding to Eigensystem Computations"<<std::endl;
uint dim = 6;
Mat_Dynamic esstest(dim,dim);
for(uint k = 0; k < dim; ++k)
{
  esstest(k,k) = 2.0;
  if(k+1<dim)
  {
  esstest(k,k+1) = 1.0;
  esstest(k+1,k) = 1.0;
  }
  if(k+2 < dim)
  {
  esstest(k,k+2) = 0.5;
  esstest(k+2,k) = 0.5;
  }
  if(k+3 < dim)
  {
  esstest(k,k+3) = 1.0;
  esstest(k+3,k) = 1.0;
  }
}
std::cout<<esstest<<std::endl;
MTLICCL::EigenSystemSolver<Mat_Dynamic> ess(esstest);
ess.calculateEVs();
Mat_Dynamic z(ess.tridiagonalize());
std::vector<MTLICCL::SolutionPair<double> > esystem = ess.eigensystem();
for(uint k = 0; k < esystem.size(); ++k)
  esystem[k].print();
Mat_Dynamic evscplx(dim, dim);
 for(uint k = 0; k < dim; ++k)
   for(uint j = 0; j < dim; ++j)
   {
     evscplx(j, k) = esystem[k].evector[j];
   }
   std::cout<<z*evscplx<<std::endl;
//exit(-1);

/*
CmplxMat testmat(dim, dim);
for(uint k = 0; k < dim; ++k)
{
  testmat(k,k) = 2.0;
  if(k+1<dim)
  {
  testmat(k,k+1) = 1.0;
  testmat(k+1,k) = 1.0;
  }
  if(k+2 < dim)
  {
  testmat(k,k+2) = std::complex<double>(0.0, 1.0)*0.5;
  testmat(k+2,k) = 0.5*std::complex<double>(0.0, -1.0);
  }
  if(k+3 < dim)
  {
  testmat(k,k+3) = std::complex<double>(1.0,+ 1.0);
  testmat(k+3,k) = std::complex<double>(1.0, -1.0);
  }
}
std::cout<<testmat<<std::endl;
MTLICCL::EigenSystemSolver<CmplxMat> ess(testmat);
MTLICCL::EigenSystemSolver<CmplxMat> ess2(testmat);
ess2.calculateEVs();
CmplxMat zcplx(ess2.tridiagonalize());
//CmplxMat zcplx(ess.tridiagonalize());
std::vector<MTLICCL::SolutionPair<double> > esystemcplx = ess2.eigensystem();
 CmplxMat evscplx(zcplx.Rows(), zcplx.Rows());
 for(uint k = 0; k < zcplx.Rows(); ++k)
   for(uint j = 0; j < zcplx.Rows(); ++j)
   {
     evscplx(j, k) = esystemcplx[k].evector[j];
   }
   std::cout<<zcplx*evscplx<<std::endl;
 CmplxMat zhc(~(zcplx*evscplx));
 for(uint k = 0; k < zhc.Rows(); ++k)
   for(uint i = 0; i < zhc.Columns(); ++i)
     zhc(k,i) = conj(zhc(k,i));
   CmplxMat id(zhc*(zcplx*evscplx));
 for(uint k = 0; k < id.Rows(); ++k)
 {
     double sum = 0.0;
   for(uint i = 0; i < id.Rows(); ++i){
     if(i!=k)
     sum += abs(id(i,k));
   }
   std::cout<<"Zeilensumme "<<k<<"    : "<<scientific<<sum<<" + "<<abs(id(k,k))<<std::endl;
 }
 std::cout.precision(12);
 ofstream file("evs1000.txt");
 file.precision(12);
 for(uint k = 0; k < esystemcplx.size(); ++k)
   file<<esystemcplx[k].evalue<<", ";*/
 std::ofstream file("energies.txt");
 std::ofstream summen("zeilensummen.txt");
 uint Nx= 512;
 uint Nb = 40;
 double lambda = 0.2;
 uint bsize = 4;
 uint ksize = Nb * bsize;
 uint size = Nx * ksize;
 std::cout<<atoi(argv[1])<< " " <<atoi(argv[2])<<std::endl;
for(uint kidx = atoi(argv[1]); kidx < atoi(argv[2]); ++kidx)
{
  double k = kidx*2.0*M_PI/Nx;
  summen<<kidx<<std::endl;
  double cosk = 2*cos(k/2);
  cosk = cosk * cosk;
  double sink = sin(k);
  std::complex<double> expk = exp(std::complex<double>(0.0,k));
for(uint s = 0; s < 2; ++s)
{
  CmplxMat hzero_k_sigma(ksize, ksize, MTLICCL::Identity<double>(0.0, 0.0));
  CmplxMat hlambda_k_sigma(ksize, ksize, MTLICCL::Identity<double>(0.0, 0.0));
  for(uint b = 0; b < Nb; ++b)
  {
      hzero_k_sigma(b*bsize + 0, b *bsize + 1) += cosk;
      hzero_k_sigma(b*bsize + 1, b *bsize + 0) += cosk;
      
      hzero_k_sigma(b*bsize + 2, b *bsize + 3) += cosk;
      hzero_k_sigma(b*bsize + 3, b *bsize + 2) += cosk;
      
      hzero_k_sigma(b*bsize + 2, b *bsize + 1) += 2;
      hzero_k_sigma(b*bsize + 1, b *bsize + 2) += 2;
      
      if(b != 0)
      {
	hzero_k_sigma(b*bsize + 0, (b-1) *bsize + 3) += 1;
	hzero_k_sigma((b-1)*bsize + 3, b *bsize + 0) += 1;
      }
      if(b!=(Nb -1))
      {
	hzero_k_sigma(b*bsize + 3, (b+1) *bsize + 0) += 1;
	hzero_k_sigma((b+1)*bsize + 0, b *bsize + 3) += 1;
      }
      int sigma = 2*s -1;
      hlambda_k_sigma(b * bsize + 0, b * bsize + 0) += -2.0*sink*sigma;
      hlambda_k_sigma(b * bsize + 1, b * bsize + 1) += 2.0*sink*sigma;
      hlambda_k_sigma(b * bsize + 2, b * bsize + 2) += -2.0*sink*sigma;
      hlambda_k_sigma(b * bsize + 3, b * bsize + 3) += 2.0*sink*sigma;
      
      hlambda_k_sigma(b * bsize + 0, b * bsize + 2) += static_cast<double>(sigma) * std::complex<double>(0.0, 1.0)*(1.0 - expk);
      hlambda_k_sigma(b * bsize + 2, b * bsize + 0) += static_cast<double>(sigma) * conj(std::complex<double>(0.0, 1.0)*(1.0 - expk));
      
      hlambda_k_sigma(b * bsize + 1, b * bsize + 3) += static_cast<double>(sigma) * std::complex<double>(0.0, -1.0)*(expk - 1.0);
      hlambda_k_sigma(b * bsize + 3, b * bsize + 1) += static_cast<double>(sigma) * conj(std::complex<double>(0.0, -1.0)*(expk - 1.0));
      
      if(b != 0)
      {
	hlambda_k_sigma(b * bsize + 0, (b-1) * bsize + 2) += std::complex<double>(0.0, sigma);
	hlambda_k_sigma((b-1) * bsize + 2, b * bsize + 0) += std::complex<double>(0.0, -sigma);
	
	hlambda_k_sigma(b * bsize + 1, (b-1) * bsize + 3) += std::complex<double>(0.0, -sigma) * conj(expk);
	hlambda_k_sigma((b-1) * bsize + 3, b * bsize + 1) += std::complex<double>(0.0, sigma) * expk;
      }
      if(b != (Nb-1) )
      {
	hlambda_k_sigma(b * bsize + 2, (b+1) * bsize + 0) += std::complex<double>(0.0, sigma) * conj(expk);
	hlambda_k_sigma((b+1) * bsize + 0, b * bsize + 2) += std::complex<double>(0.0, -sigma) * expk;
	
	hlambda_k_sigma(b * bsize + 3, (b+1) * bsize + 1) += std::complex<double>(0.0, -sigma);
	hlambda_k_sigma((b+1) * bsize + 1, b * bsize + 3) += std::complex<double>(0.0, sigma);
      }
    }
//      std::cout<<(hlambda_k_sigma-hzero_k_sigma)<<std::endl;
  CmplxMat fullmat = std::complex<double>(lambda,0.0)*hlambda_k_sigma - hzero_k_sigma;
  MTLICCL::EigenSystemSolver<CmplxMat> ess2(fullmat);
  ess2.calculateEVs();
  CmplxMat zcplx(ess2.tridiagonalize());
  std::vector<MTLICCL::SolutionPair<double> > esystemcplx = ess2.eigensystem();
  //summen<<" k = "<<k<<" "<<(ess2.testeigensystem(esystemcplx, zcplx)? "good": "false")<<std::endl;
//   for(uint j = 0; j < esystemcplx.size(); ++j)
//     file<<k<<" "<<esystemcplx[j].evalue<<std::endl;
//   
//  CmplxMat evscplx(zcplx.Rows(), zcplx.Rows());
//  for(uint i = 0; i < zcplx.Rows(); ++i)
//    for(uint j = 0; j < zcplx.Rows(); ++j)
//    {
//      evscplx(j, i) = esystemcplx[i].evector[j];
//    }
// //   std::cout<<zcplx*evscplx<<std::endl;
//  CmplxMat zhc(~(zcplx*evscplx));
//  for(uint j = 0; j < zhc.Rows(); ++j)
//    for(uint i = 0; i < zhc.Columns(); ++i)
//      zhc(j,i) = conj(zhc(j,i));
//    CmplxMat id(zhc*(zcplx*evscplx));
//  std::cout<<id<<std::endl;
//  for(uint j = 0; j < id.Rows(); ++j)
//  {
//      double sum = 0.0;
//    for(uint i = 0; i < id.Rows(); ++i){
//      if(i!=j)
//      sum += abs(id(i,j));
//    }
//     std::cout.precision(12);
//    summen<<"Zeilensumme "<<j<<"    : "<<scientific<<sum<<" + "<<abs(id(k,k))<<std::endl;
//  }
  }
}/*
CmplxMat test2(3,3);
test2(0,0) = 1;
test2(0,1) = 3;
test2(1,0) = 3;
test2(1,1) = 4;
test2(2,2) = 3.14159;
 MTLICCL::EigenSystemSolver<CmplxMat> ess2(test2);
  ess2.calculateEVs();
  CmplxMat zcplx(ess2.tridiagonalize());
  std::vector<MTLICCL::SolutionPair<double> > esystemcplx = ess2.eigensystem();
  for(uint j = 0; j < esystemcplx.size(); ++j)
    esystemcplx[j].print();*/
    return EXIT_SUCCESS;
}
