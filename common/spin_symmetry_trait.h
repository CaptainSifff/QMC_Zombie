#ifndef SPIN_SYMMETRY_TRAIT_H
#define SPIN_SYMMETRY_TRAIT_H
#include <stdexcept>
#include <utility>
#include <complex>
#include "MTL/Matrix.h"
#include "MTL/Matrix_Header.h"
#include "MTL/common.h"
#include "Alpha.h"

template <typename T>
inline T dppd(T a, T b, T c, T d)
{
    return a*c + b*d;
}

//we need gcc greater 3.1
#if( GCC_VERSION > GCC_VER(3,1,0))
#if defined(__SSE4_1__)
template <>
inline double dppd<double>(double a, double b, double c, double d)
{
    double result1;
    asm(       "movsd %[a1],%%xmm0  \n\t" //move
               "movhpd %[a2],%%xmm0  \n\t"
               "movsd %[b1],%%xmm1  \n\t"
               "movhpd %[b2],%%xmm1  \n\t"
               "dppd $49,%%xmm0,%%xmm1 \n\t" //SSE 4.1
               "movsd %%xmm1, %[res1]  \n\t"
           : [res1] "=mx" (result1)//output operands
                       : [a1] "mx" (a), [a2] "m" (b), [b1] "mx" (c), [b2] "m" (d)//input operands
                       : "xmm0", "xmm1"//list of clobbered registers
               );
    return result1;
}
#endif
#endif

/**
This is a Generator for Flo's MTL that fills a Matrix with the value of the free greensfunction
*/
template<class Configuration, class GreensFunction>
class WickMatrixFiller
{
private:
public:
    const Configuration& configuration;///< here we store a reference to the configuration with which we fill the matrix
    typedef typename GreensFunction::FreeGreensFunctionReturnValueType GFRetVal;///< a Typedef for the type of the return value of the free greensfunction
    inline WickMatrixFiller(const Configuration& c) : configuration(c)
    {}
    /**
    The actual operator that is used by the Matrix Library
    @param i the line in the matrix
    @param k the row in the Matrix
    @return the value of the free Greensfunction for the corresponding quantum numbers of the vertices
    */
    inline GFRetVal operator() (unsigned int i, unsigned int k ) const
    {
        const typename Configuration::value_type& v1 = *(configuration.begin() + i);
        const typename Configuration::value_type& v2 = *(configuration.begin() + k);
        GFRetVal t = GreensFunction::eval(v1, v2);
        return t;
    }
};

template <class MatTypeT>
class MatrixContainer_Common
{
public:
    typedef MatTypeT MatType;
    typedef typename MatTypeT::Elementtype DetType;
    static inline DetType rank1Add_Det(const DetType& mnn, const MatType& mat, const MatType& v1, const MatType& u2)
    {
        return mnn - (v1 * (mat * u2))(0,0);
    }
    static inline DetType det_of_diagonal_change(int idx, DetType g, const MatType& mat)
    {
      return static_cast<DetType>(1.0) + g * mat(idx, idx);
    }
    static inline void matrix_of_diagonal_change(uint pos, DetType r, DetType g, MatType& dst, const MatType& src)
    {
      DetType fac = g/r;
      dst = src;
      for(uint i = 0; i < src.Rows(); ++i)
	for(uint k = 0; k < src.Columns(); ++k)
	  dst(i, k) -= fac * src(i, pos) * src(pos, k);
	return;
    }
private:
};

template <class Configuration, int Has_Spin_symmetry = Configuration::Has_Spin_Symmetry>
class Matrix_Container;

template <class Configuration>
class Matrix_Container<Configuration, true> : private MatrixContainer_Common<typename Configuration::MatType>
{
public:
    typedef typename Configuration::MatType::Elementtype DetType;///< a typedef for the type of the elements and the determinant
    typedef typename Configuration::MatType MatType;///< the type of the Matrix that we use
    typedef typename Configuration::FPType FPType;
    typedef std::pair<DetType, DetType> RatioType;
    MatType up;
    MatType down;
    typename MatType::value_type upratio;
    typename MatType::value_type downratio;
    typename MatType::value_type getRatio() const {
        return upratio * downratio;
    }
    void reset()
    {
        up = MatType(0, 0, MTLICCL::None<DetType>(0.0));
        down = MatType(0, 0, MTLICCL::None<DetType>(0.0));
    }
    Matrix_Container(int size) : up(size, size, MTLICCL::None<DetType>(0.0)), down(size, size, MTLICCL::None<DetType>(0.0)), upratio(1.0), downratio(1.0) {}
    Matrix_Container() : upratio(1.0), downratio(1.0) {}
    inline DetType operator()(uint i, uint k, SPINS s1, SPINS s2) const
    {
      if(s1 != s2) return 0.0;
      else if(s1 == UP) return up(i,k);
	else
	  return down(i,k);
    }
    template<class VecType>
    void multiplyVectorbyConfiguration_left(VecType& v) const
    {
      VecType ret(v.Rows(), v.Columns());
      for(uint i = 0; i < v.Rows()/2; ++i)
      {
	DetType sumup = 0;
	DetType sumdown = 0;
	const uint len = up.Columns();
	for(uint k = 0; k < len; ++k)
	{
	  sumup += up(i, k) * v(k, 0);
	  sumdown += down(i, k) * v(len + k, 0);
	}
	ret(i, 0) = sumup;
	ret(len + i, 0) = sumdown;
      }
      v = ret;
    }
    template<class VecType>
    void multiplyVectorbyConfiguration_right(VecType& v) const
    {
      VecType ret(v.Rows(), v.Columns());
      for(uint i = 0; i < v.Rows()/2; ++i)
      {
	DetType sumup = 0;
	DetType sumdown = 0;
	const uint len = up.Columns();
	for(uint k = 0; k < len; ++k)
	{
	  sumup += up(k, i) * v(0, k);//yes I know... do the transpose before...
	  sumdown += down(k, i) * v(0, len + k);
	}
	ret(i, 0) = sumup;
	ret(len + i, 0) = sumdown;
      }
      v = ret;
    }
    template<class Greensfunction, class A>
    DetType det(const typename Configuration::Container&, const A& alpha) const;
    template <class Greensfunction, class A>
    bool needsreinit(const typename Configuration::Container& cont, const A& alpha) const;
    template <class Greensfunction, class A>
    void reevaluate(const typename Configuration::Container&, const A& alpha, DetType& incphase);
    /**
    A function for doing the fast update of the determinants
    @param srcconfiguration from this we take the data for the new weights
    @param alpha the Alpha Object we use for determining the weights
    @param u2 will save the new values of the greensfunctions
    @param v1
    @return a pair that contains in .first the ratio of the up sector and in .second the ratio of the down sector
    */
    template <class GreensFunction, class A>
    RatioType detOfOneAddedVertex(const typename Configuration::value_type& vertex, const A& srcconfiguration, const Alpha<typename Configuration::FPType>& alpha, MatType& u2, MatType& v1) const;
    RatioType detOfOneRemovedVertex(uint pos) const MTL_PURE_FUNCTION;
    template <class A>
    RatioType det_of_a_flippedIsingSpin(uint pos, const A& srcconfiguration, const Alpha<typename Configuration::FPType>& alpha, RatioType& gstore) const;
    inline void matrixOfOneAddedVertex(const Matrix_Container<Configuration, true>& src, const RatioType& ratio, const MatType& u2, const MatType& v1);
    template <class GreensFunction, class A>
    inline void matrixOfOneRemovedVertex(const A& src, uint pos, const Alpha<FPType> alpha, const A& dst);
    inline void matrix_of_a_flipped_Ising_Spin(uint pos, const Matrix_Container<Configuration, true>& src, const RatioType& ratio, const RatioType& gstore);
    template <class Greensfunction>
    inline DetType measure(const typename Configuration::Container& cont, const typename Configuration::Container::value_type creator, const typename Configuration::Container::value_type destructor) const;
private:
    template<class Greensfunction, int Spin>
    inline DetType inSector(const typename Configuration::Container& cont, const typename Configuration::Container::value_type creator, const typename Configuration::Container::value_type destructor) const MTL_PURE_FUNCTION;
    /**
    A function for a fast update of the inverse Matrix via the sherman-Morrison- formula
    @param mat the Matrix to update
    @param src the source matrix
    @param detratio The ratio of the old and new determinant r = detold/detnew
    @param u2
    @param v1
    */
    inline void fastAdd(MatType& dst, const MatType& src, DetType detratio, const MatType& u2, const MatType& v1);
    /**
    A function to do the fast update of the inverse Matrix for a Configuration
    @param u2 the outer vector to the right that gets removed. It is assumed that it is already in the right permutation for the Sherman-Morrison formula
    @param v1 the lower vector on the bottom of the matrix. It is also assumed that it is in the right order
    @param mat the target Matrix that gets overwritten
    @param src the matrix we use as input.
    */
    static inline void fastRemove(const MatType& u2, MatType v1, MatType& dst, const MatType& src);
    template <class Greensfunction, class A>
    inline void createCorrectionVectors(const A& cont, const typename Configuration::Container::value_type& creator, const typename Configuration::Container::value_type& destructor, MatType& u2, MatType& v1) const;
};

template <class Configuration>
template <class Greensfunction>
typename Matrix_Container<Configuration, true>::DetType Matrix_Container<Configuration, true>::measure(const typename Configuration::Container& cont, const typename Configuration::Container::value_type creator, const typename Configuration::Container::value_type destructor) const
{
    DetType retval;
    if (creator.spin == destructor.spin)
    {
        if (creator.spin == UP)
            retval = inSector<Greensfunction, UP>(cont, creator, destructor);
        else
            retval = inSector<Greensfunction, DOWN>(cont, creator, destructor);
    }
    else
    {
       retval = 0;
//        throw std::runtime_error("Non-diagonal Greensfunction requested!!!");
    }
    return retval;
}

template <class Configuration>
template <class GreensFunction, int Spin>
typename Matrix_Container<Configuration, true>::DetType Matrix_Container<Configuration, true>::inSector(const typename Configuration::Container& cont, const typename Configuration::Container::value_type creator, const typename Configuration::Container::value_type destructor) const
{
    const unsigned int size = static_cast<unsigned int>(cont.size());
    DetType mnn = GreensFunction::eval(creator, destructor);
    if (size == 0) return mnn;
    MatType u2;
    MatType v1;
    createCorrectionVectors<GreensFunction, typename Configuration::Container>(cont, creator, destructor, u2, v1);//now generate the entries for the Observables
    //Calculate the resulting determinant
    return this->rank1Add_Det(mnn, (Spin == UP) ? up : down, v1, u2);
}

template <class Configuration>
template <class GreensFunction, class A>
void Matrix_Container<Configuration, true>::createCorrectionVectors(const A& cont, const typename Configuration::Container::value_type& creator, const typename Configuration::Container::value_type& destructor, MatType& u2, MatType& v1) const
{
    typename A::const_iterator it = cont.begin();
    unsigned int size = static_cast<unsigned int>(cont.size());
    u2 = MatType(size, 1u, MTLICCL::None<DetType>(0.0));//reallocate the matrices in the right size
    v1 = MatType(1u, size, MTLICCL::None<DetType>(0.0));//reallocate the matrices in the right size
    for (unsigned int i = 0; i < size; ++i, ++it)
    {
        u2(i, 0) = GreensFunction::eval(*it, destructor);
        v1(0, i) = GreensFunction::eval(creator, *it);
    }
    return;
}

template <class Configuration>
template <class A>
typename Matrix_Container<Configuration, true>::RatioType Matrix_Container<Configuration, true>::det_of_a_flippedIsingSpin(uint pos, const A& srcconfiguration, const Alpha<typename Configuration::FPType>& alpha, RatioType& gstore) const
{
  RatioType retval;
  SPINS spin = srcconfiguration[pos].spin;
  gstore.first = alpha(UP, spin) - alpha(UP, !spin);
  gstore.second = alpha(DOWN, spin) - alpha(DOWN, !spin);
  retval.first = this->det_of_diagonal_change(pos, gstore.first, up);
  retval.second = this->det_of_diagonal_change(pos, gstore.second, down);
  return retval;
}

template <class Configuration>
template <class GreensFunction, class A>
typename Matrix_Container<Configuration, true>::RatioType Matrix_Container<Configuration, true>::detOfOneAddedVertex(const typename Configuration::value_type& vertex, const A& srcconfiguration, const Alpha<typename Configuration::FPType>& alpha, MatType& u2, MatType& v1) const
{
    typedef typename GreensFunction::FreeGreensFunctionReturnValueType RetType;
    RetType mnn(GreensFunction::eval(vertex, vertex));
    RatioType retval;
    retval.first = mnn - alpha(UP, vertex.spin);
    retval.second = mnn - alpha(DOWN, vertex.spin);
    if (srcconfiguration.size() > 0)//if the src has vertices do the fast updates
    {
        createCorrectionVectors<GreensFunction>(srcconfiguration, vertex, vertex, u2, v1);
        retval.first = this->rank1Add_Det(retval.first, up, v1, u2);
        retval.second = this->rank1Add_Det(retval.second, down, v1, u2);
    }
    return retval;
}

template <class Configuration>
void Matrix_Container<Configuration, true>::fastRemove(const MatType& u2, MatType v1, MatType& dst, const MatType& src)
{
    const unsigned int size = src.Rows();
    const unsigned int sizem = size - 1;
    v1(0,sizem) = v1(0,sizem) - DetType(1.0);
//determine a lot of common vectors and scalars
    MatType alphavec(src * u2);
    DetType b((v1 * alphavec)(0,0));
    DetType g(1.0);
    for (unsigned int k = 0; k < sizem; ++k)
        g -= u2(k,0) * src(sizem, k);
    g = static_cast<FPType>(1.0)/g;
    DetType d(g * b);
    DetType f(static_cast<FPType>(1.0) - d * src(sizem, sizem));
    for (unsigned int k = 0; k < size; ++k)
        f -= v1(0,k) * src(k,sizem);
    f = static_cast<FPType>(1.0)/f;
    MatType beta(v1 * src);
    DetType c(f * g * src(sizem, sizem));
    DetType c2(g + c * d);
    dst = MatType(sizem, sizem, MTLICCL::None<DetType>(0.0));
    //We do the copy back to the dst matrix also in this step
    for (unsigned int k = 0; k < sizem; ++k)//it's sufficient to iterate up to sizem because we don't need the last row
    {
        DetType t0 = f * src(k,sizem);
        DetType t1 = c2 * alphavec(k,0) + d * t0;
        t0 += c * alphavec(k,0);
        for (unsigned int j = 0; j < sizem; ++j)//it's sufficient to iterate up to "sizem" because we don't need the last columns
            dst(k, j) = src(k,j) + dppd(t1, t0, src(sizem,j), beta(0,j));//t1 * src(sizem,j) + t2 * beta(0,j);
    }
    return;
}

template <class Configuration>
void Matrix_Container<Configuration, true>::matrix_of_a_flipped_Ising_Spin(uint pos, const Matrix_Container<Configuration, true>& src, const RatioType& ratio, const RatioType& gstore)
{
   this->matrix_of_diagonal_change(pos, ratio.first, gstore.first, up, src.up);
   this->matrix_of_diagonal_change(pos, ratio.second, gstore.second, down, src.down);
  return;
}

template <class Configuration>
template <class GreensFunction, class A>
void Matrix_Container<Configuration, true>::matrixOfOneRemovedVertex(const A& src, uint pos, const Alpha<FPType> alpha, const A& dst)
{
    //if the old configuration has 3 vertices then the new one has 2
    uint size = static_cast<uint>(src.size());
    if (size > 2)//if we have two or more vertices in the new configuration... do fast updates
    {
        const unsigned int sizem = size - 1;
        MatType u2(size, 1u, MTLICCL::None<DetType>(0.0));
        MatType v1(1u, size, MTLICCL::None<DetType>(0.0));
        typename Configuration::const_iterator rem = src.begin();
        //rem is the vertex to be removed
        advance(rem, pos);
        typename Configuration::const_iterator it = dst.begin();
        for (unsigned int k = 0; k < size; ++k, ++it)
        {
            u2(k, 0) = GreensFunction::eval(*it, *rem);
            v1(0, k) = GreensFunction::eval(*rem, *it);
        }
        //now that we know the positions of the elements we can apply their corrections
        u2(sizem, 0) = 0;
        v1(0, sizem) = v1(0, sizem) - DetType(alpha(UP, rem->spin));
        //swap the Matrix to the right form
        MatType srcmat(src.matcont.up);
        swapColumns(srcmat, sizem, pos);
        swapRows(srcmat, sizem, pos);
        fastRemove(u2, v1, up, srcmat);
        v1(0, sizem) = v1(0, sizem) + DetType(alpha(UP, rem->spin) - alpha(DOWN, rem->spin));
        srcmat = src.matcont.down;
        swapColumns(srcmat, sizem, pos);
        swapRows(srcmat, sizem, pos);
        fastRemove(u2, v1, down, srcmat);
    }
    else
    {
        if (size == 2)//just write out the inverses which are the inverse element
        {
            DetType mnn(GreensFunction::eval(dst[0], dst[0]) - DetType(alpha(UP, dst[0].spin)));
            up = MatType(1, 1, MTLICCL::None<DetType>(DetType(1.0)/mnn));
            mnn += DetType(alpha(UP, dst[0].spin) - alpha(DOWN, dst[0].spin));
            down = MatType(1, 1, MTLICCL::None<DetType>(DetType(1.0)/mnn));
        }
        else
        {
            up = MatType(0, 0, MTLICCL::None<DetType>(0));
            down = MatType(0, 0, MTLICCL::None<DetType>(0));
        }
    }
}

template <class Configuration>
typename Matrix_Container<Configuration, true>::RatioType Matrix_Container<Configuration, true>::detOfOneRemovedVertex(uint pos) const
{
    return std::make_pair(up(pos, pos), down(pos, pos));
}

template <class Configuration>
void Matrix_Container<Configuration, true>::fastAdd(MatType& dst, const MatType& src, DetType detratio, const MatType& u2, const MatType& v1)
{
    /*
        Assuming the same property as in the removal of vertices:
        (A^{-1})_{n,n} = det(A(N-1))/det(A(N))       (given by detratio)
        Thus we already know this single element of the inverse matrix.
        Inserting this we can simplify some expressions that result from the Sherman-Morrison-formula
        */
    detratio = 1.0/detratio;//to be consistent with the NxN case
    const unsigned int sizem = src.Rows();
    dst = MatType(src.Rows()+1, src.Columns()+1, MTLICCL::None<DetType>(0.0));
    dst(sizem, sizem) = detratio;//the entry in the lower right corner
    //fill up the outer right column of the matrix. It's the src matrix multiplied by the new outer right vector times (-detratio)
    for (unsigned int k = 0; k < sizem; ++k)
    {
        for (unsigned int j = 0; j < sizem; ++j)
            dst(k, sizem) += src(k,j) * u2(j,0);
        dst(k, sizem) *= (-detratio);
    }
//fill up the a vector. This is the one at the bottom
    MatType aR(v1 * src);
//fill up the new inverse Matrix
    for (unsigned int k = 0; k < sizem; ++k)
    {
        dst(sizem, k) = -detratio * aR(0, k);
        //do the rest that's necessary for the outer product
        for (unsigned int j = 0; j < sizem; ++j)
            dst(k, j) = src(k,j) - dst(k, sizem) * aR(0,j);
    }
    return;
}

template <class Configuration>
void Matrix_Container<Configuration, true>::matrixOfOneAddedVertex(const Matrix_Container<Configuration, true>& src, const RatioType& ratio, const MatType& u2, const MatType& v1)
{
    if (src.up.Rows() == 0)
    {
        MatType mat(1,1, MTLICCL::None<DetType>(0.0));
        mat(0,0) = DetType(1.0)/ratio.first;
        up = mat;
        mat(0,0) = DetType(1.0)/ratio.second;
        down = mat;
    }
    else
    {
        this->fastAdd(up, src.up, ratio.first, u2, v1);
        this->fastAdd(down, src.down, ratio.second, u2, v1);
    }
    return;
}

template <typename GreensFunctionType>
struct Reevaluate_Spin_Symmetric
    {};

template <>
struct Reevaluate_Spin_Symmetric<double>
{
    template<class GF, class Container, class MatCont>
    static void reevaluate(const Container& container, const Alpha<double>& alpha, MatCont& matcont, double& incphase)
    {
        typedef double FPType;
        const unsigned int size = static_cast<uint>(container.size());
        typedef typename MatCont::MatType MatType;
        //Construct the Matrix with Greensfunctions only
        MatType matrix(size, size, WickMatrixFiller<Container, GF>(container));
        //Fix the diagonal values
        //first determine the up sector
        typename Container::const_iterator it = container.begin();
        for (unsigned int i = 0; i < size; ++i, ++it)
            matrix(i, i) -= alpha(UP, it->spin);
        unsigned int idx[size];//to store the permutation vector of the LU Decomposition
        MatType lu(matrix);//temporal storage for the lu decomposed matrices
        incphase = ludecompose(lu, idx);//store the sign of the permutation
        for (unsigned int k = 0; k < lu.Rows(); ++k)
        {
            if (lu(k,k) < 0.0)
                incphase *= -1.0;
        }
        matcont.up = inversefromLU(lu, idx);
        //Now for the positive sigmas/down sector
        it = container.begin();
        for (unsigned int i = 0; i < size; ++i, ++it)
        {
            //add the correction for the positive sigmas again and then subtract the negative Spin-sector
            matrix(i,i) += (alpha(UP, it->spin) - alpha(DOWN, it->spin));
        }
        lu = matrix;
        incphase *= ludecompose(lu, idx);//store the sign of the permutation
        for (unsigned int k = 0; k < lu.Rows(); ++k)
        {
            if (lu(k,k) < 0.0)
                incphase *= -1.0;
        }
        matcont.down = inversefromLU(lu, idx);
    }
};

template <>
struct Reevaluate_Spin_Symmetric<float>
{
    template<class GF, class Container, class MatCont>
    static void reevaluate(const Container& container, const Alpha<float>& alpha, MatCont& matcont, float& incphase)
    {
        typedef float FPType;
        const unsigned int size = container.size();
        typedef typename MatCont::MatType MatType;
        //Construct the Matrix with Greensfunctions only
        MatType matrix(size, size, WickMatrixFiller<Container, GF>(container));
        //Fix the diagonal values
        //first determine the up sector
        typename Container::const_iterator it = container.begin();
        for (unsigned int i = 0; i < size; ++i, ++it)
            matrix(i, i) -= alpha(UP, it->spin);
        unsigned int idx[size];//to store the permutation vector of the LU Decomposition
        MatType lu(matrix);//temporal storage for the lu decomposed matrices
        incphase = ludecompose(lu, idx);//store the sign of the permutation
        for (unsigned int k = 0; k < lu.Rows(); ++k)
        {
            if (lu(k,k) < 0.0)
                incphase *= -1.0;
        }
        matcont.up = inversefromLU(lu, idx);
        //Now for the positive sigmas/down sector
        it = container.begin();
        for (unsigned int i = 0; i < size; ++i, ++it)
        {
            //add the correction for the positive sigmas again and then subtract the negative Spin-sector
            matrix(i,i) += (alpha(UP, it->spin) - alpha(DOWN, it->spin));
        }
        lu = matrix;
        incphase *= ludecompose(lu, idx);//store the sign of the permutation
        for (unsigned int k = 0; k < lu.Rows(); ++k)
        {
            if (lu(k,k) < 0.0)
                incphase *= -1.0;
        }
        matcont.down = inversefromLU(lu, idx);
    }
};

template <typename FPType>
struct Reevaluate_Spin_Symmetric<std::complex<FPType> >
{
    template<class GF, class Container, class MatCont>
    static void reevaluate(const Container& container, const Alpha<FPType>& alpha, MatCont& matcont, std::complex<FPType>& incphase)
    {
//std::cout<<"old phase: "<<this->phase<<std::endl;
//std::cout<<"old incphase: "<<this->incphase<<std::endl;
        typedef typename MatCont::MatType MatType;
        const unsigned int size = static_cast<unsigned int>(container.size());
        //Construct the Matrix with Greensfunctions only
        MatType matrix(size, size, WickMatrixFiller<Container, GF>(container));
        //Fix the diagonal values
        //first determine the up sector
        typename Container::const_iterator it = container.begin();
        for (unsigned int i = 0; i < size; ++i, ++it)
            matrix(i,i) -= alpha(UP, it->spin);
        unsigned int idx[size];//to store the permutation vector of the LU Decomposition
        MatType lu(matrix);//temporal storage for the lu decomposed matrices
        incphase = ludecompose(lu, idx);//store the sign of the permutation
        FPType tempphase = 0.0;
        for (unsigned int k = 0; k < lu.Rows(); ++k)
            tempphase += arg(lu(k, k));
        matcont.up = inversefromLU(lu, idx);
        //Now for the positive sigmas/down sector
        it = container.begin();
        for (unsigned int i = 0; i < size; ++i, ++it)
        {
            //add the correction for the positive sigmas again and then subtract the negative Spin-sector
            matrix(i,i) += (alpha(UP, it->spin) - alpha(DOWN, it->spin));
        }
        lu = matrix;
        incphase *= ludecompose(lu, idx);//store the sign of the permutation
        //let's determine in a numerical stable way the phase
        for (unsigned int k = 0; k < lu.Rows(); ++k)//we can just add the phases of the downsector to the phases of the up sector
            tempphase += std::arg(lu(k,k));
        matcont.down = inversefromLU(lu, idx);
        incphase *= std::polar(1.0, tempphase);
    }
};

template <class Configuration>
template <class Greensfunction, class A>
void Matrix_Container<Configuration, true>::reevaluate(const typename Configuration::Container& cont, const A& alpha, DetType& incphase)
{
    Reevaluate_Spin_Symmetric<DetType>::template reevaluate<Greensfunction, typename Configuration::Container, Matrix_Container<Configuration, true> >(cont, alpha, *this, incphase);
    return;
}

template <class Configuration>
template <class Greensfunction, class A>
bool Matrix_Container<Configuration, true>::needsreinit(const typename Configuration::Container& cont, const A& alpha) const
{
    typedef typename Configuration::FPType FPType;
    //determine the scalar product we are looking at.
    const unsigned int i = rand()%up.Rows();
    const unsigned int j = rand()%up.Rows();
    //its i'th line times j'th column
    DetType resup = 0;
    DetType resdown = 0;
    const typename Configuration::Container::value_type& v1 = *(cont.begin() + i);
    typename Configuration::Container::const_iterator v2i = cont.begin();
    for (unsigned int k = 0; k < static_cast<unsigned int>(up.Rows()); ++k, ++v2i)
    {
        typename Greensfunction::FreeGreensFunctionReturnValueType gf = Greensfunction::eval(v1, *v2i);
        resup += (gf - ( k == i ? alpha(UP, v2i->spin) : 0)) * up(k,j);
        resdown += (gf - ( k == i ? alpha(DOWN, v2i->spin) : 0)) * down(k,j);
    }
    if (i == j)
    {
        if ((std::abs(resup - static_cast<FPType>(1.0)) < 1000.0*tiny<FPType>()) && (std::abs(resdown - static_cast<FPType>(1.0)) < 1000.0*tiny<FPType>())) return false;//rebuild not necessary
    }
    else
    {
        if ((std::abs(resup) < 1000.0*tiny<FPType>()) && (std::abs(resdown) < 1000.0*tiny<FPType>())) return false;//rebuild not necessary
    }
    return true;//rebuild is necessary!
}

template <class Configuration>
template <class Greensfunction, class A>
typename Configuration::MatType::Elementtype Matrix_Container<Configuration, true>::det(const typename Configuration::Container& cont, const A& alpha) const
{
    DetType retval;
    const unsigned int size = cont.size();
    //Construct the Matrix with Greensfunctions only
    MatType matrix(size, size, WickMatrixFiller<typename Configuration::Container, Greensfunction>(cont));
    //Fix the diagonal values
    //first determine the up sector
    typename Configuration::Container::const_iterator it = cont.begin();
    for (unsigned int i = 0; i < size; ++i, ++it)
        matrix(i,i) -= alpha(UP, it->spin);
    unsigned int idx[size];//to store the permutation vector of the LU Decomposition
    MatType lu(matrix);//temporal storage for the lu decomposed matrices
    retval = ludecompose(lu, idx);//store the sign of the permutation
    for (unsigned int k = 0; k < lu.Rows(); ++k)
        retval *= lu(k,k);
    //Now for the positive sigmas/down sector
    it = cont.begin();
    for (unsigned int i = 0; i < size; ++i, ++it)
    {
        //add the correction for the positive sigmas again and then subtract the negative Spin-sector
        matrix(i,i) += (alpha(UP, it->spin) - alpha(DOWN, it->spin));
    }
    lu = matrix;
    retval *= ludecompose(lu, idx);//store the sign of the permutation
    //let's determine in a numerical stable way the phase
    for (unsigned int k = 0; k < lu.Rows(); ++k)//we can just add the phases of the downsector to the phases of the up sector
        retval *= lu(k,k);
    return retval;
}

template <class Configuration>
class Matrix_Container<Configuration, false> : private MatrixContainer_Common<typename Configuration::MatType>
{
public:
    typedef typename Configuration::MatType::Elementtype DetType;///< a typedef for the type of the elements and the determinant
    typedef typename Configuration::MatType MatType;///< the type of the Matrix that we use
    typedef typename Configuration::FPType FPType;
    typedef MTLICCL::Matrix<MTLICCL::Static<MTLICCL::Config<DetType, 2, 2> > > RatioType;
    MatType mat;
    void reset()
    {
        mat = MatType(0, 0, MTLICCL::None<DetType>(0.0));
    }
    inline DetType& operator()(uint i, uint k, SPINS s1, SPINS s2) const
    {
      if(s1 == s2)
      {
	if(s1 == UP)
	  return mat(2*i, 2*k);
	else
	  return mat(2*i + 1, 2*k + 1);
      }
      else
      {
	if(s1 == UP)
	  return mat(2*i + 1, 2*k);
	else
	  return mat(2*i, 2*k + 1);
      }
    }
    template<class VecType>
    void multiplyVectorbyConfiguration_left(VecType& v) const
    {
      v = mat * v;
    }
    template<class VecType>
    void multiplyVectorbyConfiguration_right(VecType& v) const
    {
      v =  v * mat;
    }
    Matrix_Container(int size) : mat(size, size, MTLICCL::None<DetType>(0.0)) {}
    Matrix_Container() {}
    template<class Greensfunction, class A>
    DetType det(const typename Configuration::Container&, const A& alpha) const;
    template <class Greensfunction, class A>
    bool needsreinit(const typename Configuration::Container& cont, const A& alpha) const;
    template <class Greensfunction, class A>
    void reevaluate(const typename Configuration::Container&, const A& alpha, DetType& incphase);
    /**
    This function will corespond to fastupdateweight used for Vertex addition
    @param u2 will save the new values of the greensfunctions
    @param v1
    */
    template <class GreensFunction, class A>
    RatioType detOfOneAddedVertex(const typename Configuration::value_type& vertex, const A& srcconfiguration, const Alpha<typename Configuration::FPType>& alpha, MatType& u2, MatType& v1) const;
    RatioType detOfOneRemovedVertex(uint pos) const;
    template <class A>
    RatioType det_of_a_flippedIsingSpin(uint pos, const A& srcconfiguration, const Alpha<typename Configuration::FPType>& alpha, RatioType& gstore) const;
    inline void matrixOfOneAddedVertex(const Matrix_Container<Configuration, false>& src, const RatioType& ratio, const MatType& u2, const MatType& v1);
    template <class GreensFunction, class A>
    inline void matrixOfOneRemovedVertex(const A& src, uint pos, const Alpha<FPType> alpha, const A& dst);
    inline void matrix_of_a_flipped_Ising_Spin(uint pos, const Matrix_Container<Configuration, false>& src, RatioType ratio, const RatioType& gstore);
    template <class Greensfunction>
    inline DetType measure(const typename Configuration::Container& cont, const typename Configuration::Container::value_type creator, const typename Configuration::Container::value_type destructor) const;
    template <class T> friend struct Reevaluate_Not_Spin_Symmetric;
private:
    template <class GreensFunction, typename T>
    static inline void create_All_four_Spin_Possiblilities(typename Configuration::Container::value_type v1, typename Configuration::Container::value_type v2, T (&retval) [4]);
    template <class Greensfunction, class A>
    inline void createCorrectionVectors(const A& cont, const typename Configuration::Container::value_type& creator, const typename Configuration::Container::value_type& destructor, MatType& u2, MatType& v1) const;
    static inline void rank2InverseUpdate(const uint blocksize, MatType& dst, const MatType& src, const MatType& u2, const MatType& v1, RatioType ratio);
    static inline void rank2FastRemove(const uint blocksize, MatType& dst, const MatType& src, const MatType& u2, const MatType& v1);
};

template <class Configuration>
template <class GreensFunction>
typename Matrix_Container<Configuration, false>::DetType Matrix_Container<Configuration, false>::measure(const typename Configuration::Container& cont, const typename Configuration::Container::value_type creator, const typename Configuration::Container::value_type destructor) const
{
    const uint blocksize = 2;
    const unsigned int size = static_cast<unsigned int>(cont.size());
    DetType mnn = GreensFunction::eval(creator, destructor);
    if (size == 0) return mnn;
    MatType u2(blocksize * size, 1u, MTLICCL::None<DetType>(0.0));
    MatType v1(1u, blocksize * size, MTLICCL::None<DetType>(0.0));//reallocate the matrices in the right size
    typename Configuration::const_iterator it = cont.begin();
    for (unsigned int i = 0; i < size; ++i, ++it)
    {
        typename Configuration::Container::value_type vert(*it);
        vert.spin = /*UP*/DOWN;//Not quite sure why, but it works...
        u2(blocksize * i, 0) = GreensFunction::eval(vert, destructor);
        v1(0, blocksize * i) = GreensFunction::eval(creator, vert);
        vert.spin = /*DOWN*/UP;
        u2(blocksize * i + 1, 0) = GreensFunction::eval(vert, destructor);
        v1(0, blocksize * i + 1) = GreensFunction::eval(creator, vert);
    }
    return this->rank1Add_Det(mnn, mat, v1, u2);
}

template <typename MatType>
void copySquareBlock(MatType& dst, const MatType& src, uint blocksize, uint src_lpos, uint src_cpos)
{
    for (uint i = 0; i < blocksize; ++i)
        for (uint k = 0; k < blocksize; ++k)
            dst(i, k) = src(src_lpos + i, src_cpos + k);
    return;
}

template <class Configuration>
inline void Matrix_Container<Configuration, false>::rank2FastRemove(const uint blocksize, MatType& dst, const MatType& src, const MatType& u2, const MatType& v1)
{//about 40% of time is spent here
    MatType Z(v1*src);//f in the notes
    const MatType e(src * u2);
    MatType Q(blocksize, blocksize, MTLICCL::Identity<DetType>(1.0));//S in the notes
    MatType X(Q);
    for (uint i = 0; i < blocksize; ++i)
        for (uint k = 0; k < blocksize; ++k)
            Q(i,k) -= Z(i, src.Rows() - blocksize + k);
    Q = inverse(Q);//Q now contains S^-1
    Z = Q * Z;//this is R in the notes we calculate the last element of this vector, too
    const MatType p(Z*u2);// the last entry of u2 is zero, so this could be saved
    Q = Q * p;//this is Q in the notes
    MatType A_dd(blocksize, blocksize, MTLICCL::None<DetType>(0.0));
    copySquareBlock(A_dd, src, blocksize, src.Rows() - blocksize, src.Columns() - blocksize);
    X = X - A_dd*Q;
    for (uint i = 0; i < blocksize; ++i)
    {
        for (uint k = 0; k < blocksize; ++k)
        {
            for (uint l = 0; l  < (u2.Rows() - 1); ++l)//remember u2's last element is zero
                X(i,k) -= src(src.Rows() - blocksize + i, l ) * u2(l, k);
        }
    }
    X = inverse(X);
    for (uint i = 0; i < dst.Rows(); ++i)
    {
      const uint maxk = (dst.Columns()/2)*2;
        for (uint k = 0; k < maxk; k += 2)
        {
            typename MatType::value_type reg0 = src(i,k);
	    typename MatType::value_type reg1 = src(i,k + 1);
            for (uint l = 0; l < blocksize; ++l)
	    {
                reg0 += src(i, src.Columns() - blocksize + l)*Z(l, k);
		reg1 += src(i, src.Columns() - blocksize + l)*Z(l, k+1);
	    }
            for (uint l = 0; l < blocksize; ++l)
            {
                typename MatType::value_type rega = e(i,l);
                typename MatType::value_type regb0 = 0;
		typename MatType::value_type regb1 = 0;
                for (uint p = 0; p < blocksize; ++p)
                    rega += src(i, src.Columns() - blocksize + p) * Q(p, l);
//			        asm volatile("nop;nop;nop;");
                for (uint q = 0; q < blocksize; ++q)
                {
                    typename MatType::value_type regc0 = src(src.Rows() - blocksize + q,k);
		    typename MatType::value_type regc1 = src(src.Rows() - blocksize + q,k + 1);
                    for (uint s = 0; s < blocksize; ++s)
		    {
                        regc0 += A_dd(q,s) * Z(s,k);
			regc1 += A_dd(q,s) * Z(s,k+1);
		    }
                    regb0 += X(l,q) * regc0;
		    regb1 += X(l,q) * regc1;
                }
                reg0 += rega * regb0;
		reg1 += rega * regb1;
            }
            dst(i, k) = reg0;
	    dst(i, k + 1) = reg1;
        }
        
        for (uint k = maxk; k < dst.Columns(); ++k)
        {
            typename MatType::value_type reg = src(i,k);
            for (uint l = 0; l < blocksize; ++l)
                reg += src(i, src.Columns() - blocksize + l)*Z(l, k);
            for (uint l = 0; l < blocksize; ++l)
            {
                typename MatType::value_type rega = e(i,l);
                typename MatType::value_type regb = 0;
                for (uint p = 0; p < blocksize; ++p)
                    rega += src(i, src.Columns() - blocksize + p) * Q(p, l);
//			        asm volatile("nop;nop;nop;");//time is wasted around here, mostly because gcc doesn't propagate blocksize
                for (uint q = 0; q < blocksize; ++q)
                {
                    typename MatType::value_type regc = src(src.Rows() - blocksize + q,k);
                    for (uint s = 0; s < blocksize; ++s)
                        regc += A_dd(q,s) * Z(s,k);
                    regb += X(l,q) * regc;
                }
                reg += rega * regb;
            }
            dst(i, k) = reg;
        }
        
    }
    return;
}

template <class Configuration>
template <class GreensFunction, class A>
inline void Matrix_Container<Configuration, false>::matrixOfOneRemovedVertex(const A& src, uint pos, const Alpha<FPType> alpha, const A& dst)
{
    constexpr uint deltan = 2;
    //if the old configuration has 3 vertices then the new one has 2
    uint size = static_cast<uint>(src.size());
    if (size > 2)//if we have two or more vertices in the new configuration... do fast updates
    {
        mat = MatType(deltan*(size - 1), deltan*(size - 1) , MTLICCL::None<DetType>(0.0));//reallocate matrix in the right size
        MatType u2 = MatType(2*size, deltan, MTLICCL::None<DetType>(0.0));//reallocate the matrices in the right size
        MatType v1 = MatType(deltan, 2*size, MTLICCL::None<DetType>(0.0));//reallocate the matrices in the right size
        typename Configuration::const_iterator rem = src.begin();
        //rem is the vertex to be removed it's at position pos in the old configuration
        advance(rem, pos);
        createCorrectionVectors<GreensFunction>(dst, *rem, *rem, u2, v1);//create correction vectors. does not write to the last position, so we have to update v1
        DetType *const block[4] = {&v1(0, 2*size - deltan), &v1(0, 2*size - deltan + 1), &v1(1, 2*size - deltan), &v1(1,2*size - deltan +1)};
        create_All_four_Spin_Possiblilities<GreensFunction>(*rem, *rem, block);
        v1(0, 2*size - deltan) -= alpha(UP, rem->spin) - 1.0;
        v1(1, 2*size - deltan + 1) -= alpha(DOWN, rem->spin) - 1.0;
        //let's swap the matrix into the right form. we move the to-be-removed-lines to the end
        MatType srcmat(src.matcont.mat);
        for (uint i = 0; i < deltan; ++i)
            swapRows(srcmat, deltan*(size - 1) + i, deltan * pos + i);//FIXME: Works for Matrices of square blocks
        for (uint i = 0; i < deltan; ++i)
            swapColumns(srcmat, deltan*(size - 1) + i, deltan * pos + i);//FIXME: Works for Matrices of square blocks
        rank2FastRemove(deltan, mat, srcmat, u2, v1);
    }
    else
    {
        if (size == 2)//just write out the inverses which are the inverse element
        {
            mat = MatType(deltan*1, deltan*1 , MTLICCL::None<DetType>(0.0));//reallocate matrix in the right size
            DetType *const block[4] = {&mat(0,0), &mat(0,1), &mat(1,0), &mat(1,1)};
            create_All_four_Spin_Possiblilities<GreensFunction>(dst[0], dst[0], block);
            mat(0,0) -= alpha(UP, dst[0].spin);
            mat(1,1) -= alpha(DOWN, dst[0].spin);
            mat = inverse(mat);
        }
        else
        {
            mat = MatType(0, 0, MTLICCL::None<DetType>(0));
        }
    }
    return;
}

template <class Configuration>
template <class A>
typename Matrix_Container<Configuration, false>::RatioType Matrix_Container<Configuration, false>::det_of_a_flippedIsingSpin(uint pos, const A& srcconfiguration, const Alpha<typename Configuration::FPType>& alpha, RatioType& gstore) const
{
  /*
  The general formula is:
  det(A + e_i * v^T) = det(A) *det(1 + a * E)
  a is the correction matrix to the old entry E.
  the case of an empty matrix should be catched beforehand
  */
  const uint deltan = 2;
  RatioType ret;
  SPINS spin = srcconfiguration[pos].spin;
  for (uint i = 0; i < deltan; ++i)
  {
      SPINS truespin = (i == 0? UP : DOWN);//FIX THIS FOR THE GENERAL CASE!
      gstore(i,i) = alpha(truespin, spin) - alpha(truespin, !spin);
      for (uint k = 0; k < deltan; ++k)
      {
          ret(i, k) = (i == k? 1.0 : 0.0) +  gstore(i,i) * srcconfiguration.matcont.mat(deltan * pos + i, deltan * pos + k);
      }
  }
  return ret;
}

template <class Configuration>
typename Matrix_Container<Configuration, false>::RatioType Matrix_Container<Configuration, false>::detOfOneRemovedVertex(uint pos) const
{
    const uint deltan = 2;
    RatioType ret;
    for (uint i = 0; i < deltan; ++i)
    {
        for (uint k = 0; k < deltan; ++k)
            ret(i, k) = mat(deltan*pos + i, deltan*pos + k);
    }
    return ret;
}

template <class Configuration>
void Matrix_Container<Configuration, false>::matrix_of_a_flipped_Ising_Spin(uint pos, const Matrix_Container<Configuration, false>& src, RatioType ratio, const RatioType& gstore)
{
  //the target is mat
  const uint n = 2;
  uint size = src.mat.Rows();
  ratio = inverse(ratio);
  ratio = ratio * gstore;
  MatType avec(size, n, MTLICCL::InitFromMatrix<MatType>(src.mat,0 ,pos * n));
  MatType bvec(n, size, MTLICCL::InitFromMatrix<MatType>(src.mat, pos *n, 0));
  mat = src.mat - avec * ratio * bvec;
  return;
}

template <class Configuration>
void Matrix_Container<Configuration, false>::rank2InverseUpdate(const uint blocksize, MatType& dst, const MatType& src, const MatType& u2, const MatType& v1, RatioType ratio)
{
    const uint oldsiz = src.Rows();
    const uint newsiz = oldsiz + blocksize;
    dst = MatType(newsiz, newsiz, MTLICCL::None<DetType>(0.0));
    ratio = inverse(ratio);
    for (uint i = 0; i < blocksize; ++i)
        for (uint k = 0; k < blocksize; ++k)
            dst(oldsiz + i, oldsiz + k) = ratio(i, k);
    MatType alpha(v1 * src);
    alpha = ratio * alpha;
    for (uint i = 0; i < blocksize; ++i)
        for (uint k = 0; k < oldsiz; ++k)
            dst(oldsiz + i, k) = - alpha(i, k);//the minus sign is necessary for consistency with the 1D case
    MatType beta(src * u2);
    //let's create the main part of the matrix
    for (uint i = 0; i < oldsiz; ++i)
        for (uint k = 0; k < oldsiz; ++k)
        {
            dst(i,k) = src(i,k);
            for (uint l = 0; l < blocksize; ++l)
                dst(i,k) += beta(i,l) * alpha(l,k);
        }
    beta = beta * ratio;
    for (uint i = 0; i < oldsiz; ++i)
        for (uint k = 0; k < blocksize; ++k)
            dst(i, oldsiz + k) = - beta(i, k);
    return;
}

template <class Configuration>
inline void Matrix_Container<Configuration, false>::matrixOfOneAddedVertex(const Matrix_Container<Configuration, false>& src, const RatioType& ratio, const MatType& u2, const MatType& v1)
{
    if (src.mat.Rows() == 0)
    {
        mat = inverse(ratio);
    }
    else
        rank2InverseUpdate(2, mat, src.mat, u2, v1, ratio);
    return;
}

template <class T>
struct Access
{
    static inline T& get(T& t) {
        return t;
    }
};

template <class T>
struct Access<T *const>
{
    static inline T& get(T *const  t) {
        return *t;
    }
};

template <class Configuration>
template <class GreensFunction, typename T>
inline void Matrix_Container<Configuration, false>::create_All_four_Spin_Possiblilities(typename Configuration::Container::value_type v1, typename Configuration::Container::value_type v2, T (&retval) [4])
{
    v1.spin = UP;
    v2.spin = UP;
    Access<T>::get(retval[0]) = GreensFunction::eval(v1, v2);

    v1.spin = DOWN;
    v2.spin = UP;
    Access<T>::get(retval[1]) = GreensFunction::eval(v1, v2);

    v1.spin = UP;
    v2.spin = DOWN;
    Access<T>::get(retval[2]) = GreensFunction::eval(v1, v2);

    v1.spin = DOWN;
    v2.spin = DOWN;
    Access<T>::get(retval[3]) = GreensFunction::eval(v1, v2);
    return;
}

template <class Configuration>
template <class Greensfunction, class A>
inline void Matrix_Container<Configuration, false>::createCorrectionVectors(const A& cont, const typename Configuration::Container::value_type& creator, const typename Configuration::Container::value_type& destructor, MatType& u2, MatType& v1) const
{
    typename A::const_iterator it = cont.begin();
    unsigned int size = static_cast<unsigned int>(cont.size());
    for (unsigned int i = 0; i < size; ++i, ++it)
    {
        DetType *const block[4] = {&u2(2*i,0), &u2(2*i,1), &u2(2*i+1,0), &u2(2*i+1,1)};
        create_All_four_Spin_Possiblilities<Greensfunction>(*it, destructor, block);

        DetType *const block2[4] = {&v1(0,2*i), &v1(0,2*i+1), &v1(1,2*i), &v1(1,2*i+1)};
        create_All_four_Spin_Possiblilities<Greensfunction>(creator, *it, block2);
    }
    return;
}

template <class Configuration>
template <class GreensFunction, class A>
typename Matrix_Container<Configuration, false>::RatioType Matrix_Container<Configuration, false>::detOfOneAddedVertex(const typename Configuration::value_type& vertex, const A& srcconfiguration, const Alpha<typename Configuration::FPType>& alpha, MatType& u2, MatType& v1) const
{
    typedef typename GreensFunction::FreeGreensFunctionReturnValueType RetType;
    const uint size = static_cast<uint>(srcconfiguration.size());
    RatioType retval;
    DetType *const block[4] = {&retval(0,0), &retval(0,1), &retval(1,0), &retval(1,1)};
    create_All_four_Spin_Possiblilities<GreensFunction>(vertex, vertex, block);
    retval(0,0) -= alpha(UP, vertex.spin);
    retval(1,1) -= alpha(DOWN, vertex.spin);
    if (size > 0)//if the src has vertices do the fast updates
    {
        u2 = MatType(2*size, 2u, MTLICCL::None<DetType>(0.0));//reallocate the matrices in the right size
        v1 = MatType(2u, 2*size, MTLICCL::None<DetType>(0.0));//reallocate the matrices in the right size
        createCorrectionVectors<GreensFunction>(srcconfiguration, vertex, vertex, u2, v1);
        retval = retval - MTLICCL::Matrix<MTLICCL::Static<MTLICCL::Config<RetType, 2, 2> > >(v1 * (mat * u2));
    }
    return retval;
}

template <typename GreensFunctionType>
struct Reevaluate_Not_Spin_Symmetric
    {};

template <>
struct Reevaluate_Not_Spin_Symmetric<double>
{
    template<class GreensFunction, class Container, class MatCont>
    static void reevaluate(const Container& cont, const Alpha<double>& alpha, MatCont& matcont, double& incphase)
    {
        const uint blocksize = 2u;
        const unsigned int size = static_cast<uint>(cont.size());
        //Construct the Matrix with Greensfunctions only
        typename MatCont::MatType matrix(blocksize*size, blocksize*size, MTLICCL::None<double>(0.0));
        for (unsigned int k = 0; k < size; k++)
            for (unsigned int i = 0; i < size; i++)
            {
                double* const block[4] = {&matrix(2*i, 2*k), &matrix(2*i, 2*k + 1), &matrix(2*i + 1, 2*k), &matrix(2*i + 1, 2*k + 1)};
                MatCont::template create_All_four_Spin_Possiblilities<GreensFunction>(*(cont.begin() + i), *(cont.begin() + k), block);
            }
        for (unsigned int k = 0; k < size; k++)
        {
            matrix(blocksize*k, blocksize*k) -= alpha(UP, (cont.begin() + k)->spin);
            matrix(blocksize*k + 1, blocksize*k + 1) -= alpha(DOWN, (cont.begin() + k)->spin);
        }
        unsigned int idx[matrix.Rows()];//to store the permutation vector of the LU Decomposition
        incphase = ludecompose(matrix, idx);//store the sign of the permutation
        for (unsigned int k = 0; k < matrix.Rows(); ++k)
            if (matrix(k,k) < 0.0)
                incphase *= -1.0;
        matcont.mat = inversefromLU(matrix, idx);
    }
};

template <typename FPType>
struct Reevaluate_Not_Spin_Symmetric<std::complex<FPType> >
{
    template<class GreensFunction, class Container, class MatCont>
    static void reevaluate(const Container& cont, const Alpha<FPType>& alpha, MatCont& matcont, std::complex<FPType>& incphase)
    {
        const unsigned int size = static_cast<unsigned int>(cont.size());
        //Construct the Matrix with Greensfunctions only
        typename MatCont::MatType matrix(2*size, 2*size, MTLICCL::None<std::complex<FPType> >(0.0));
        for (unsigned int k = 0; k < size; k++)
            for (unsigned int i = 0; i < size; i++)
            {
                std::complex<FPType>* const block[4] = {&matrix(2*i, 2*k), &matrix(2*i, 2*k + 1), &matrix(2*i + 1, 2*k), &matrix(2*i + 1, 2*k + 1)};
                MatCont::template create_All_four_Spin_Possiblilities<GreensFunction>(*(cont.begin() + i), *(cont.begin() + k), block);
            }
        for (unsigned int k = 0; k < size; k++)
        {
            matrix(2*k, 2*k) -= alpha(UP, (cont.begin() + k)->spin);
            matrix(2*k + 1, 2*k + 1) -= alpha(DOWN, (cont.begin() + k)->spin);
        }
        unsigned int idx[matrix.Rows()];//to store the permutation vector of the LU Decomposition
        incphase = ludecompose(matrix, idx);//store the sign of the permutation
        FPType tempphase = 0.0;
        for (unsigned int k = 0; k < matrix.Rows(); ++k)
            tempphase += std::arg(matrix(k,k));
        incphase *= std::polar(1.0, tempphase);
        matcont.mat = inversefromLU(matrix, idx);
    }
};

template <class Configuration>
template <class Greensfunction, class A>
void Matrix_Container<Configuration, false>::reevaluate(const typename Configuration::Container& cont, const A& alpha, DetType& incphase)
{
    Reevaluate_Not_Spin_Symmetric<DetType>::template reevaluate<Greensfunction, typename Configuration::Container, Matrix_Container<Configuration, false> >(cont, alpha, *this, incphase);
}

template <class Configuration>
template <class GreensFunction, class A>
typename Configuration::MatType::Elementtype Matrix_Container<Configuration, false>::det(const typename Configuration::Container& cont, const A& alpha) const
{//Note the creative mixing of Ising spin and electron spin
    const unsigned int size = cont.size();
    //Construct the Matrix with Greensfunctions only
    MatType matrix(2*size, 2*size,  MTLICCL::None<DetType>(0.0));
    for (unsigned int k = 0; k < size; k++)
        for (unsigned int i = 0; i < size; i++)
        {
            DetType* const block[4] = {&matrix(2*i, 2*k), &matrix(2*i, 2*k + 1), &matrix(2*i + 1, 2*k), &matrix(2*i + 1, 2*k + 1)};
            create_All_four_Spin_Possiblilities<GreensFunction>(*(cont.begin() + i), *(cont.begin() + k), block);
        }
    for (unsigned int k = 0; k < size; k++)
    {
        matrix(2*k, 2*k) -= alpha(UP, (cont.begin() + k)->spin);
        matrix(2*k + 1, 2*k + 1) -= alpha(DOWN, (cont.begin() + k)->spin);
    }
    unsigned int idx[matrix.Rows()];//to store the permutation vector of the LU Decomposition
    DetType retval = ludecompose(matrix, idx);//store the sign of the permutation
    for (unsigned int k = 0; k < matrix.Rows(); ++k)
        retval *= matrix(k,k);
    return retval;
}

template <class Configuration>
template <class GreensFunction, class A>
bool Matrix_Container<Configuration, false>::needsreinit(const typename Configuration::Container& cont, const A& alpha) const
{
    //determine the scalar product we are looking at.
    const unsigned int i = rand()%static_cast<unsigned int>(cont.size());
    const unsigned int j = rand()%static_cast<unsigned int>(cont.size());//Not sure whether i==j can be generated
    const typename Configuration::Container::value_type& v1 = *(cont.begin() + i);
    typename Configuration::Container::const_iterator v2i = cont.begin();
    //its i'th line times j'th column
    DetType mat22[4] = {0};
    for (unsigned int k = 0; k < cont.size(); k++, ++v2i)
    {
        typename Configuration::Container::value_type v2 = *v2i;
        DetType block[4];
        create_All_four_Spin_Possiblilities<GreensFunction>(v1, v2, block);
        if (k == i)
        {
            SPINS isingspin = v1.spin;
            block[0] -= alpha(UP, isingspin);
            block[3] -= alpha(DOWN, isingspin);
        }
        mat22[0] += block[0] * mat(2*k, 2*j) + block[1] * mat(2*k + 1, 2*j);
        mat22[1] += block[0] * mat(2*k, 2*j + 1) + block[1] * mat(2*k + 1, 2*j + 1);
        mat22[2] += block[2] * mat(2*k,2*j) + block[3] * mat(2*k+1, 2*j);
        mat22[3] += block[2] * mat(2*k, 2*j + 1) + block[3] * mat(2*k+1, 2*j + 1);
    }
//    std::cout<<i<<" "<<j<<std::endl;
//    std::cout<<mat22[0]<<" "<<mat22[1]<<" "<<mat22[2]<<" "<<mat22[3]<<std::endl;
    DetType expected[4] = {0};
    if (i == j)
    {
        expected[0] = DetType(1.0);
        expected[3] = DetType(1.0);
    }
    bool norebuildnecessary = true;
    for (uint k = 0; k < 4; ++k)
    {
        norebuildnecessary = norebuildnecessary && (abs(mat22[k] - expected[k]) < 1000.0*tiny<FPType>());
    }
    return !norebuildnecessary;
}

#endif
