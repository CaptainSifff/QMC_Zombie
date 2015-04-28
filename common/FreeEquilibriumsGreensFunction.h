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
#ifndef FREEEQUILIBRIUMSGREENFUNCTION_H
#define FREEEQUILIBRIUMSGREENFUNCTION_H
#include <limits>
#include "ddqmc.h"
#include "Vertex.h"
#include "Greensfunction.h"

#ifdef __SSE3__
#include "emmintrin.h"
#include "pmmintrin.h"
#if defined( __SSE4_1__ )
#include "smmintrin.h"
#endif
#endif

/**
An implementation of a free-particle greens-function.
I assume that the lattice constant is a = 1
*/
template<typename FPType = float>
class FreeHubbardGreensFunction
{
public:
    enum {timeevolution = 0,
    has_Spin_symmetry = true,
    has_Giomega = false,
    Is_Impurity_model = false,
    has_TRS = true,
    };
    typedef Hubbard_Vertex<FPType> Vertex;
    typedef FPType FreeGreensFunctionReturnValueType;///< a typedef of the data type of the free greensfunction
    /**
    This evaluates the value of the free particle Green-function for the two given vertices
    @param v1 the first vertex
    @param v2 the second vertex
    @return the value of the free greensfunction evaluated with the given vertices
    */
    static inline FreeGreensFunctionReturnValueType eval(const Vertex& v1, const Vertex& v2) throw();
    /**
    A function for initializing the Tables that Make up the Greensfunction. Beta and the nr of sites are read from Config Files.
    @param CFG a class that contains the necessary information to extract the parameters.
    */
    template <class CFG>
    static inline void init(CFG&);
    /**
    This function frees the memory used by the tables for the Greensfunction
    */
    static inline void tidyup();
    /**
    To access the nr of Atoms of the chain
    @return the nr of atoms in the chain
    */
    static unsigned int getLen() throw()
    {
        return N;
    }
    /**
    Get the Length of the interactioninterval.
    @return the inverse Temperature beta
    */
    static FPType getContourLen() throw()
    {
        return beta;
    }
private:
    static FreeGreensFunctionReturnValueType** g;///< the free particle greensfunction
    static unsigned int slices;///< the number of timeslices. this means the resolution of the free greensfunction in tau space
    /**
    the length of one timeslice on the tau axis
    */
    static FPType betaslice;
    /**
    the inverse temperature beta
    */
    static FPType beta;
    /**
    the length of the chain
    */
    static unsigned int N;
    /**
    the Dispersionrelation -2t cos(k) -my . Note that the chemical potential is included here.
    @param k the momentum
    @param t the kinetic energy
    @param my the chemical potential
    @return the energy
    */
    static inline FPType disp(FPType k, FPType t, FPType my) throw();
    /**
    The free particle Greens function in k-space
    */
    static inline FPType GreenNullKSpace(FPType, FPType, FPType, FPType) throw();
};

template<typename FPType>
typename FreeHubbardGreensFunction<FPType>::FreeGreensFunctionReturnValueType** FreeHubbardGreensFunction<FPType>::g = NULL;

template<typename FPType>
FPType FreeHubbardGreensFunction<FPType>::betaslice;

template<typename FPType>
FPType FreeHubbardGreensFunction<FPType>::beta;

template<typename FPType>
unsigned int FreeHubbardGreensFunction<FPType>::N;

template<typename FPType>
unsigned int FreeHubbardGreensFunction<FPType>::slices;

template <typename FPType>
FPType FreeHubbardGreensFunction<FPType>::disp(FPType k, FPType t, FPType my) throw()//some system function seems to be also called epsilon
{
    return -2.0 * t * std::cos(k) - my;
}

template <typename FPType>
void FreeHubbardGreensFunction<FPType>::tidyup()
{
    for (unsigned int j = 0; j < (slices + 1); ++j)
        delete [] g[j];
    delete [] g;
}

template <typename FPType>
FPType FreeHubbardGreensFunction<FPType>::GreenNullKSpace(FPType k, FPType tau, FPType t, FPType my) throw()
{
    return std::exp( disp(k, t, my) * tau) * fermi(beta * (disp(k, t, my) ));
}

template <typename FPType>
template <class CFG>
void FreeHubbardGreensFunction<FPType>::init(CFG& curparams)
{
    //lattice sites
    N = curparams.N;
    beta = curparams.beta;
    slices = 100000;//Number of TimeSlices
    const unsigned int slicesp = slices + 1;
    betaslice = beta / static_cast<FPType>(slices);
    g = new FreeGreensFunctionReturnValueType*[slicesp];
    if (N == 1)//take special care of the 1 site hubbard model
    {
        for (unsigned int j = 0; j < slicesp; ++j)
        {
            g[j] = new FreeGreensFunctionReturnValueType[1];
            g[j][0] = GreenNullKSpace(0, j * betaslice, curparams.t, curparams.mu);
//            std::cout<<j * betaslice<<" "<<g[j][0]<<std::endl;
        }
        return;
    }
#pragma omp parallel for
    for (int i = 0; i < static_cast<int>(slicesp); ++i)//for every timeslice
    {
        g[i] = new FreeGreensFunctionReturnValueType[N];
        for (int j = 0; j < static_cast<int>(N); ++j)//for every realspacepoint
        {
            std::complex<FPType> tempsum = 0;
            for (FPType k = 0; k < static_cast<FPType>(N); ++k)//sum over all k-space values
                tempsum += exp(std::complex<FPType>(0.0, 1.0 * 2.0 * M_PI/N * k * j) ) * GreenNullKSpace(k * 2.0 * M_PI /  N, i * betaslice, curparams.t, curparams.mu);//R_j = j, because the lattice constant is assumed to be equal to one.
            if (imag(tempsum)> 0.05)//some sanity check, that the imaginary part is really small
            {
                std::cout<<"Error imaginary Part of the free Greenfunction seems not to vanish!"<<std::endl;
                exit(-1);
            }
            g[i][j] = real(tempsum)/N;// the Fourier transform has 1/N as normalization
        }
    }
    //Debugging Output
    /*    for(unsigned int i = 0; i < slicesp; ++i)
        {
    std::cout<<(i * betaslice)<<" ";
          for(int j = 0; j < static_cast<int>(N); ++j)//for every realspacepoint
            std::cout<<g[i][j]<<" ";
          std::cout<<std::endl;
        }*/
    return;
}

template <typename T>
std::complex<T> interpolate(std::complex<T>& f0, std::complex<T>& fx, std::complex<T>& fy, T& delx, T& dely, T dx, T dy) MTL_PURE_FUNCTION;

template <typename T>
std::complex<T> interpolate(std::complex<T>& f0, std::complex<T>& fx, std::complex<T>& fy, T& delx, T& dely, T dx, T dy)
{
    std::complex<T> ablx = (fx - f0)/delx;
    std::complex<T> ably = (fy - f0)/dely;
    return f0 + dx*ablx + dy * ably;
}

#ifdef __SSE3__

template <>
std::complex<double> interpolate<double>(std::complex<double>& f0, std::complex<double>& fx, std::complex<double>& fy, double& delx, double& dely, double dx, double dy) MTL_PURE_FUNCTION;

template <>
std::complex<double> interpolate<double>(std::complex<double>& f0, std::complex<double>& fx, std::complex<double>& fy, double& delx, double& dely, double dx, double dy)
{
    struct Local
    {
        union
        {
            double data[2];
            __v2df vec __attribute__ ((aligned (16)));
        };
    } f0v;
    f0v.vec = _mm_set_pd(f0.imag(), f0.real());

    __v2df fxv = _mm_set_pd(fx.imag(), fx.real());
    __v2df fyv = _mm_set_pd(fy.imag(), fy.real());
    __v2df delxv = _mm_set1_pd(delx);
    __v2df delyv = _mm_set1_pd(dely);
    __v2df dv = _mm_set_pd(dy, dx);

    fxv = _mm_sub_pd(fxv, f0v.vec);
    fyv = _mm_sub_pd(fyv, f0v.vec);
    fxv = _mm_div_pd(fxv, delxv);
    fyv = _mm_div_pd(fyv, delyv);
    delxv = fxv;
    fxv = _mm_shuffle_pd(fxv, fyv, 0);
    delxv = _mm_shuffle_pd(delxv, fyv, 3);
#ifdef __SSE4_1__
    fxv = _mm_dp_pd(dv, fxv, 49);//result stored in high part
    delxv = _mm_dp_pd(dv, delxv, 50);//result stored in low part
    delxv = _mm_move_sd(delxv, fxv);//combine high and low part of fxv, fyv in fxv
    f0v.vec = _mm_add_pd(f0v.vec, delxv);
#elif defined(__SSE3__)
    delxv = _mm_mul_pd(dv, delxv);
    fxv = _mm_mul_pd(dv, fxv);
    fxv = _mm_hadd_pd(fxv, delxv);
    f0v.vec = _mm_add_pd(f0v.vec, fxv);
#endif
    return std::complex<double>(f0v.data[0], f0v.data[1]);
}
#endif

template<typename FPType>
typename FreeHubbardGreensFunction<FPType>::FreeGreensFunctionReturnValueType FreeHubbardGreensFunction<FPType>::eval(const Vertex& v1, const Vertex& v2) throw()
{
    //determine the Differences between the two
    FPType delta_tau = v1.tau - v2.tau;
    int delta = v1.site - v2.site;
    //Take care of negative values
    if (delta < 0) delta = N + delta;//because delta is in this case already a negative number
    if (delta_tau == 0)
    {
        //return only the particle number
        return g[0][delta];
    }
    if (delta_tau == beta)
        return g[slices][delta];
    FPType sign = 1.0;
    if (delta_tau < 0)//Again take care of negative values
    {
        delta_tau = beta + delta_tau;
        sign = -1.0;//The Matsubara GreensFunction changes sign if we add a beta!!!!!!!!!!
    }
    FPType fptau_idx0;
    FPType rem = std::modf(delta_tau/betaslice, &fptau_idx0);//round to the smaller index and determine how far we're of
    long int tau_idx0 = lround(fptau_idx0);
    return lerp(rem, g[tau_idx0][delta], g[tau_idx0 + 1][delta]) * sign;//return the value of the greensfunction
}

/**
An implementation of a free-particle greens-function for real time evolution
I assume that the lattice constant is a = 1
*/
template<typename FPType = double>
class FreeRealTimeHubbardGreensFunction
{
public:
    enum {timeevolution = 1,
    has_Spin_symmetry = true,
    has_Giomega = false,
    Is_Impurity_model = false,
    has_TRS = false,
    };//a symbol that signifies that we need realtime evolution
    template <class GF, class Configuration, unsigned int>
    friend class PhaseEvaluator;
    typedef Hubbard_Vertex<FPType> Vertex;
    typedef std::complex<FPType> FreeGreensFunctionReturnValueType;///< a typedef of the data type of the free greensfunction
    /**
    This evaluates the value of the free particle Green-function for the two given vertices
    @param v1 the first vertex
    @param v2 the second vertex
    @return the value of the free greensfunction evaluated with the given vertices
    */
    static inline FreeGreensFunctionReturnValueType eval(const Vertex& v1, const Vertex& v2) throw();
    /**
    A function for initializing the Tables that make up the Greensfunction. The length  of the interaction contour and the nr of sites are read from Config Files
     */
    template <class CFG>
    static inline void init(CFG&);
    /**
    This function frees the memory used by the tables for the Greensfunction
     */
    static inline void tidyup();
    /**
    To access the nr of atoms of the chain
    @return the nr of atoms in the chain
     */
    static inline unsigned int getLen() throw()
    {
        return N;
    }
    /**
    Get the Length of the interactioninterval.
    @return the length from the interaction interval specified by t_M and beta
     */
    static inline FPType getContourLen() throw()
    {
        return contourlen;
    }
    /**
    This function takes a Vertex and switches it to the other branch
    */
    static inline Vertex switchBranch(const Vertex& v) throw()
    {
    Vertex retval(v);
    if(v.tau < 2.0 * t_exp)
    {
    retval.tau = 2.0*t_exp - v.tau;
    }
    return retval;
    }
private:
    static FreeGreensFunctionReturnValueType** greal;///< the free particle greensfunction, restricted to values on the real time axis
    static FreeGreensFunctionReturnValueType** grealac;///< the free particle greensfunction, restricted to values on the real time axis
    static FreeGreensFunctionReturnValueType*** gmixed;///< the free particle greensfunction, where it can be fully complex
    static unsigned int slices;///< the number of timeslices. this means the resolution of the free greensfunction in contour space
    static unsigned int realslices;///< the number of timeslices. this means the resolution of the free greensfunction on the realtime axis
    /**
    The lengths of the timeslices on the contour
     */
    static FPType realslice;
    static FPType mixedslicere;
    static FPType mixedsliceim;
    /**
    the inverse temperature beta
     */
    static FPType beta;
    static FPType t_exp;///< the time up to which the expansion is done
    static FPType contourlen;///< this is the length of the contour on which we simulate
    /**
    the length of the chain
     */
    static unsigned int N;
    /**
    The Dispersionrelation -2t cos(k) -my of the free Hubbard model. Note that the chemical potential is included here.
    @param k the momentum
    @param t the kinetic energy
    @param my the chemical potential
    @return the energy
     */
    static inline FPType disp(FPType k, FPType beta, FPType my) throw();
    /**
    The free particle Greens function in k-space
    @param s the Difference of the time - values. it can be arbitrarily complex
     */
    static inline FreeGreensFunctionReturnValueType GreenNullKSpace(FPType, FreeGreensFunctionReturnValueType s, FPType, FPType) throw();
    static inline FreeGreensFunctionReturnValueType GreenNullKSpaceAC(FPType, FreeGreensFunctionReturnValueType s, FPType, FPType) throw();
    /**
    The gamma Function that maps contour times to complex numbers
    @param s the contourtime
    @return the time whether it is real or imaginary
    */
    static inline FreeGreensFunctionReturnValueType gamma(FPType s) throw();
};

template<typename FPType>
typename FreeRealTimeHubbardGreensFunction<FPType>::FreeGreensFunctionReturnValueType** FreeRealTimeHubbardGreensFunction<FPType>::greal = NULL;

template<typename FPType>
typename FreeRealTimeHubbardGreensFunction<FPType>::FreeGreensFunctionReturnValueType** FreeRealTimeHubbardGreensFunction<FPType>::grealac = NULL;

template<typename FPType>
typename FreeRealTimeHubbardGreensFunction<FPType>::FreeGreensFunctionReturnValueType*** FreeRealTimeHubbardGreensFunction<FPType>::gmixed = NULL;

template<typename FPType>
FPType FreeRealTimeHubbardGreensFunction<FPType>::realslice;

template<typename FPType>
FPType FreeRealTimeHubbardGreensFunction<FPType>::mixedslicere;

template<typename FPType>
FPType FreeRealTimeHubbardGreensFunction<FPType>::mixedsliceim;

template<typename FPType>
FPType FreeRealTimeHubbardGreensFunction<FPType>::beta;

template<typename FPType>
FPType FreeRealTimeHubbardGreensFunction<FPType>::t_exp;

template<typename FPType>
FPType FreeRealTimeHubbardGreensFunction<FPType>::contourlen;

template<typename FPType>
unsigned int FreeRealTimeHubbardGreensFunction<FPType>::N;

template<typename FPType>
unsigned int FreeRealTimeHubbardGreensFunction<FPType>::slices;

template<typename FPType>
unsigned int FreeRealTimeHubbardGreensFunction<FPType>::realslices;

template<typename FPType>
typename FreeRealTimeHubbardGreensFunction<FPType>::FreeGreensFunctionReturnValueType FreeRealTimeHubbardGreensFunction<FPType>::gamma(FPType s) throw()
{
    if (s < t_exp) //forward branch
        return s;
    if (s < 2.0 * t_exp) //backward branch
        return 2.0 * t_exp - s;
    return FreeGreensFunctionReturnValueType(0.0, -(s - 2.0 * t_exp));//imaginary branch
}

template <typename FPType>
FPType FreeRealTimeHubbardGreensFunction<FPType>::disp(FPType k, FPType t, FPType my) throw()//some system function seems to be also called epsilon
{
    return -2.0 * t * std::cos(k) - my;
}

template <typename FPType>
void FreeRealTimeHubbardGreensFunction<FPType>::tidyup()
{
    for (unsigned int j = 0; j < realslices+1; ++j)
    {
        delete [] greal[j];
        delete [] grealac[j];
    }
    delete [] greal;
    delete [] grealac;
    unsigned int slicesp = slices + 1;
    for (unsigned int re = 0; re < slicesp; ++re)
    {
        for (unsigned int im = 0; im < slicesp; ++im)
        {
            delete [] gmixed[re][im];
        }
        delete [] gmixed[re];
    }
    delete [] gmixed;
}

template <typename FPType>
typename FreeRealTimeHubbardGreensFunction<FPType>::FreeGreensFunctionReturnValueType FreeRealTimeHubbardGreensFunction<FPType>::GreenNullKSpace(FPType k, std::complex<FPType> z, FPType t, FPType my) throw()
{
    FPType ek = disp(k, t, my);
    return std::exp(std::complex<FPType>(0.0, ek * z.real())) * std::exp(-ek*z.imag()) * fermi(beta * ek);
}

template <typename FPType>
typename FreeRealTimeHubbardGreensFunction<FPType>::FreeGreensFunctionReturnValueType FreeRealTimeHubbardGreensFunction<FPType>::GreenNullKSpaceAC(FPType k, std::complex<FPType> z, FPType t, FPType my) throw()
{
    FPType ek = disp(k, t, my);
    return -std::exp(std::complex<FPType>(0.0, ek * z.real())) * (std::exp(ek*(beta - z.imag()))/(1.0 + std::exp(ek*beta)));
}

template <typename FPType>
template <class CFG>
void FreeRealTimeHubbardGreensFunction<FPType>::init(CFG& curparams)
{
    //lattice sites
    N = curparams.N;
    beta = curparams.beta;
    t_exp = curparams.t_exp;
    FPType t = curparams.t;
    FPType mu = curparams.mu;
    contourlen = beta + 2.0 * t_exp;
    realslices = 5000;//Number of timeslices for the real axis
    slices = 3000;//Number of timeslices for slicing up the complex plane
    const unsigned int slicesp = slices + 1;
    realslice = 2.0 * t_exp / static_cast<FPType>(realslices);//the length of one timeslice on the real axis part
    mixedslicere = t_exp / static_cast<FPType>(slices);//the length of the real part of one slice in the complex part of the plane
    mixedsliceim = beta / static_cast<FPType>(slices);//the length of the imaginary part of one slice in the complex plane
    greal = new FreeGreensFunctionReturnValueType*[realslices/*+1*/];//the array for the Greensfunction evaluated on the purely real part
    grealac = new FreeGreensFunctionReturnValueType*[realslices/*+1*/];//the array for the Greensfunction evaluated on the purely real anti causal part
    gmixed = new FreeGreensFunctionReturnValueType**[slicesp];//the array for Greensfunction evaluated for fully complex arguments
    if (N == 1)//take special care of the 1 site hubbard model
    {
        //set up the table for the real time part
        for (unsigned int j = 0; j < (realslices/*+1*/); ++j)//creates a table with the values from [0, 2*t_exp]
        {
            greal[j] = new FreeGreensFunctionReturnValueType[1];
            grealac[j] = new FreeGreensFunctionReturnValueType[1];
            greal[j][0] = GreenNullKSpace(0, j * realslice, t, mu);
            grealac[j][0] = GreenNullKSpaceAC(0, j * realslice, t, mu);
        }
        //set up the table with the mixed parts
        for (unsigned int i_re = 0; i_re < slicesp; ++i_re)//tabulates from [0, t_exp]
        {
            gmixed[i_re] = new FreeGreensFunctionReturnValueType*[slicesp];
            for (unsigned int i_im = 0; i_im < slicesp; ++i_im)//tabulates from [0, beta]
            {
                gmixed[i_re][i_im] = new FreeGreensFunctionReturnValueType[1];
                gmixed[i_re][i_im][0] = GreenNullKSpace(0, std::complex<FPType>(i_re * mixedslicere, - static_cast<FPType>(i_im) * mixedsliceim), t, mu);
            }
        }
        return;
    }
// a static omp for loop is sufficient here
#pragma omp parallel for
    for (int j = 0; j < static_cast<int>(realslices/*+1*/); ++j)//creates a table with the values from [0, 2*t_exp]
    {
        greal[j] = new FreeGreensFunctionReturnValueType[N];
        grealac[j] = new FreeGreensFunctionReturnValueType[N];
        for (int p = 0; p < static_cast<int>(N); ++p)//for every realspacepoint
        {
            std::complex<FPType> tempsum = 0;
            std::complex<FPType> tempsumac = 0;
            for (FPType k = 0; k < static_cast<FPType>(N); ++k)//sum over all k-space values
            {
                tempsum += exp(std::complex<FPType>(0.0, 2.0 * M_PI/N * k * p) ) * GreenNullKSpace(k * 2.0 * M_PI /  N, j * realslice, t, mu);//R_p = p, because the lattice constant is assumed to be equal to one.
                tempsumac += exp(std::complex<FPType>(0.0, 2.0 * M_PI/N * k * p) ) * GreenNullKSpaceAC(k * 2.0 * M_PI /  N, j * realslice, t, mu);
            }
            greal[j][p] = tempsum / static_cast<FPType>(N);// the Fourier transform has 1/N as normalization
            grealac[j][p] = tempsumac / static_cast<FPType>(N);// the Fourier transform has 1/N as normalization
        }
    }
    //set up the table with the mixed parts. A static omp for loop is sufficient here
#pragma omp parallel for
    for (int i_re = 0; i_re < static_cast<int>(slicesp); ++i_re)//tabulates from [0, t_exp]
    {
        gmixed[i_re] = new FreeGreensFunctionReturnValueType*[slicesp];
        for (unsigned int i_im = 0; i_im < slicesp; ++i_im)//tabulates from [0, beta]
        {
            gmixed[i_re][i_im] = new FreeGreensFunctionReturnValueType[N];
            for (int p = 0; p < static_cast<int>(N); ++p)
            {
                std::complex<FPType> tempsum = 0;
                std::complex<FPType> tempsumac = 0;
                for (FPType k = 0; k < static_cast<FPType>(N); ++k)//sum over all k-space values
                {
                    tempsum += exp(std::complex<FPType>(0.0, 2.0 * M_PI/N * k * p) ) * GreenNullKSpace(k * 2.0 * M_PI /  N, std::complex<FPType>(i_re * mixedslicere, - static_cast<FPType>(i_im) * mixedsliceim), t, mu);//R_p = p, because the lattice constant is assumed to be equal to one.
                }
                gmixed[i_re][i_im][p] = tempsum / static_cast<FPType>(N);// the Fourier transform has 1/N as normalization
            }
        }
    }
    return;
}

template<typename FPType>
typename FreeRealTimeHubbardGreensFunction<FPType>::FreeGreensFunctionReturnValueType FreeRealTimeHubbardGreensFunction<FPType>::eval(const Vertex& v1, const Vertex& v2) throw()
{
    const FPType tiny = 0.00000001;
    //determine the Differences between the two
    std::complex<FPType> delta_gamma = gamma(v1.tau) - gamma(v2.tau);//the difference between two different complex arguments
    int delta = v1.site - v2.site;//the difference in sites
    typename FreeRealTimeHubbardGreensFunction<FPType>::FreeGreensFunctionReturnValueType retval;
    if ( std::fabs(v1.tau - v2.tau) < tiny)//if the difference in CONTOUR-TIME is zero we just return the particle number
    {
        //Take care of negative values
        if (delta < 0) delta = N + delta;//because delta is in this case already a negative number
        //return only the particle number
        return greal[0][delta];//we take the particle number from the time ordered realtime data
    }
    const FPType imaginaryPart(imag(delta_gamma));
    if (v1.tau >= 2.0 * t_exp)
    {
        //vertex 1 is on the imaginary axis
        FPType re_idx, im_idx;
        if (v2.tau >= 2.0 *t_exp)
        {
            //both vertices are on the imaginary axis
            if (delta < 0) delta += N;
            if (imaginaryPart < 0)
            {//v1.s > v2.s --> right ordering
                //the conventional greensfunction
                //+-
                FPType im_rem(std::modf((-imaginaryPart)/mixedsliceim, &im_idx));
                retval = gmixed[0][lround(im_idx)][delta];
            }
            else
            {//v1.s < v2.s --> wrong ordering
                //we add a beta for the ordering
                //++
                std::modf(-(imaginaryPart-beta)/mixedsliceim, &im_idx);
                retval = -gmixed[0][lround(im_idx)][delta];
            }
        }
        else
        {
            //we know we have to return the "lesser" part of the greensfunction, because v2 is on the real axis => v1.s > v2.s'
            //indepent of wether v2 is on the forward or on the backward contour we always have to return the lesser part of the greens function and the realpart of delta_gamma is always negative. The imaginary part is also always negative. The tables for the lesser greensfunction are made for negative imaginary-parts
            FPType im_rem(std::modf((-imaginaryPart)/mixedsliceim, &im_idx));
            delta = -delta;
            if (delta < 0) delta += N;
            FPType re_rem(std::modf(-real(delta_gamma)/mixedslicere, &re_idx));
            return conj(gmixed[lround(re_idx)][lround(im_idx)][delta]);//by a conj  we're fixing up the fact, that we don't have tables for positive real part
        }
    }
    else
    {
        if (v2.tau >= 2.0 * t_exp)//v1 in real part v2 in imaginary part => v2.s > v1.s
        {
            //we return the "greater" part of the greensfunction
            FPType re_idx, im_idx;
            FPType re_rem(std::modf(real(delta_gamma)/mixedslicere, &re_idx));
            FPType im_rem(std::modf((beta-imaginaryPart)/mixedsliceim, &im_idx));
            //indepent of wether v1 is on the forward or on the backward contour we always have to return the greater part of the greens function and the realpart of delta_gamma is always positive. The imaginary part is also always positive. the tables for the greater greensfunction are made for positive imaginary-parts
            if (delta < 0) delta += N;
            return -gmixed[lround(re_idx)][lround(im_idx)][delta];
        }
        else
        {
            //both vertices are in the real part of the contour
            const FPType realpart = real(delta_gamma);
            if (v1.tau < t_exp)
            {
                if (v2.tau < t_exp)
                {
                    //both on the forward contour, enforce contourordering
                    FPType re_idx;
                    if (v1.tau < v2.tau)
                    {
                        //bad order, that means we return the "greater" greensfunctions
                        FPType re_rem(std::modf(-realpart/realslice, &re_idx));//the time is a negative value
                        delta = -delta;
                        if (delta < 0) delta += N;
                        retval = conj(grealac[lround(re_idx)][delta]);
                    }
                    else
                    {
                        //good order, this means we return the "lesser" greensfunction
                        FPType re_rem(std::modf(realpart/realslice, &re_idx));//real part is greater than zero
                        if (delta < 0) delta += N;
                        retval = greal[lround(re_idx)][delta];
                    }
                }
                else
                {
                    //v1 on forward v2 on backward => v2.tau > v1.tau => enforce contourordering
                    FPType re_idx;
                    FPType re_rem(std::modf(realpart/realslice, &re_idx));
		    if(realpart >= 0)
		    {
		    if (delta < 0) delta += N;
                    retval = grealac[lround(re_idx)][delta];//the real part is always lesser than zero. we have to compensate this with a complex conjugation
		    }
		    else
		    {
                    delta = -delta;
                    if (delta < 0) delta += N;
                    retval = conj(grealac[lround(-re_idx)][delta]);//the real part is always lesser than zero. we have to compensate this with a complex conjugation
		    }
                }
            }
            else
            {
                if (v2.tau < t_exp)
                {
                    //v1 on backward v2 on forward contour => v1.tau > v2.tau => return lesser greensfunction
                    FPType re_idx;//the real part has an arbitrary sign
		    FPType re_rem(std::modf(realpart/realslice, &re_idx));
		    if(realpart >= 0.0)
		    {
                    if (delta < 0) delta += N;
                    return greal[lround(re_idx)][delta];
		    }
		    else
		    {
		    delta = -delta;
                    if (delta < 0) delta += N;
                    return conj(greal[lround(-re_idx)][delta]);
		    }
                }
                else
                {
                    //v1, v2 on backward contour. check contourordering
                    FPType re_idx;
                    FPType re_rem(std::modf(realpart/realslice, &re_idx));
                    if (v1.tau < v2.tau)
                    {
                        //wrong Contour-Order, thus we return the "greater" Greensfunction, aber der realteil von deltagamma ist positiv
                        if (delta < 0) delta += N;
                        return grealac[lround(re_idx)][delta];
                    }
                    else
                    {
                        //v1 is after v2
                        //right contour-order, thus we return the "lesser" greensfunction, aber der realteil von delta_gamma ist negativ
                        delta = -delta;
                        if (delta < 0) delta += N;
                        retval = conj(greal[lround(-re_idx)][delta]);
                    }
                }
            }
        }
    }
    return retval;
}

/**
An implementation of a free-particle greens-function for the Cold-Atom 
simulations that doesn't tabulate the realtime branches.
I assume that the lattice constant is a = 1
*/
template<typename FPType = double>
class FreeRealTimeHubbardGreensFunctionLessMem
{
public:
    enum {timeevolution = 1,//a symbol that signifies that we need realtime evolution
    has_Spin_symmetry = 1,//a symbol to denote that this greensfunction has the symmetry in its spin-sectors.
    has_Giomega = false,
    Is_Impurity_model = false,
    has_TRS = false,
    };
    template <class GF, class Configuration, unsigned int>
    friend class PhaseEvaluator;
    typedef Hubbard_Vertex<FPType> Vertex;
    typedef std::complex<FPType> FreeGreensFunctionReturnValueType;///< a typedef of the data type of the free greensfunction
    /**
    This evaluates the value of the free particle Green-function for the two given vertices
    @param v1 the first vertex
    @param v2 the second vertex
    @return the value of the free greensfunction evaluated with the given vertices
    */
    static inline FreeGreensFunctionReturnValueType eval(const Vertex& v1, const Vertex& v2) throw();
    /**
    A function for initializing the Tables that make up the Greensfunction. The length  of the interaction contour and the nr of sites are read from Config Files
     */
    template <class CFG>
    static inline void init(CFG&);
    /**
    This function frees the memory used by the tables for the Greensfunction
     */
    static inline void tidyup();
    /**
    To access the nr of atoms of the chain
    @return the nr of atoms in the chain
     */
    static inline unsigned int getLen() throw()
    {
        return N;
    }
    /**
    This function takes a Vertex and switches it to the other branch
    */
    static inline Vertex switchBranch(const Vertex& v) throw()
    {
    Vertex retval(v);
    if(v.tau < 2.0 * t_exp)
    {
    retval.tau = 2.0*t_exp - v.tau;
    }
    return retval;
    }
private:
    static FreeGreensFunctionReturnValueType** gimag;///< the free particle greensfunction, restricted to values on the imaginary axis
    static unsigned int slices;///< the number of timeslices. this means the resolution of the free greensfunction on the imaginary axis
    static FPType imagslice;///< the length of a timeslice on the imaginary axis
    static FPType beta;///< the inverse temperature beta
    static FPType t_exp;///< the time up to which the expansion is done
    static unsigned int N;///< the length of the chain
    /**
    The Dispersionrelation -2t cos(k) -my of the free Hubbard model. Note that the chemical potential is included here.
    The lattice constant is assumed to be one.
    @param k the momentum
    @param t the kinetic energy
    @param my the chemical potential
    @return the energy
     */
    static inline FPType disp(FPType k, FPType beta, FPType my) throw();
    /**
    The free particle Greens function in k-space
    @param s the Difference of the time - values. it can be arbitrarily complex
    */
    static inline FreeGreensFunctionReturnValueType GreenNullKSpace(int, FreeGreensFunctionReturnValueType s) throw();
    /**
    The gamma function that maps contour times to complex numbers
    @param s the contourtime
    @return the time whether it is a real time or an imaginary time
    */
    static inline FreeGreensFunctionReturnValueType gamma(FPType s) throw();
    static inline FreeGreensFunctionReturnValueType GZero(int d, FreeGreensFunctionReturnValueType z) throw();
    static inline FreeGreensFunctionReturnValueType GZero(int d, FPType re) throw();
    static inline FreeGreensFunctionReturnValueType GZeroAC(int d, FPType re) throw();
    static FPType** imagtimeEV;
    static FPType realslice;
    static FreeGreensFunctionReturnValueType** exponentials;///< this contains the exponential prefactors of the Fouriertransform
    static FreeGreensFunctionReturnValueType** realtimeEVf;///< this is supposed to contain the realtime evolution multplied by the fermi function
    static FreeGreensFunctionReturnValueType** realtimeEVACf;///< this is supposed to contain the realtime evolution multplied by the 1 - fermi function
};

template <typename FPType>
typename FreeRealTimeHubbardGreensFunctionLessMem<FPType>::FreeGreensFunctionReturnValueType FreeRealTimeHubbardGreensFunctionLessMem<FPType>::GZero(int d, std::complex<FPType> z) throw()
{
    FreeGreensFunctionReturnValueType tempsum = 0;
    if (d < 0) d += N;
    for (unsigned int k = 0; k < N; ++k)//sum over all k-space values
    {
        tempsum += exponentials[k][d] * GreenNullKSpace(k, z);//R_p = p, because the lattice constant is assumed to be equal to one.
    }
//    if(imag(tempsum) < 0.1* std::numeric_limits<FPType>::epsilon()) tempsum.imag() = 0;
    return tempsum / static_cast<FPType>(N);// the Fourier transform has 1/N as normalization
}

template <typename FPType>
typename FreeRealTimeHubbardGreensFunctionLessMem<FPType>::FreeGreensFunctionReturnValueType FreeRealTimeHubbardGreensFunctionLessMem<FPType>::GZero(int d, const FPType re) throw()
{
    FreeGreensFunctionReturnValueType tempsum = 0;
    FPType re_idx;
    FPType re_rem(std::modf(re/realslice, &re_idx));
    bool doconj = false;
    if(re_idx < 0)
    {
    d = -d;
    re_idx = -re_idx;
    doconj = true;
    }
    if (d < 0) d += N;
    for (unsigned int k = 0; k < N; ++k)//sum over all k-space values
    {
        tempsum += exponentials[k][d] * realtimeEVf[lround(re_idx)][k];//R_p = p, because the lattice constant is assumed to be equal to one.
    }
    if (doconj) tempsum = conj(tempsum);
//    if(imag(tempsum) < 0.1* std::numeric_limits<FPType>::epsilon()) tempsum.imag() = 0;
    return tempsum / static_cast<FPType>(N);// the Fourier transform has 1/N as normalization
}

template <typename FPType>
typename FreeRealTimeHubbardGreensFunctionLessMem<FPType>::FreeGreensFunctionReturnValueType FreeRealTimeHubbardGreensFunctionLessMem<FPType>::GZeroAC(int d, const FPType re) throw()
{
    FreeGreensFunctionReturnValueType tempsum = 0;
    FPType re_idx;
    FPType re_rem(std::modf(-re/realslice, &re_idx));
    bool doconj = false;
    if(re_idx < 0)
    {
    d = -d;
    re_idx = -re_idx;
    doconj = true;
    }
    if (d < 0) d += N;
    for (unsigned int k = 0; k < N; ++k)//sum over all k-space values
    {
       tempsum += -exponentials[k][d] * realtimeEVACf[lround(re_idx)][k];//R_p = p, because the lattice constant is assumed to be equal to one.
    }
//    if(imag(tempsum) < 0.1* std::numeric_limits<FPType>::epsilon()) tempsum.imag() = 0;
    if (doconj) tempsum = conj(tempsum);
    return tempsum / static_cast<FPType>(N);// the Fourier transform has 1/N as normalization
}

template<typename FPType>
typename FreeRealTimeHubbardGreensFunctionLessMem<FPType>::FreeGreensFunctionReturnValueType** FreeRealTimeHubbardGreensFunctionLessMem<FPType>::gimag = NULL;

template<typename FPType>
typename FreeRealTimeHubbardGreensFunctionLessMem<FPType>::FreeGreensFunctionReturnValueType** FreeRealTimeHubbardGreensFunctionLessMem<FPType>::exponentials = NULL;

template<typename FPType>
typename FreeRealTimeHubbardGreensFunctionLessMem<FPType>::FreeGreensFunctionReturnValueType** FreeRealTimeHubbardGreensFunctionLessMem<FPType>::realtimeEVf = NULL;

template<typename FPType>
typename FreeRealTimeHubbardGreensFunctionLessMem<FPType>::FreeGreensFunctionReturnValueType** FreeRealTimeHubbardGreensFunctionLessMem<FPType>::realtimeEVACf = NULL;

template<typename FPType>
FPType FreeRealTimeHubbardGreensFunctionLessMem<FPType>::imagslice;

template<typename FPType>
FPType FreeRealTimeHubbardGreensFunctionLessMem<FPType>::realslice;

template<typename FPType>
FPType** FreeRealTimeHubbardGreensFunctionLessMem<FPType>::imagtimeEV = NULL;

template<typename FPType>
FPType FreeRealTimeHubbardGreensFunctionLessMem<FPType>::beta;

template<typename FPType>
FPType FreeRealTimeHubbardGreensFunctionLessMem<FPType>::t_exp;

template<typename FPType>
unsigned int FreeRealTimeHubbardGreensFunctionLessMem<FPType>::N;

template<typename FPType>
unsigned int FreeRealTimeHubbardGreensFunctionLessMem<FPType>::slices;

template<typename FPType>
typename FreeRealTimeHubbardGreensFunctionLessMem<FPType>::FreeGreensFunctionReturnValueType FreeRealTimeHubbardGreensFunctionLessMem<FPType>::gamma(FPType s) throw()
{
    if (s < t_exp) //forward branch
        return s;
    if (s < 2.0 * t_exp) //backward branch
        return 2.0 * t_exp - s;
    return FreeGreensFunctionReturnValueType(0.0, -(s - 2.0 * t_exp));//imaginary branch
}

template <typename FPType>
FPType FreeRealTimeHubbardGreensFunctionLessMem<FPType>::disp(FPType k, FPType t, FPType my) throw()//some system function seems to be also called epsilon
{
    return -2.0 * t * std::cos(k) - my;
}

template <typename FPType>
void FreeRealTimeHubbardGreensFunctionLessMem<FPType>::tidyup()
{
    for (unsigned int j = 0; j < slices; ++j)
    {
        delete [] gimag[j];
        delete [] imagtimeEV[j];
        delete [] realtimeEVf[j];
        delete [] realtimeEVACf[j];
    }
    delete [] imagtimeEV[slices];
    delete [] realtimeEVf[slices];
    delete [] realtimeEVACf[slices];
    for(unsigned int r = 0; r < N; ++r)
        delete [] exponentials[r];
    delete [] exponentials;
    delete [] gimag;
    delete [] imagtimeEV;
    delete [] realtimeEVf;
    delete [] realtimeEVACf;
}

template <typename FPType>
typename FreeRealTimeHubbardGreensFunctionLessMem<FPType>::FreeGreensFunctionReturnValueType FreeRealTimeHubbardGreensFunctionLessMem<FPType>::GreenNullKSpace(int n, std::complex<FPType> z) throw()
{
    FPType im_idx;
    FPType im_rem(std::modf(-z.imag()/imagslice, &im_idx));
    FPType re_idx;
    FPType re_rem(std::modf(-z.real()/realslice, &re_idx));
    return realtimeEVf[lround(re_idx)][n] * imagtimeEV[lround(im_idx)][n];
}

template <typename FPType>
template <class CFG>
void FreeRealTimeHubbardGreensFunctionLessMem<FPType>::init(CFG& curparams)
{
    //lattice sites
    N = curparams.N;
    beta = curparams.beta;
    t_exp = curparams.t_exp;
    FPType t = curparams.t;
    FPType mu = curparams.mu;
    slices = 10000;//Number of timeslices for the imaginary axis
    const unsigned int slicesp = slices + 1;
    imagslice =  beta/ static_cast<FPType>(slices);//the length of one timeslice on the imaginary axis part
    realslice =  t_exp/ static_cast<FPType>(slices);//the length of one timeslice on the real axis part
    gimag = new FreeGreensFunctionReturnValueType*[slices];//the array for the Greensfunction evaluated on the purely imaginary part
    FPType* disparr = new FPType[N];
    FPType* fermiarr = new FPType[N];
    FPType* fermiarrAC = new FPType[N];
    exponentials = new FreeGreensFunctionReturnValueType*[N];
    imagtimeEV = new FPType*[slicesp];
    realtimeEVf = new FreeGreensFunctionReturnValueType*[slicesp];
    realtimeEVACf = new FreeGreensFunctionReturnValueType*[slicesp];
    if (N == 1)//take special care of the 1 site hubbard model
    {
        disparr[0] = disp(0,t,mu);
        fermiarr[0] = fermi(beta*disparr[0]);
        fermiarrAC[0] = std::exp(disparr[0]*beta)*fermiarr[0];
        exponentials[0] = new FreeGreensFunctionReturnValueType[1];
        exponentials[0][0] = 1.0;
        //set up the table for the imaginary time part
        for (unsigned int j = 0; j < slices; ++j)//creates a table with the values from [0, beta]
        {
            gimag[j] = new FreeGreensFunctionReturnValueType[1];
            imagtimeEV[j] = new FPType[1];
            realtimeEVf[j] = new FreeGreensFunctionReturnValueType[1];
            realtimeEVACf[j] = new FreeGreensFunctionReturnValueType[1];
        }
        imagtimeEV[slices] = new FPType[1];
        realtimeEVf[slices] = new FreeGreensFunctionReturnValueType[1];
        realtimeEVACf[slices] = new FreeGreensFunctionReturnValueType[1];
	delete [] disparr;
	delete [] fermiarr;
	delete [] fermiarrAC;
        return;
    }
    //fill up the array with the dispersion-relation and anything else
    for(unsigned int k = 0; k < N; ++k)
    {
        FPType kf = static_cast<FPType>(k);
        disparr[k] = disp(kf * 2.0 * M_PI / N, t,mu);
        fermiarr[k] = fermi(beta*disparr[k]);
        fermiarrAC[k] = std::exp(disparr[k]*beta)*fermiarr[k];
        exponentials[k] = new FreeGreensFunctionReturnValueType[N];
        for(unsigned int r = 0; r < N; ++r)
           exponentials[k][r] = std::exp(std::complex<FPType>(0.0, kf * 2.0 * M_PI / N * r) );
    }
    //a static omp for loop is sufficient here
#pragma omp parallel for
    for (int j = 0; j < static_cast<int>(slices); ++j)//creates a table with the values from [0, beta]
    {
        gimag[j] = new FreeGreensFunctionReturnValueType[N];
        imagtimeEV[j] = new FPType[N];
        realtimeEVf[j] = new FreeGreensFunctionReturnValueType[N];
        realtimeEVACf[j] = new FreeGreensFunctionReturnValueType[N];
        for (int p = 0; p < static_cast<int>(N); ++p)//for every realspacepoint
        {
            imagtimeEV[j][p] = std::exp(disparr[p]*static_cast<FPType>(j)*imagslice);
            realtimeEVf[j][p] = std::exp(std::complex<FPType>(0.0, disparr[p]*static_cast<FPType>(j)*realslice))*fermiarr[p];
            realtimeEVACf[j][p] = std::exp(std::complex<FPType>(0.0, -disparr[p]*static_cast<FPType>(j)*realslice))*fermiarrAC[p];
        }
    }
#pragma omp parallel for
    for (int j = 0; j < static_cast<int>(slices); ++j)//fills the table with the values from [0, beta]
    {
        for (int p = 0; p < static_cast<int>(N); ++p)//for every realspacepoint
        {
            std::complex<FPType> tempsum = 0;
            for (int k = 0; k < static_cast<int>(N); ++k)//sum over all k-space values
            {
                tempsum += exponentials[k][p] * GreenNullKSpace(k, std::complex<FPType>(0.0, -j * imagslice) );//R_p = p, because the lattice constant is assumed to be equal to one.
            }
            gimag[j][p] = tempsum / static_cast<FPType>(N);// the Fourier transform has 1/N as normalization
        }
    }
    imagtimeEV[slices] = new FPType[N];
    realtimeEVf[slices] = new FreeGreensFunctionReturnValueType[N];
    realtimeEVACf[slices] = new FreeGreensFunctionReturnValueType[N];
    for (int p = 0; p < static_cast<int>(N); ++p)//for every realspacepoint
    {
         imagtimeEV[slices][p] = std::exp(disparr[p]*static_cast<FPType>(slices)*imagslice);
         realtimeEVf[slices][p] = std::exp(std::complex<FPType>(0.0, disparr[p]*static_cast<FPType>(slices)*realslice))*fermiarr[p];
         realtimeEVACf[slices][p] = std::exp(std::complex<FPType>(0.0, -disparr[p]*static_cast<FPType>(slices)*realslice))*fermiarrAC[p];
    }
    delete [] disparr;
    delete [] fermiarr;
    delete [] fermiarrAC;
    return;
}

template<typename FPType>
typename FreeRealTimeHubbardGreensFunctionLessMem<FPType>::FreeGreensFunctionReturnValueType FreeRealTimeHubbardGreensFunctionLessMem<FPType>::eval(const Vertex& v1, const Vertex& v2) throw()
{
    const FPType tiny = 0.00000001;
    //determine the Differences between the two
    std::complex<FPType> delta_gamma = gamma(v1.tau) - gamma(v2.tau);//the difference between two different complex arguments
    int delta = v1.site - v2.site;//the difference in sites
    typename FreeRealTimeHubbardGreensFunction<FPType>::FreeGreensFunctionReturnValueType retval;
    if ( std::fabs(v1.tau - v2.tau) < tiny)//if the difference in CONTOUR-TIME is zero we just return the particle density
    {
        //Take care of negative values
        if (delta < 0) delta = N + delta;//because delta is in this case already a negative number
        //return only the particle density
        return gimag[0][delta];//we take the particle number from the imaginary time data
    }
    const FPType imaginaryPart(imag(delta_gamma));
    if (v1.tau >= 2.0 * t_exp)
    {
        //vertex 1 is on the imaginary axis
        if (v2.tau >= 2.0 *t_exp)
        {
            FPType im_idx;
            //both vertices are on the imaginary axis
            if (delta < 0) delta += N;
            if (imaginaryPart < 0)
            {
                //the conventional greensfunction
                //+-
                FPType im_rem(std::modf((-imaginaryPart)/imagslice, &im_idx));
		long int im_idx_int = lround(im_idx);
                retval = gimag[im_idx_int][delta];
            }
            else
            {
                //we add a beta for the ordering
                //++
                FPType im_rem(std::modf((beta-imaginaryPart)/imagslice, &im_idx));
		long int im_idx_int = lround(im_idx);
                retval = -gimag[im_idx_int][delta];
            }
        }
        else
        {
            //we know we have to return the "lesser" part of the greensfunction, because v2 is on the real axis => v1.s > v2.s'
            //indepent of wether v2 is on the forward or on the backward contour we always have to return the lesser part of the greens function and the realpart of delta_gamma is always negative. The imaginary part is also always negative. The tables for the lesser greensfunction are made for negative imaginary-parts
	    delta = -delta;
            retval = conj(GZero(delta, delta_gamma)); //Re < 0, Im < 0
        }
    }
    else
    {
        if (v2.tau >= 2.0 * t_exp)//v1 in real part v2 in imaginary part => v2.s > v1.s
        {
            //we return the "greater" part of the greensfunction
            //indepent of wether v1 is on the forward or on the backward contour we always have to return the greater part of the greens function.
            return -GZero(delta, std::complex<FPType>(-real(delta_gamma),imaginaryPart-beta));//we still need the minus sign for the time-ordering!
        }
        else
        {
            //both vertices are in the real part of the contour
            const FPType realpart = real(delta_gamma);
            if (v1.tau < t_exp)
            {
                if (v2.tau < t_exp)
                {
                    //both on the forward contour, enforce contourordering
                    if (v1.tau < v2.tau)
                    {
                        //bad order, that means we return the "greater" greensfunctions
                        retval = GZeroAC(delta, realpart);//Im = 0, Re < 0
                    }
                    else
                    {
                        //good order, this means we return the "lesser" greensfunction
                        retval = GZero(delta, realpart);//Re > 0, Im = 0
                    }
                }
                else
                {
                    //v1 on forward v2 on backward => v2.tau > v1.tau => enforce contourordering
                    retval = GZeroAC(delta, realpart);
                }
            }
            else
            {
                if (v2.tau < t_exp)
                {
                    //v1 on backward v2 on forward contour => v1.tau > v2.tau => return lesser greensfunction
                    retval = GZero(delta, realpart);
                }
                else
                {
                    //v1, v2 on backward contour. check contourordering
                    if (v1.tau < v2.tau)
                    {
                        //wrong Contour-Order, thus we return the "greater" Greensfunction, aber der realteil von deltagamma ist positiv
                        retval = conj(GZeroAC(delta, -realpart));// Re > 0
                    }
                    else
                    {
                        //v1 is after v2
                        //right contour-order, thus we return the "lesser" greensfunction, aber der realteil von delta_gamma ist negativ
                        retval = GZero(delta, realpart); // Re < 0
                    }
                }
            }
        }
    }
    return retval;
}
//#endif


#endif
